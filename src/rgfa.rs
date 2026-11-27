//! Incremental `minigraph`-compatible rGFA (`SN:Z:`, `SO:i:`, `SR:i:`) SV-graph builder and parser.

use crate::mdl::{AlignConfig, AlignedRun, align_contig_pair_runs};
use crate::netting::net_one_to_one_runs_scored;
use crate::seeding::{ContigIndex, prepare_contig_index_owned, select_top_k_pairs};
use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use std::collections::HashSet;
use std::io::{BufRead, Write};

/// Default upper bound (bp) for a single grafted SV insertion/deletion bubble to prevent
/// translocation or paralogous jumps from grafting multi-Mb bogus segments.
pub const DEFAULT_MAX_SV_BP: usize = 1_000_000;

/// A single rGFA segment with stable linear coordinate provenance (`SN`, `SO`, `SR`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RgfaSegment {
    pub id: String,
    pub seq: Vec<u8>,
    /// Stable origin sequence/contig name (`SN:Z:`).
    pub origin_name: String,
    /// 0-based start offset on `origin_name` (`SO:i:`).
    pub origin_offset: usize,
    /// Graph rank (`SR:i:`): `0` for the initial reference backbone, `>= 1` for grafted assemblies.
    pub rank: u32,
}

/// A directed GFA link `(from, from_fwd) -> (to, to_fwd)`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RgfaLink {
    pub from: String,
    pub from_fwd: bool,
    pub to: String,
    pub to_fwd: bool,
}

/// Incremental `minigraph`-style rGFA structural-variation graph.
#[derive(Debug, Clone, Default)]
pub struct RgfaGraph {
    pub segments: Vec<RgfaSegment>,
    pub links: Vec<RgfaLink>,
    pub min_sv_bp: usize,
    pub max_sv_bp: usize,
    pub skip_softmasked: bool,
    pub next_rank: u32,
    next_seg_num: usize,
    link_set: HashSet<RgfaLink>,
}

impl RgfaGraph {
    /// Initialize a rank-0 linear reference backbone from `(contig_name, seq)` pairs.
    #[must_use]
    pub fn from_reference(ref_contigs: &[(&str, &[u8])], min_sv_bp: usize) -> Self {
        let mut segments = Vec::with_capacity(ref_contigs.len());
        let mut next_seg_num = 1usize;
        for &(name, seq) in ref_contigs {
            if seq.is_empty() {
                continue;
            }
            segments.push(RgfaSegment {
                id: format!("s{next_seg_num}"),
                seq: seq.to_vec(),
                origin_name: name.to_owned(),
                origin_offset: 0,
                rank: 0,
            });
            next_seg_num += 1;
        }
        Self {
            segments,
            links: Vec::new(),
            min_sv_bp: min_sv_bp.max(1),
            max_sv_bp: DEFAULT_MAX_SV_BP,
            skip_softmasked: false,
            next_rank: 1,
            next_seg_num,
            link_set: HashSet::new(),
        }
    }

    /// Configure whether softmasked (lowercase) bases are skipped during seed indexing and alignment.
    #[must_use]
    pub const fn with_respect_softmask(mut self, enable: bool) -> Self {
        self.skip_softmasked = enable;
        self
    }

    fn alloc_seg_id(&mut self) -> String {
        let id = format!("s{}", self.next_seg_num);
        self.next_seg_num += 1;
        id
    }

    fn add_link_dedup(&mut self, from: &str, from_fwd: bool, to: &str, to_fwd: bool) {
        let link = RgfaLink {
            from: from.to_owned(),
            from_fwd,
            to: to.to_owned(),
            to_fwd,
        };
        if self.link_set.insert(link.clone()) {
            self.links.push(link);
        }
    }

    /// Split rank-0 segment covering `(origin_name, split_pos)` strictly interior to a segment.
    fn split_rank0_at(&mut self, origin_name: &str, split_pos: usize) {
        let Some(idx) = self.segments.iter().position(|s| {
            s.rank == 0
                && s.origin_name == origin_name
                && s.origin_offset < split_pos
                && s.origin_offset
                    .checked_add(s.seq.len())
                    .is_some_and(|end| split_pos < end)
        }) else {
            return;
        };
        let old_id = self.segments[idx].id.clone();
        let rel = split_pos - self.segments[idx].origin_offset;
        let tail_seq = self.segments[idx].seq.split_off(rel);
        // Release spare capacity on the left piece so K left-to-right splits stay O(L) RAM
        // instead of leaking O(K * L) bytes (~30+ GB on mega-chromosomes).
        self.segments[idx].seq.shrink_to_fit();
        let right_id = self.alloc_seg_id();
        let tail_node = RgfaSegment {
            id: right_id.clone(),
            seq: tail_seq,
            origin_name: origin_name.to_owned(),
            origin_offset: split_pos,
            rank: 0,
        };
        self.segments.insert(idx + 1, tail_node);

        // Redirect any outgoing `+` links from old_id AND any incoming `-` links to old_id
        // (independent `if`s so self-loops `L old + old -` update both ends).
        for l in &mut self.links {
            let mut changed = false;
            if (l.from == old_id && l.from_fwd) || (l.to == old_id && !l.to_fwd) {
                self.link_set.remove(l);
                if l.from == old_id && l.from_fwd {
                    l.from.clone_from(&right_id);
                }
                if l.to == old_id && !l.to_fwd {
                    l.to.clone_from(&right_id);
                }
                changed = true;
            }
            if changed {
                self.link_set.insert(l.clone());
            }
        }
        self.add_link_dedup(&old_id, true, &right_id, true);
    }

    /// Find the rank-0 segment ID ending at `pos` on `origin_name`.
    fn rank0_ending_at(&self, origin_name: &str, pos: usize) -> Option<String> {
        self.segments
            .iter()
            .find(|s| {
                s.rank == 0
                    && s.origin_name == origin_name
                    && s.origin_offset.checked_add(s.seq.len()) == Some(pos)
            })
            .map(|s| s.id.clone())
    }

    /// Find the rank-0 segment ID starting at `pos` on `origin_name`.
    fn rank0_starting_at(&self, origin_name: &str, pos: usize) -> Option<String> {
        self.segments
            .iter()
            .find(|s| s.rank == 0 && s.origin_name == origin_name && s.origin_offset == pos)
            .map(|s| s.id.clone())
    }

    /// Reconstruct full linear sequence of each rank-0 reference contig.
    fn rank0_contigs(&self) -> Vec<(String, Vec<u8>)> {
        let mut seen = HashSet::new();
        let mut names = Vec::new();
        for s in &self.segments {
            if s.rank == 0 && seen.insert(s.origin_name.clone()) {
                names.push(s.origin_name.clone());
            }
        }
        names
            .into_iter()
            .map(|name| {
                let mut parts: Vec<&RgfaSegment> = self
                    .segments
                    .iter()
                    .filter(|s| s.rank == 0 && s.origin_name == name)
                    .collect();
                parts.sort_by_key(|s| s.origin_offset);
                let seq: Vec<u8> = parts
                    .into_iter()
                    .flat_map(|s| s.seq.iter().copied())
                    .collect();
                (name, seq)
            })
            .collect()
    }

    /// Incrementally graft a multi-contig query assembly under a single rGFA rank (`self.next_rank`),
    /// indexing rank-0 reference contigs once for the entire assembly, netting runs 1-to-1, and
    /// preserving forward rGFA `SN`/`SO` coordinates (using `-` links for reverse-strand SVs).
    #[allow(clippy::too_many_lines)]
    pub fn graft_assembly_contigs(&mut self, query_contigs: &[(&str, &[u8])]) -> usize {
        if query_contigs.is_empty() {
            return 0;
        }
        let skip_soft = self.skip_softmasked;
        let ref_contigs = self.rank0_contigs();
        let ref_indexed: Vec<(String, ContigIndex)> = ref_contigs
            .into_par_iter()
            .filter_map(|(name, seq)| {
                prepare_contig_index_owned(seq, skip_soft)
                    .ok()
                    .map(|idx| (name, idx))
            })
            .collect();
        if ref_indexed.is_empty() {
            self.next_rank = self.next_rank.saturating_add(1);
            return 0;
        }
        let q_indexed: Vec<(String, ContigIndex)> = query_contigs
            .par_iter()
            .filter(|(_, s)| !s.is_empty())
            .filter_map(|&(name, seq)| {
                prepare_contig_index_owned(seq.to_vec(), skip_soft)
                    .ok()
                    .map(|idx| (name.to_owned(), idx))
            })
            .collect();
        if q_indexed.is_empty() {
            self.next_rank = self.next_rank.saturating_add(1);
            return 0;
        }

        let pair_jobs: Vec<(usize, usize)> = if ref_indexed.len() > 1 && q_indexed.len() > 1 {
            let mask = select_top_k_pairs(&ref_indexed, &q_indexed, 1, 8);
            let mut jobs = Vec::new();
            for (qi, _) in q_indexed.iter().enumerate() {
                for (ri, row) in mask.iter().enumerate() {
                    if row[qi] {
                        jobs.push((qi, ri));
                    }
                }
            }
            jobs
        } else {
            let mut jobs = Vec::with_capacity(q_indexed.len() * ref_indexed.len());
            for qi in 0..q_indexed.len() {
                for ri in 0..ref_indexed.len() {
                    jobs.push((qi, ri));
                }
            }
            jobs
        };

        let rank = self.next_rank;
        let max_sv = if self.max_sv_bp == 0 {
            DEFAULT_MAX_SV_BP
        } else {
            self.max_sv_bp
        };
        let cfg = AlignConfig {
            skip_softmasked: skip_soft,
            ..AlignConfig::default()
        };
        let aligned_pairs: Vec<(usize, usize, Vec<AlignedRun>)> = pair_jobs
            .par_iter()
            .filter_map(|&(qi, ri)| {
                let r_idx = &ref_indexed[ri].1;
                let q_idx = &q_indexed[qi].1;
                let aln = align_contig_pair_runs(r_idx, q_idx, &cfg).ok()?;
                let mut runs =
                    net_one_to_one_runs_scored(&aln.aligned_runs, r_idx.seq(), q_idx.seq());
                if runs.len() < 2 {
                    return None;
                }
                runs.sort_by_key(|r| r.x_start);
                Some((qi, ri, runs))
            })
            .collect();

        let existing_other_rank_names: HashSet<String> = self
            .segments
            .iter()
            .filter(|s| s.rank != rank)
            .map(|s| s.origin_name.clone())
            .collect();

        let mut sv_grafted = 0usize;
        for (qi, ri, runs) in aligned_pairs {
            let (raw_query_name, q_idx) = &q_indexed[qi];
            let query_origin = if existing_other_rank_names.contains(raw_query_name) {
                format!("r{rank}#{raw_query_name}")
            } else {
                raw_query_name.clone()
            };
            let (ref_name, r_idx) = &ref_indexed[ri];
            let query_seq = q_idx.seq();
            let ref_seq = r_idx.seq();
            {
                for pair in runs.windows(2) {
                    let prev = pair[0];
                    let cur = pair[1];
                    if prev.is_rc != cur.is_rc {
                        continue;
                    }
                    let prev_ref_end = prev.x_start + prev.len;
                    let mut cur_ref_start = cur.x_start;
                    let mut ref_trim = 0usize;
                    if cur_ref_start < prev_ref_end {
                        let ov = prev_ref_end - cur_ref_start;
                        if ov >= cur.len || ov > 64 {
                            continue;
                        }
                        cur_ref_start = prev_ref_end;
                        ref_trim = ov;
                    }

                    // Extract forward query interval [q_start, q_end) per rGFA spec.
                    let (q_start, q_end, sv_fwd) = if prev.is_rc {
                        // On reverse strand, y decreases as x increases; ref_trim shifts cur_q_hi down.
                        let prev_q_lo = prev.y_start.saturating_add(1).saturating_sub(prev.len);
                        let mut cur_q_hi = cur.y_start.saturating_add(1).saturating_sub(ref_trim);
                        if prev_q_lo < cur_q_hi {
                            let ov = cur_q_hi - prev_q_lo;
                            if ov >= cur.len.saturating_sub(ref_trim) || ov > 64 {
                                continue;
                            }
                            cur_q_hi = prev_q_lo;
                            cur_ref_start = cur_ref_start.saturating_add(ov).min(ref_seq.len());
                        }
                        if prev_q_lo > query_seq.len() {
                            continue;
                        }
                        (cur_q_hi, prev_q_lo, false)
                    } else {
                        let prev_q_end = prev.y_start + prev.len;
                        let mut cur_q_start = cur.y_start.saturating_add(ref_trim);
                        if cur_q_start < prev_q_end {
                            let ov = prev_q_end - cur_q_start;
                            if ov >= cur.len.saturating_sub(ref_trim) || ov > 64 {
                                continue;
                            }
                            cur_q_start = prev_q_end;
                            cur_ref_start = cur_ref_start.saturating_add(ov).min(ref_seq.len());
                        }
                        if cur_q_start > query_seq.len() {
                            continue;
                        }
                        (prev_q_end, cur_q_start, true)
                    };

                    let ref_gap = cur_ref_start.saturating_sub(prev_ref_end);
                    let q_gap = q_end.saturating_sub(q_start);
                    if (ref_gap < self.min_sv_bp && q_gap < self.min_sv_bp)
                        || ref_gap > max_sv
                        || q_gap > max_sv
                    {
                        continue;
                    }

                    // Split rank-0 reference backbone at both SV breakpoints.
                    self.split_rank0_at(ref_name, prev_ref_end);
                    self.split_rank0_at(ref_name, cur_ref_start);

                    let (Some(left_id), Some(right_id)) = (
                        self.rank0_ending_at(ref_name, prev_ref_end),
                        self.rank0_starting_at(ref_name, cur_ref_start),
                    ) else {
                        continue;
                    };

                    if q_gap >= self.min_sv_bp {
                        let sv_id = self.alloc_seg_id();
                        // Always store forward query_seq[q_start..q_end] to match SN/SO per rGFA spec.
                        self.segments.push(RgfaSegment {
                            id: sv_id.clone(),
                            seq: query_seq[q_start..q_end].to_vec(),
                            origin_name: query_origin.clone(),
                            origin_offset: q_start,
                            rank,
                        });
                        self.add_link_dedup(&left_id, true, &sv_id, sv_fwd);
                        self.add_link_dedup(&sv_id, sv_fwd, &right_id, true);
                    } else {
                        // Pure deletion >= min_sv_bp: direct bypass edge across deleted ref interval.
                        self.add_link_dedup(&left_id, true, &right_id, true);
                    }
                    sv_grafted += 1;
                }
                sv_grafted += self.graft_inversion_spans(ref_name, &runs, max_sv);
            }
        }

        self.next_rank = self.next_rank.saturating_add(1);
        sv_grafted
    }

    /// Detect balanced strand-flip inversion spans bracketed by opposite-strand flanking runs and
    /// graft canonical rGFA orientation-flip links (`left+ -> inv_last-` and `inv_first- -> right+`).
    fn graft_inversion_spans(
        &mut self,
        ref_name: &str,
        runs: &[AlignedRun],
        max_sv: usize,
    ) -> usize {
        if runs.len() < 3 {
            return 0;
        }
        // Group consecutive same-strand runs (sorted by `x_start`) into maximal strand blocks.
        let mut blocks: Vec<(usize, usize, bool)> = Vec::new();
        for r in runs {
            let r_end = r.x_start.saturating_add(r.len);
            if let Some(last) = blocks.last_mut()
                && last.2 == r.is_rc
                && r.x_start <= last.1.saturating_add(max_sv)
            {
                last.1 = last.1.max(r_end);
            } else {
                blocks.push((r.x_start, r_end, r.is_rc));
            }
        }

        // Identify the contig's majority strand so intervening forward spacers between two
        // consecutive inversions are never treated as inverted blocks.
        let (fwd_bp, rc_bp) = runs.iter().fold((0usize, 0usize), |(f, r), run| {
            if run.is_rc {
                (f, r.saturating_add(run.len))
            } else {
                (f.saturating_add(run.len), r)
            }
        });
        let majority_rc = rc_bp > fwd_bp;

        let mut grafted = 0usize;
        for win in blocks.windows(3) {
            let (_, left_end, left_rc) = win[0];
            let (raw_inv_start, raw_inv_end, inv_rc) = win[1];
            let (right_start, _, right_rc) = win[2];
            if left_rc != majority_rc || right_rc != majority_rc || inv_rc == majority_rc {
                continue;
            }
            let inv_start = raw_inv_start.max(left_end);
            let inv_end = raw_inv_end.min(right_start);
            if inv_end <= inv_start {
                continue;
            }
            let inv_len = inv_end - inv_start;
            let left_gap = inv_start.saturating_sub(left_end);
            let right_gap = right_start.saturating_sub(inv_end);
            if inv_len < self.min_sv_bp
                || inv_len > max_sv
                || left_gap > max_sv
                || right_gap > max_sv
            {
                continue;
            }
            self.split_rank0_at(ref_name, left_end);
            self.split_rank0_at(ref_name, inv_start);
            self.split_rank0_at(ref_name, inv_end);
            self.split_rank0_at(ref_name, right_start);

            let (Some(left_id), Some(inv_first_id), Some(inv_last_id), Some(right_id)) = (
                self.rank0_ending_at(ref_name, left_end),
                self.rank0_starting_at(ref_name, inv_start),
                self.rank0_ending_at(ref_name, inv_end),
                self.rank0_starting_at(ref_name, right_start),
            ) else {
                continue;
            };
            // Canonical bidirected rGFA inversion links: entering inverted block on `-` enters its
            // highest-offset rank-0 piece (`inv_last_id-`) and exits its lowest-offset piece (`inv_first_id-`).
            self.add_link_dedup(&left_id, true, &inv_last_id, false);
            self.add_link_dedup(&inv_first_id, false, &right_id, true);
            grafted += 1;
        }
        grafted
    }

    /// Incrementally graft a single query contig (`rank = self.next_rank`) onto the rGFA graph.
    pub fn graft_assembly(&mut self, query_name: &str, query_seq: &[u8]) -> usize {
        if query_seq.is_empty() {
            return 0;
        }
        self.graft_assembly_contigs(&[(query_name, query_seq)])
    }

    /// Write standard rGFA (`H`, `S` with `SN:Z:`, `SO:i:`, `SR:i:`, and `L` lines).
    ///
    /// # Errors
    /// Propagates any write error on `out`.
    pub fn write_rgfa<W: Write>(&self, out: &mut W) -> Result<()> {
        writeln!(out, "H\tVN:Z:1.0")?;
        for seg in &self.segments {
            writeln!(
                out,
                "S\t{}\t{}\tSN:Z:{}\tSO:i:{}\tSR:i:{}",
                seg.id,
                String::from_utf8_lossy(&seg.seq),
                seg.origin_name,
                seg.origin_offset,
                seg.rank
            )?;
        }
        for link in &self.links {
            writeln!(
                out,
                "L\t{}\t{}\t{}\t{}\t0M",
                link.from,
                if link.from_fwd { '+' } else { '-' },
                link.to,
                if link.to_fwd { '+' } else { '-' }
            )?;
        }
        Ok(())
    }

    /// Parse an rGFA stream into an [`RgfaGraph`].
    ///
    /// # Errors
    /// Returns an error on malformed rGFA lines, missing `SN`/`SO`/`SR` tags, coordinate overflow,
    /// invalid link orientations, or lines exceeding the 256 MiB safety cap.
    pub fn read_rgfa<R: BufRead>(mut reader: R) -> Result<Self> {
        const MAX_RGFA_LINE_BYTES: usize = 256 * 1024 * 1024;
        let mut segments = Vec::new();
        let mut links = Vec::new();
        let mut link_set = HashSet::new();
        let mut max_rank = 0u32;
        let mut max_seg_id = 0usize;
        let mut raw_line = String::new();
        let mut line_no = 0usize;

        loop {
            let n = crate::io_util::read_bounded_line(
                &mut reader,
                &mut raw_line,
                MAX_RGFA_LINE_BYTES,
                "rGFA",
            )?;
            if n == 0 {
                break;
            }
            line_no += 1;
            let line = raw_line.trim_end_matches(['\r', '\n']);
            if line.is_empty() || line == "H" || line.starts_with("H\t") || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            match fields[0] {
                "S" => {
                    anyhow::ensure!(fields.len() >= 3, "malformed rGFA S line at line {line_no}");
                    let mut origin_name = None;
                    let mut origin_offset = None;
                    let mut rank = None;
                    for &tag in &fields[3..] {
                        if let Some(v) = tag.strip_prefix("SN:Z:") {
                            origin_name = Some(v.to_owned());
                        } else if let Some(v) = tag.strip_prefix("SO:i:") {
                            origin_offset =
                                Some(v.parse::<usize>().context("invalid SO:i: tag in rGFA")?);
                        } else if let Some(v) = tag.strip_prefix("SR:i:") {
                            rank = Some(v.parse::<u32>().context("invalid SR:i: tag in rGFA")?);
                        }
                    }
                    let origin_name =
                        origin_name.with_context(|| format!("missing SN:Z: on line {line_no}"))?;
                    let origin_offset = origin_offset
                        .with_context(|| format!("missing SO:i: on line {line_no}"))?;
                    let rank = rank.with_context(|| format!("missing SR:i: on line {line_no}"))?;
                    let seq_bytes = fields[2].as_bytes().to_vec();
                    anyhow::ensure!(
                        origin_offset.checked_add(seq_bytes.len()).is_some(),
                        "rGFA SO:i: coordinate overflow on line {line_no}"
                    );
                    if let Some(num) = fields[1]
                        .strip_prefix('s')
                        .and_then(|s| s.parse::<usize>().ok())
                    {
                        max_seg_id = max_seg_id.max(num);
                    }
                    max_rank = max_rank.max(rank);
                    segments.push(RgfaSegment {
                        id: fields[1].to_owned(),
                        seq: seq_bytes,
                        origin_name,
                        origin_offset,
                        rank,
                    });
                }
                "L" => {
                    anyhow::ensure!(fields.len() >= 5, "malformed rGFA L line at line {line_no}");
                    let parse_orient = |s: &str| match s {
                        "+" => Ok(true),
                        "-" => Ok(false),
                        other => bail!("invalid rGFA link orientation `{other}` on line {line_no}"),
                    };
                    let edge = RgfaLink {
                        from: fields[1].to_owned(),
                        from_fwd: parse_orient(fields[2])?,
                        to: fields[3].to_owned(),
                        to_fwd: parse_orient(fields[4])?,
                    };
                    if link_set.insert(edge.clone()) {
                        links.push(edge);
                    }
                }
                other => bail!("unsupported rGFA record `{other}` on line {line_no}"),
            }
        }
        let next_seg_num = max_seg_id.max(segments.len()).saturating_add(1);
        Ok(Self {
            segments,
            links,
            min_sv_bp: 50,
            max_sv_bp: DEFAULT_MAX_SV_BP,
            skip_softmasked: false,
            next_rank: max_rank.saturating_add(1),
            next_seg_num,
            link_set,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeding::revcomp;
    use crate::splitmix64;

    fn synth_dna(seed: u64, len: usize) -> Vec<u8> {
        let bases = *b"ACGT";
        let mut st = seed;
        (0..len)
            .map(|_| {
                st = splitmix64(st);
                bases[(st & 3) as usize]
            })
            .collect()
    }

    #[test]
    fn rgfa_grafts_insertion_and_deletion_svs_with_exact_sn_so_sr_tags_and_roundtrips() {
        let ref_seq = synth_dna(0xABCD_1234, 2_000);
        // Query 1: 120-bp insertion at ref pos 600 (rank 1) + tiny 3-bp micro-indel at pos 1400 (ignored when min_sv_bp=50!)
        let ins_payload = synth_dna(0x9999_0001, 120);
        let mut q1 = Vec::new();
        q1.extend_from_slice(&ref_seq[..600]);
        q1.extend_from_slice(&ins_payload);
        q1.extend_from_slice(&ref_seq[600..1400]);
        q1.extend_from_slice(b"ACG");
        q1.extend_from_slice(&ref_seq[1400..]);

        // Query 2: 150-bp deletion of ref [1100..1250) (rank 2)
        let mut q2 = Vec::new();
        q2.extend_from_slice(&ref_seq[..1100]);
        q2.extend_from_slice(&ref_seq[1250..]);

        let mut g = RgfaGraph::from_reference(&[("chr1", &ref_seq)], 50);
        let sv1 = g.graft_assembly("HG002#1#chr1", &q1);
        let sv2 = g.graft_assembly("HG003#1#chr1", &q2);

        assert!(
            sv1 >= 1,
            "expected >= 1 grafted SV for 120bp insertion, got {sv1}"
        );
        assert!(
            sv2 >= 1,
            "expected >= 1 grafted SV for 150bp deletion, got {sv2}"
        );

        let mut rank0: Vec<&RgfaSegment> = g.segments.iter().filter(|s| s.rank == 0).collect();
        rank0.sort_by_key(|s| s.origin_offset);
        let ref_recon: Vec<u8> = rank0.iter().flat_map(|s| s.seq.iter().copied()).collect();
        assert_eq!(
            ref_recon, ref_seq,
            "rank-0 segments must reconstruct reference exactly"
        );

        let rank1: Vec<&RgfaSegment> = g.segments.iter().filter(|s| s.rank == 1).collect();
        assert_eq!(rank1.len(), 1, "expected 1 rank-1 SV segment (120bp ins)");
        assert_eq!(rank1[0].origin_name, "HG002#1#chr1");
        assert!(
            (110..=130).contains(&rank1[0].seq.len()),
            "unexpected SV segment length {}",
            rank1[0].seq.len()
        );

        let mut buf = Vec::new();
        g.write_rgfa(&mut buf).expect("write rgfa");
        let parsed = RgfaGraph::read_rgfa(buf.as_slice()).expect("read rgfa");
        assert_eq!(parsed.segments, g.segments);
        assert_eq!(parsed.links, g.links);
    }

    #[test]
    fn rgfa_multi_contig_assembly_uses_single_rank_and_revcomp_sv_preserves_forward_sn_so() {
        let ref_seq = synth_dna(0xFEED_FACE, 2_000);
        let ins = synth_dna(0xBEEF_1111, 100);
        let mut fwd_q = Vec::new();
        fwd_q.extend_from_slice(&ref_seq[..800]);
        fwd_q.extend_from_slice(&ins);
        fwd_q.extend_from_slice(&ref_seq[800..]);
        let rc_q = revcomp(&fwd_q);

        let mut g = RgfaGraph::from_reference(&[("chr1", &ref_seq)], 50);
        let svs = g.graft_assembly_contigs(&[("ctg_rc", &rc_q), ("ctg_empty", b"")]);
        assert!(svs >= 1, "expected reverse-strand insertion SV");
        assert_eq!(
            g.next_rank, 2,
            "multi-contig assembly must increment rank exactly once"
        );
        let sv_seg = g
            .segments
            .iter()
            .find(|s| s.rank == 1)
            .expect("rank-1 SV segment");
        assert_eq!(
            sv_seg.seq,
            rc_q[sv_seg.origin_offset..sv_seg.origin_offset + sv_seg.seq.len()],
            "rGFA SV segment sequence must match forward slice [SO, SO+len) of SN"
        );
        assert!(
            g.links
                .iter()
                .any(|l| l.to == sv_seg.id && !l.to_fwd && l.from_fwd),
            "reverse-strand SV segment must be linked with '-' orientation"
        );
    }

    #[test]
    fn read_rgfa_rejects_overflow_bad_orient_and_avoids_seg_id_collision() {
        let bad_so = format!("S\ts1\tACGT\tSN:Z:c1\tSO:i:{}\tSR:i:0\n", usize::MAX);
        assert!(RgfaGraph::read_rgfa(bad_so.as_bytes()).is_err());

        let bad_or = "S\ts1\tACGT\tSN:Z:c1\tSO:i:0\tSR:i:0\nL\ts1\t?\ts1\t+\t0M\n";
        assert!(RgfaGraph::read_rgfa(bad_or.as_bytes()).is_err());

        let gapped_ids =
            "S\ts1\tACGT\tSN:Z:c1\tSO:i:0\tSR:i:0\nS\ts5\tTGCA\tSN:Z:c1\tSO:i:4\tSR:i:0\n";
        let mut g = RgfaGraph::read_rgfa(gapped_ids.as_bytes()).expect("valid rgfa");
        assert_eq!(g.alloc_seg_id(), "s6");
    }

    #[test]
    fn rgfa_grafts_balanced_inversion_strand_flip_edges() {
        let ref_seq = synth_dna(0xCAFE_BEEF, 3_000);
        let mut q_inv = Vec::with_capacity(ref_seq.len());
        q_inv.extend_from_slice(&ref_seq[..1_200]);
        q_inv.extend_from_slice(&revcomp(&ref_seq[1_200..1_800]));
        q_inv.extend_from_slice(&ref_seq[1_800..]);

        let mut g = RgfaGraph::from_reference(&[("chr1", &ref_seq)], 50);
        let svs = g.graft_assembly("inv_hap", &q_inv);
        assert!(
            svs >= 1,
            "expected >= 1 grafted SV for 600 bp balanced inversion, got {svs}"
        );

        let mut rank0: Vec<&RgfaSegment> = g.segments.iter().filter(|s| s.rank == 0).collect();
        rank0.sort_by_key(|s| s.origin_offset);
        let ref_recon: Vec<u8> = rank0.iter().flat_map(|s| s.seq.iter().copied()).collect();
        assert_eq!(ref_recon, ref_seq);

        assert!(
            g.links.iter().any(|l| l.from_fwd && !l.to_fwd),
            "expected forward-to-reverse inversion entry link (+ -> -), links: {:?}",
            g.links
        );
        assert!(
            g.links.iter().any(|l| !l.from_fwd && l.to_fwd),
            "expected reverse-to-forward inversion exit link (- -> +), links: {:?}",
            g.links
        );
    }

    #[test]
    fn rgfa_many_left_to_right_splits_do_not_leak_quadratic_vec_capacity() {
        let ref_seq = synth_dna(0x1234_5678, 100_000);
        let mut g = RgfaGraph::from_reference(&[("chr1", &ref_seq)], 50);
        for pos in (200..99_800).step_by(200) {
            g.split_rank0_at("chr1", pos);
        }
        let total_len: usize = g.segments.iter().map(|s| s.seq.len()).sum();
        let total_cap: usize = g.segments.iter().map(|s| s.seq.capacity()).sum();
        assert_eq!(total_len, ref_seq.len());
        assert!(
            total_cap <= total_len.saturating_mul(2).saturating_add(4096),
            "quadratic Vec::split_off capacity leak: total_cap={total_cap} vs total_len={total_len}"
        );
    }

    #[test]
    fn rgfa_two_consecutive_inversions_do_not_flip_intervening_spacer() {
        let ref_seq = synth_dna(0xDEAF_BEEF, 6_000);
        let mut q = Vec::with_capacity(ref_seq.len());
        q.extend_from_slice(&ref_seq[..1_000]);
        q.extend_from_slice(&revcomp(&ref_seq[1_000..1_600]));
        q.extend_from_slice(&ref_seq[1_600..3_600]);
        q.extend_from_slice(&revcomp(&ref_seq[3_600..4_200]));
        q.extend_from_slice(&ref_seq[4_200..]);

        let mut g = RgfaGraph::from_reference(&[("chr1", &ref_seq)], 50);
        let svs = g.graft_assembly("inv2_hap", &q);
        assert!(svs >= 2, "expected >= 2 inversions grafted, got {svs}");

        // Exactly two (+ -> -) entry links and two (- -> +) exit links; the 2,000 bp forward
        // spacer [1600..3600) must NEVER be linked with `-` orientation!
        let seg_by_id: std::collections::HashMap<&str, &RgfaSegment> =
            g.segments.iter().map(|s| (s.id.as_str(), s)).collect();
        for l in &g.links {
            if !l.to_fwd {
                let target = seg_by_id[l.to.as_str()];
                assert!(
                    target.seq.len() <= 800,
                    "forward spacer falsely linked on '-' strand: {:?} (len={}, so={})",
                    l,
                    target.seq.len(),
                    target.origin_offset
                );
            }
            if !l.from_fwd {
                let src = seg_by_id[l.from.as_str()];
                assert!(
                    src.seq.len() <= 800,
                    "forward spacer falsely exited on '-' strand: {:?} (len={}, so={})",
                    l,
                    src.seq.len(),
                    src.origin_offset
                );
            }
        }
    }
}
