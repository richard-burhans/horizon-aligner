//! Multi-level streaming HAL Newick parsing, chunked BED liftover, and MAF / block projection
//! for [`StreamingHalFile`].

use super::stream_hal::{HalLiftoverHit, HalMafStats, LiftedFrag, StreamingHalFile, TreeStep};
use std::collections::HashSet;
use std::io::{Read, Seek, Write};

pub(super) struct TreeNode {
    pub(super) name: String,
    pub(super) parent: Option<usize>,
    pub(super) children: Vec<usize>,
}

const MAX_NEWICK_DEPTH: usize = 4_096;
const MAX_NEWICK_NODES: usize = 65_536;

fn skip_ws(b: &[u8], p: &mut usize) {
    while *p < b.len() && b[*p].is_ascii_whitespace() {
        *p += 1;
    }
}

fn parse_label(b: &[u8], p: &mut usize) -> String {
    skip_ws(b, p);
    let start = *p;
    while *p < b.len() && !matches!(b[*p], b'(' | b')' | b',' | b':' | b';') {
        *p += 1;
    }
    let name = String::from_utf8_lossy(&b[start..*p]).trim().to_owned();
    skip_ws(b, p);
    if *p < b.len() && b[*p] == b':' {
        *p += 1;
        while *p < b.len() && !matches!(b[*p], b'(' | b')' | b',' | b';') {
            *p += 1;
        }
    }
    name
}

/// Parse a HAL Newick tree string (with optional branch lengths `:len` and internal node names)
/// into a flat node list preserving exact Newick child order, enforcing depth/node safety caps
/// iteratively on the heap so deep inputs never overflow the call stack.
pub(super) fn parse_hal_newick(newick: &str) -> std::io::Result<Vec<TreeNode>> {
    enum ParseState {
        ExpectNode,
        AfterNode,
    }

    let bytes = newick.trim().as_bytes();
    let mut nodes: Vec<TreeNode> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    let mut pos = 0usize;
    if bytes.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "empty HAL Phylogeny Newick string",
        ));
    }

    let mut state = ParseState::ExpectNode;
    loop {
        match state {
            ParseState::ExpectNode => {
                if stack.len() > MAX_NEWICK_DEPTH {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "HAL Newick tree exceeds max nesting depth",
                    ));
                }
                if nodes.len() >= MAX_NEWICK_NODES {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "HAL Newick tree exceeds max node count",
                    ));
                }
                skip_ws(bytes, &mut pos);
                let idx = nodes.len();
                let parent = stack.last().copied();
                nodes.push(TreeNode {
                    name: String::new(),
                    parent,
                    children: Vec::new(),
                });
                if let Some(p_idx) = parent {
                    nodes[p_idx].children.push(idx);
                }
                if pos < bytes.len() && bytes[pos] == b'(' {
                    pos += 1;
                    stack.push(idx);
                } else {
                    let name = parse_label(bytes, &mut pos);
                    nodes[idx].name = if name.is_empty() {
                        format!("Anc_{idx}")
                    } else {
                        name
                    };
                    state = ParseState::AfterNode;
                }
            }
            ParseState::AfterNode => {
                skip_ws(bytes, &mut pos);
                let Some(&parent_idx) = stack.last() else {
                    break;
                };
                if pos < bytes.len() && bytes[pos] == b',' {
                    pos += 1;
                    state = ParseState::ExpectNode;
                } else if pos < bytes.len() && bytes[pos] == b')' {
                    pos += 1;
                    stack.pop();
                    let name = parse_label(bytes, &mut pos);
                    nodes[parent_idx].name = if name.is_empty() {
                        format!("Anc_{parent_idx}")
                    } else {
                        name
                    };
                } else {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "malformed HAL Newick tree: missing ')'",
                    ));
                }
            }
        }
    }

    skip_ws(bytes, &mut pos);
    if pos < bytes.len() && bytes[pos] == b';' {
        pos += 1;
    }
    skip_ws(bytes, &mut pos);
    if pos < bytes.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "trailing garbage in HAL Phylogeny Newick string",
        ));
    }
    Ok(nodes)
}

/// One gap-free projected row emitted by [`StreamingHalFile::for_each_projected_block`].
#[derive(Debug, Clone)]
pub struct HalProjectedRow {
    /// Genome / species name.
    pub species: String,
    /// Sequence / contig name within `species`.
    pub contig: String,
    /// 0-based UCSC MAF start coordinate (relative to `strand`).
    pub start: u64,
    /// Aligned ungapped length in base pairs.
    pub size: u64,
    /// Strand (`'+'` or `'-'`).
    pub strand: char,
    /// Total length of `contig`.
    pub src_size: u64,
    /// Oriented ASCII DNA sequence of length `size`.
    pub seq: Vec<u8>,
}

/// One projected elementary gap-free alignment block (reference row + `>= 1` target rows).
#[derive(Debug, Clone)]
pub struct HalProjectedBlock {
    /// Reference genome row (always on `'+'` strand).
    pub ref_row: HalProjectedRow,
    /// Projected target genome rows spanning the exact same reference sub-interval.
    pub target_rows: Vec<HalProjectedRow>,
}

impl<R: Read + Seek> StreamingHalFile<R> {
    /// Liftover a 0-based half-open interval `[start, end)` on `(src_genome, src_seq)` to
    /// `dst_genome` across arbitrary multi-level tree paths (chunked in 64 kb windows so large
    /// BED intervals never silently truncate).
    ///
    /// # Errors
    /// Returns an `std::io::Error` if either genome/sequence is missing or HDF5 reading fails.
    pub fn liftover_interval(
        &mut self,
        src_genome: &str,
        src_seq: &str,
        start: u64,
        end: u64,
        dst_genome: &str,
    ) -> std::io::Result<Vec<HalLiftoverHit>> {
        const LIFTOVER_CHUNK_BP: u64 = 65_536;
        let src_idx = *self.name_to_idx.get(src_genome).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("source genome `{src_genome}` not found in HAL"),
            )
        })?;
        let dst_idx = *self.name_to_idx.get(dst_genome).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("destination genome `{dst_genome}` not found in HAL"),
            )
        })?;
        self.ensure_genome_sequences(src_idx)?;
        self.ensure_genome_sequences(dst_idx)?;
        let seq_idx = *self.genomes[src_idx]
            .seq_name_to_idx
            .get(src_seq)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    format!("sequence `{src_seq}` not found in genome `{src_genome}`"),
                )
            })?;
        let seq_meta = &self.genomes[src_idx].sequences[seq_idx];
        let clamped_s = start.min(seq_meta.length);
        let clamped_e = end.min(seq_meta.length);
        if clamped_e <= clamped_s {
            return Ok(Vec::new());
        }
        let g_s = seq_meta.offset.saturating_add(clamped_s);
        let g_e = seq_meta.offset.saturating_add(clamped_e);
        let path = self.find_tree_path(src_idx, dst_idx);
        let mut all_lifted: Vec<LiftedFrag> = Vec::new();
        let mut cur_s = g_s;
        while cur_s < g_e {
            let cur_e = cur_s.saturating_add(LIFTOVER_CHUNK_BP).min(g_e);
            all_lifted.extend(self.liftover_global(&path, cur_s, cur_e)?);
            cur_s = cur_e;
        }
        let dst_g = &self.genomes[dst_idx];
        let mut hits: Vec<HalLiftoverHit> = Vec::with_capacity(all_lifted.len());
        for f in all_lifted {
            let mut ds = f.dst_s;
            while ds < f.dst_e {
                let Some((dseq, ls, le)) = dst_g.map_global_to_seq(ds, f.dst_e) else {
                    break;
                };
                let strand = if f.rev { '-' } else { '+' };
                if let Some(last) = hits.last_mut()
                    && last.seq_name == dseq.name
                    && last.strand == strand
                    && ((!f.rev && last.end == ls) || (f.rev && le == last.start))
                {
                    if f.rev {
                        last.start = ls;
                    } else {
                        last.end = le;
                    }
                } else {
                    hits.push(HalLiftoverHit {
                        seq_name: dseq.name.clone(),
                        start: ls,
                        end: le,
                        seq_len: dseq.length,
                        strand,
                    });
                }
                ds = ds.saturating_add(le - ls);
            }
        }
        Ok(hits)
    }

    /// Stream elementary multi-level projected alignment blocks onto `ref_genome`, partitioning
    /// each reference segment at all lifted target `[src_s, src_e)` breakpoints so multi-level
    /// indels/rearrangements never silently drop target alignments.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if any genome/sequence is missing, HDF5 reading fails, or
    /// `on_block` returns an error.
    #[allow(clippy::too_many_lines)]
    pub fn for_each_projected_block<F>(
        &mut self,
        ref_genome: &str,
        seq_slice: Option<(&str, u64, Option<u64>)>,
        targets: Option<&[String]>,
        mut on_block: F,
    ) -> std::io::Result<HalMafStats>
    where
        F: FnMut(HalProjectedBlock) -> std::io::Result<()>,
    {
        let ref_idx = *self.name_to_idx.get(ref_genome).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("reference genome `{ref_genome}` not found in HAL"),
            )
        })?;
        let target_indices: Vec<usize> = if let Some(req) = targets {
            let mut v = Vec::with_capacity(req.len());
            let mut seen = HashSet::with_capacity(req.len());
            for name in req {
                let &idx = self.name_to_idx.get(name).ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("target genome `{name}` not found in HAL"),
                    )
                })?;
                if idx != ref_idx && seen.insert(idx) {
                    v.push(idx);
                }
            }
            v
        } else {
            let mut leaves: Vec<usize> = self
                .genomes
                .iter()
                .enumerate()
                .filter_map(|(i, g)| {
                    (i != ref_idx && g.children.is_empty() && g.group_addr != u64::MAX).then_some(i)
                })
                .collect();
            if leaves.is_empty() {
                leaves = (0..self.genomes.len())
                    .filter(|&i| i != ref_idx && self.genomes[i].group_addr != u64::MAX)
                    .collect();
            }
            leaves
        };

        self.ensure_genome_sequences(ref_idx)?;
        for &ti in &target_indices {
            self.ensure_genome_sequences(ti)?;
        }

        let paths: Vec<(usize, Vec<TreeStep>)> = target_indices
            .iter()
            .map(|&ti| (ti, self.find_tree_path(ref_idx, ti)))
            .collect();

        let mut ref_windows: Vec<(u64, u64)> = Vec::new();
        if let Some((seq_name, start, len_opt)) = seq_slice {
            let seq_idx = *self.genomes[ref_idx]
                .seq_name_to_idx
                .get(seq_name)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("sequence `{seq_name}` not found in `{ref_genome}`"),
                    )
                })?;
            let sm = &self.genomes[ref_idx].sequences[seq_idx];
            let s = start.min(sm.length);
            let e = len_opt.map_or(sm.length, |l| s.saturating_add(l).min(sm.length));
            if e > s {
                ref_windows.push((sm.offset + s, sm.offset + e));
            }
        } else {
            for sm in &self.genomes[ref_idx].sequences {
                if sm.length > 0 {
                    ref_windows.push((sm.offset, sm.offset + sm.length));
                }
            }
        }

        let mut stats = HalMafStats::default();
        self.ensure_genome_header(ref_idx)?;
        // Project in 16 kb windows directly via `liftover_global` instead of calling
        // `liftover_global` once per tiny reference segment; `liftover_global` already splits
        // across every segment & target breakpoint along the tree path.
        let max_project_window_bp = 16_384u64;

        let mut per_target_frags: Vec<(usize, Vec<LiftedFrag>)> = Vec::with_capacity(paths.len());
        let mut breaks: Vec<u64> = Vec::new();

        for (win_s, win_e) in ref_windows {
            let mut cur_s = win_s;
            while cur_s < win_e {
                let cur_e = cur_s.saturating_add(max_project_window_bp).min(win_e);
                if cur_e <= cur_s {
                    break;
                }
                per_target_frags.clear();
                breaks.clear();
                breaks.push(cur_s);
                breaks.push(cur_e);

                for &(t_idx, ref path) in &paths {
                    let lifted = self.liftover_global(path, cur_s, cur_e)?;
                    if !lifted.is_empty() {
                        for f in &lifted {
                            if f.src_s > cur_s && f.src_s < cur_e {
                                breaks.push(f.src_s);
                            }
                            if f.src_e > cur_s && f.src_e < cur_e {
                                breaks.push(f.src_e);
                            }
                        }
                        per_target_frags.push((t_idx, lifted));
                    }
                }

                if per_target_frags.is_empty() {
                    cur_s = cur_e;
                    continue;
                }
                breaks.sort_unstable();
                breaks.dedup();

                for win in breaks.windows(2) {
                    let sub_s = win[0];
                    let sub_e = win[1];
                    let sub_len = usize::try_from(sub_e - sub_s).unwrap_or(0);
                    if sub_len == 0 {
                        continue;
                    }
                    let Some((r_species, r_contig, r_local_s, r_seq_len)) = self.genomes[ref_idx]
                        .map_global_to_seq(sub_s, sub_e)
                        .and_then(|(s, ls, le)| {
                            (le - ls == sub_len as u64).then(|| {
                                (
                                    self.genomes[ref_idx].name.clone(),
                                    s.name.clone(),
                                    ls,
                                    s.length,
                                )
                            })
                        })
                    else {
                        continue;
                    };

                    let mut target_rows: Vec<HalProjectedRow> = Vec::new();
                    for (t_idx, frags) in &per_target_frags {
                        for f in frags {
                            if f.src_s > sub_s || f.src_e < sub_e {
                                continue;
                            }
                            let off_lo = sub_s - f.src_s;
                            let off_hi = sub_e - f.src_s;
                            let (ts, te) = if f.rev {
                                (f.dst_e - off_hi, f.dst_e - off_lo)
                            } else {
                                (f.dst_s + off_lo, f.dst_s + off_hi)
                            };
                            let Some((t_species, t_contig, t_local_s, t_seq_len)) = self.genomes
                                [*t_idx]
                                .map_global_to_seq(ts, te)
                                .and_then(|(s, ls, le)| {
                                    (le - ls == sub_len as u64).then(|| {
                                        (
                                            self.genomes[*t_idx].name.clone(),
                                            s.name.clone(),
                                            ls,
                                            s.length,
                                        )
                                    })
                                })
                            else {
                                continue;
                            };
                            let mut dna = self.read_dna_slice(*t_idx, ts, sub_len)?;
                            if dna.len() != sub_len {
                                continue;
                            }
                            let (strand, maf_start) = if f.rev {
                                dna = crate::seeding::revcomp(&dna);
                                ('-', t_seq_len.saturating_sub(t_local_s + sub_len as u64))
                            } else {
                                ('+', t_local_s)
                            };
                            target_rows.push(HalProjectedRow {
                                species: t_species,
                                contig: t_contig,
                                start: maf_start,
                                size: sub_len as u64,
                                strand,
                                src_size: t_seq_len,
                                seq: dna,
                            });
                        }
                    }

                    if !target_rows.is_empty() {
                        let ref_dna = self.read_dna_slice(ref_idx, sub_s, sub_len)?;
                        if ref_dna.len() == sub_len {
                            let n_trows = target_rows.len();
                            on_block(HalProjectedBlock {
                                ref_row: HalProjectedRow {
                                    species: r_species,
                                    contig: r_contig,
                                    start: r_local_s,
                                    size: sub_len as u64,
                                    strand: '+',
                                    src_size: r_seq_len,
                                    seq: ref_dna,
                                },
                                target_rows,
                            })?;
                            stats.blocks_emitted += 1;
                            stats.ref_bp_emitted += sub_len as u64;
                            stats.target_rows_emitted += n_trows;
                        }
                    }
                }
                cur_s = cur_e;
            }
        }
        Ok(stats)
    }

    /// Stream multi-level MAF v1 blocks projected onto `ref_genome` (optionally restricted to a
    /// single sequence slice and target genome subset).
    ///
    /// # Errors
    /// Returns an `std::io::Error` if any genome/sequence is missing or I/O fails.
    pub fn project_to_maf<W: Write>(
        &mut self,
        ref_genome: &str,
        seq_slice: Option<(&str, u64, Option<u64>)>,
        targets: Option<&[String]>,
        out: &mut W,
    ) -> std::io::Result<HalMafStats> {
        writeln!(out, "##maf version=1 scoring=horizon_hal")?;
        self.for_each_projected_block(ref_genome, seq_slice, targets, |blk| {
            writeln!(out, "a score={}.0", blk.target_rows.len() + 1)?;
            writeln!(
                out,
                "s {}.{} {} {} {} {} {}",
                blk.ref_row.species,
                blk.ref_row.contig,
                blk.ref_row.start,
                blk.ref_row.size,
                blk.ref_row.strand,
                blk.ref_row.src_size,
                String::from_utf8_lossy(&blk.ref_row.seq),
            )?;
            for t in &blk.target_rows {
                writeln!(
                    out,
                    "s {}.{} {} {} {} {} {}",
                    t.species,
                    t.contig,
                    t.start,
                    t.size,
                    t.strand,
                    t.src_size,
                    String::from_utf8_lossy(&t.seq),
                )?;
            }
            writeln!(out)
        })
    }
}
