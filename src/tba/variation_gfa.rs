//! Native base-level variation GFA v1.1 emission from multi-species TBA/MAF blocks.
//!
//! Decomposes each alignment block into maximal column-equivalence runs: columns sharing the same
//! presence mask and nucleotide equivalence partition compact into shared sequence-bearing `S`
//! segments (with real DNA bases, `SN:Z:`, `SO:i:`, `SR:i:`), while substitutions and indels
//! branch into clean Tarjan snarl bubbles. Emits deduplicated `L` links (`0M`) and exact-span
//! GFA v1.1 `W` walks (`sum(LN) == W_end - W_start`) directly consumable by `horizon_clean`,
//! `horizon_resolve --vcf`, `horizon_lift`, `horizon_wga --graph`, `odgi`, and `vg`.

use super::model::{TbaBlock, for_each_multi_tba_block_reader};
use super::writer::write_maf_block;
use crate::formats::strip_sample_prefix;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as _;
use std::io::{BufRead, Write};

/// Maximum unaligned intra-contig insertion span (bp) spliced as a sequence-bearing `S` segment
/// when source FASTA sequences are provided to [`VariationGfaBuilder`].
pub const MAX_UNALIGNED_INSERTION_SPAN_BP: usize = 200_000;

/// Summary telemetry for native base-level variation GFA emission.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct VariationGfaStats {
    /// Multi-species alignment blocks decomposed (`rows.len() >= 2`).
    pub blocks_processed: usize,
    /// Total sequence-bearing `S` segments written.
    pub segments_emitted: usize,
    /// Total deduplicated `L` links written.
    pub links_emitted: usize,
    /// Total contiguous GFA v1.1 `W` walks written.
    pub walks_emitted: usize,
    /// Number of polymorphic/gapped column runs that branched into `>= 2` alleles or a deletion.
    pub polymorphic_runs: usize,
    /// Total DNA bases across all emitted `S` segments.
    pub total_segment_bp: usize,
}

#[derive(Debug, Clone, Copy)]
struct VarWalkVisit {
    fwd_start: usize,
    fwd_end: usize,
    /// Packed `(seg_id << 1) | u64::from(is_fwd)` so each walk visit occupies exactly 24 bytes.
    seg_dir: u64,
}

impl VarWalkVisit {
    #[inline]
    fn new(fwd_start: usize, fwd_end: usize, seg_id: u64, is_fwd: bool) -> Self {
        Self {
            fwd_start,
            fwd_end,
            seg_dir: (seg_id << 1) | u64::from(is_fwd),
        }
    }

    #[inline]
    fn seg_id(self) -> u64 {
        self.seg_dir >> 1
    }

    #[inline]
    fn is_fwd(self) -> bool {
        (self.seg_dir & 1) != 0
    }
}

/// Split an optional `PanSN` `sample#hap` species label or leading `<hap>#<seq_id>` contig prefix
/// when the haplotype token is a valid non-empty ASCII integer; otherwise default haplotype to `"0"`.
fn split_pansn_fields<'a>(species: &'a str, contig: &'a str) -> (&'a str, &'a str, &'a str) {
    if let Some((sample, hap)) = species.split_once('#')
        && !sample.is_empty()
        && !hap.is_empty()
        && hap.bytes().all(|b| b.is_ascii_digit())
    {
        return (sample, hap, contig);
    }
    if let Some((hap, rest)) = contig.split_once('#')
        && !hap.is_empty()
        && !rest.is_empty()
        && hap.bytes().all(|b| b.is_ascii_digit())
    {
        return (species, hap, rest);
    }
    (species, "0", contig)
}

/// Streaming builder converting [`TbaBlock`]s into a sequence-bearing GFA v1.1 variation graph.
///
/// Emits `S` segments immediately per block while accumulating compact 24-byte walk visits per
/// `(species, contig)` so multi-GB MAF files convert with low peak memory.
pub struct VariationGfaBuilder<W: Write> {
    out: W,
    next_seg_id: u64,
    stats: VariationGfaStats,
    walks: BTreeMap<(String, String), Vec<VarWalkVisit>>,
    /// Deterministic global rGFA rank (`SR:i:`) per species so every `SN:Z:` has a unique rank.
    species_ranks: HashMap<String, usize>,
    /// Optional source FASTA sequences keyed by `(species, stripped_contig)` used to splice
    /// unaligned/singleton insertion spans (`<= MAX_UNALIGNED_INSERTION_SPAN_BP`) into walks.
    source_seqs: HashMap<(String, String), Vec<u8>>,
    // Reusable per-block scratch buffers (zero per-column heap allocation).
    row_offsets: Vec<usize>,
    prev_gids: Vec<u16>,
    cur_gids: Vec<u16>,
    base_to_gid: [u16; 256],
    touched_bases: Vec<u8>,
    group_first_row: Vec<usize>,
    group_counts: Vec<usize>,
    group_seg_ids: Vec<u64>,
    seq_buf: Vec<u8>,
}

/// Compute canonical equivalence group IDs for column `col`: `0` for gap `'-'`, else `1..=K`
/// assigned in row-appearance order by uppercased ASCII base. Returns `K` (`0` if all-gap).
fn classify_column(
    base_to_gid: &mut [u16; 256],
    touched_bases: &mut Vec<u8>,
    blk: &TbaBlock,
    col: usize,
    out_gids: &mut Vec<u16>,
) -> u16 {
    out_gids.clear();
    let mut next_gid = 0u16;
    for row in &blk.rows {
        let raw = row.seq.get(col).copied().unwrap_or(b'-');
        if raw == b'-' {
            out_gids.push(0);
        } else {
            let b = raw.to_ascii_uppercase();
            let mut gid = base_to_gid[b as usize];
            if gid == 0 {
                next_gid += 1;
                gid = next_gid;
                base_to_gid[b as usize] = gid;
                touched_bases.push(b);
            }
            out_gids.push(gid);
        }
    }
    for &b in touched_bases.iter() {
        base_to_gid[b as usize] = 0;
    }
    touched_bases.clear();
    next_gid
}

impl<W: Write> VariationGfaBuilder<W> {
    /// Create a new streaming variation-GFA v1.1 builder writing segments directly to `out`.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if writing the GFA header fails.
    pub fn new(mut out: W) -> std::io::Result<Self> {
        out.write_all(b"H\tVN:Z:1.1\n")?;
        Ok(Self {
            out,
            next_seg_id: 1,
            stats: VariationGfaStats::default(),
            walks: BTreeMap::new(),
            species_ranks: HashMap::new(),
            source_seqs: HashMap::new(),
            row_offsets: Vec::new(),
            prev_gids: Vec::new(),
            cur_gids: Vec::new(),
            base_to_gid: [0u16; 256],
            touched_bases: Vec::with_capacity(8),
            group_first_row: Vec::with_capacity(8),
            group_counts: Vec::with_capacity(8),
            group_seg_ids: Vec::with_capacity(8),
            seq_buf: Vec::new(),
        })
    }

    /// Register source contig sequences `(species, stripped_contig) -> dna_bytes` so unaligned
    /// insertion spans between consecutive aligned blocks stay spliced into continuous walks.
    pub fn add_source_contigs(&mut self, species: &str, contigs: &[(String, Vec<u8>)]) {
        self.rank_for_species(species);
        for (ctg, seq) in contigs {
            let stripped = strip_sample_prefix(ctg, species);
            self.source_seqs
                .insert((species.to_owned(), stripped.to_owned()), seq.clone());
        }
    }

    fn rank_for_species(&mut self, species: &str) -> usize {
        if let Some(&r) = self.species_ranks.get(species) {
            return r;
        }
        let next = self.species_ranks.len();
        self.species_ranks.insert(species.to_owned(), next);
        next
    }

    /// Feed one aligned [`TbaBlock`] into the variation graph, immediately streaming `S` lines.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if writing segment lines fails.
    #[allow(clippy::too_many_lines)]
    pub fn push_block(&mut self, blk: &TbaBlock) -> std::io::Result<()> {
        let n_rows = blk.rows.len();
        let n_cols = blk.n_cols();
        if n_rows < 2 || n_cols == 0 || blk.rows.iter().any(|r| r.seq.len() != n_cols) {
            return Ok(());
        }
        self.stats.blocks_processed += 1;
        self.row_offsets.clear();
        self.row_offsets.resize(n_rows, 0);

        let mut col = 0usize;
        while col < n_cols {
            let num_groups = usize::from(classify_column(
                &mut self.base_to_gid,
                &mut self.touched_bases,
                blk,
                col,
                &mut self.prev_gids,
            ));
            let mut end_col = col + 1;
            while end_col < n_cols {
                let _ = classify_column(
                    &mut self.base_to_gid,
                    &mut self.touched_bases,
                    blk,
                    end_col,
                    &mut self.cur_gids,
                );
                if self.cur_gids != self.prev_gids {
                    break;
                }
                end_col += 1;
            }
            let run_len = end_col - col;
            if num_groups == 0 {
                col = end_col;
                continue;
            }

            let has_gaps = self.prev_gids.contains(&0);
            if num_groups > 1 || has_gaps {
                self.stats.polymorphic_runs += 1;
            }

            // Record first representative row and member count for each allele group `1..=num_groups`.
            self.group_first_row.clear();
            self.group_first_row.resize(num_groups, usize::MAX);
            self.group_counts.clear();
            self.group_counts.resize(num_groups, 0);
            for (ri, &gid) in self.prev_gids.iter().enumerate() {
                if gid > 0 {
                    let g = usize::from(gid - 1);
                    if self.group_first_row[g] == usize::MAX {
                        self.group_first_row[g] = ri;
                    }
                    self.group_counts[g] += 1;
                }
            }

            // Emit one sequence-bearing `S` segment per distinct allele group.
            self.group_seg_ids.clear();
            self.group_seg_ids.resize(num_groups, 0);
            for g in 0..num_groups {
                let rep_ri = self.group_first_row[g];
                let rep = &blk.rows[rep_ri];
                let seg_id = self.next_seg_id;
                self.next_seg_id += 1;
                self.group_seg_ids[g] = seg_id;

                let sr_rank = self.rank_for_species(&rep.species);
                self.seq_buf.clear();
                self.seq_buf
                    .extend(rep.seq[col..end_col].iter().map(u8::to_ascii_uppercase));
                let seq_str = std::str::from_utf8(&self.seq_buf).map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("non-UTF-8 DNA segment {seg_id}: {e}"),
                    )
                })?;
                let rep_ctg = strip_sample_prefix(&rep.contig, &rep.species);
                let (a_fwd_s, a_fwd_e) = rep.fwd_interval();
                let a_so = if rep.strand == '-' {
                    a_fwd_e.saturating_sub(self.row_offsets[rep_ri] + run_len)
                } else {
                    a_fwd_s + self.row_offsets[rep_ri]
                };
                writeln!(
                    self.out,
                    "S\t{seg_id}\t{seq_str}\tLN:i:{run_len}\tSP:i:{}\tSN:Z:{}.{rep_ctg}\tSO:i:{a_so}\tSR:i:{sr_rank}",
                    self.group_counts[g], rep.species,
                )?;
                self.stats.segments_emitted += 1;
                self.stats.total_segment_bp += run_len;
            }

            // Record exact forward-coordinate walk visits for every present row.
            for (ri, row) in blk.rows.iter().enumerate() {
                let gid = self.prev_gids[ri];
                if gid == 0 {
                    continue;
                }
                let seg_id = self.group_seg_ids[usize::from(gid - 1)];
                let stripped = strip_sample_prefix(&row.contig, &row.species);
                let (fwd_s, fwd_e) = row.fwd_interval();
                let is_fwd = row.strand != '-';
                let (sub_s, sub_e) = if is_fwd {
                    let s = fwd_s + self.row_offsets[ri];
                    (s, s + run_len)
                } else {
                    let e = fwd_e.saturating_sub(self.row_offsets[ri]);
                    let s = e.saturating_sub(run_len);
                    (s, e)
                };
                self.row_offsets[ri] += run_len;
                self.walks
                    .entry((row.species.clone(), stripped.to_owned()))
                    .or_default()
                    .push(VarWalkVisit::new(sub_s, sub_e, seg_id, is_fwd));
            }

            col = end_col;
        }
        Ok(())
    }

    /// Finalize the variation graph by emitting deduplicated `L` links and `PanSN` `W` walks.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if writing links or walks fails.
    pub fn finish(mut self) -> std::io::Result<VariationGfaStats> {
        // Pack oriented link endpoints `(seg_id << 1) | is_fwd` into `(u64, u64)` (16 B/entry).
        let mut emitted_links: HashSet<(u64, u64)> = HashSet::new();
        let mut step_buf = String::new();
        for ((species, contig), mut visits) in self.walks {
            visits.sort_unstable_by_key(|v| (v.fwd_start, v.fwd_end, v.seg_dir));
            let sr_rank = if let Some(&r) = self.species_ranks.get(&species) {
                r
            } else {
                let next = self.species_ranks.len();
                self.species_ranks.insert(species.clone(), next);
                next
            };
            if let Some(src_seq) = self.source_seqs.get(&(species.clone(), contig.clone()))
                && visits.len() >= 2
            {
                let mut spliced = Vec::with_capacity(visits.len() + 8);
                spliced.push(visits[0]);
                for &v in &visits[1..] {
                    let prev_end = spliced[spliced.len() - 1].fwd_end;
                    if v.fwd_start > prev_end
                        && v.fwd_start - prev_end <= MAX_UNALIGNED_INSERTION_SPAN_BP
                        && v.fwd_start <= src_seq.len()
                    {
                        let gap_len = v.fwd_start - prev_end;
                        let seg_id = self.next_seg_id;
                        self.next_seg_id += 1;
                        self.seq_buf.clear();
                        self.seq_buf.extend(
                            src_seq[prev_end..v.fwd_start]
                                .iter()
                                .map(u8::to_ascii_uppercase),
                        );
                        let seq_str = std::str::from_utf8(&self.seq_buf).map_err(|e| {
                            std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                format!("non-UTF-8 spliced insertion segment {seg_id}: {e}"),
                            )
                        })?;
                        writeln!(
                            self.out,
                            "S\t{seg_id}\t{seq_str}\tLN:i:{gap_len}\tSP:i:1\tSN:Z:{species}.{contig}\tSO:i:{prev_end}\tSR:i:{sr_rank}",
                        )?;
                        self.stats.segments_emitted += 1;
                        self.stats.total_segment_bp += gap_len;
                        spliced.push(VarWalkVisit::new(prev_end, v.fwd_start, seg_id, true));
                    }
                    spliced.push(v);
                }
                visits = spliced;
            }
            let (sample, hap, seq_id) = split_pansn_fields(&species, &contig);

            let mut run_start = 0usize;
            while run_start < visits.len() {
                let mut run_end = run_start + 1;
                while run_end < visits.len()
                    && visits[run_end].fwd_start == visits[run_end - 1].fwd_end
                {
                    run_end += 1;
                }
                let run = &visits[run_start..run_end];
                let w_start = run[0].fwd_start;
                let w_end = run[run.len() - 1].fwd_end;

                for win in run.windows(2) {
                    let a = win[0].seg_dir;
                    let b = win[1].seg_dir;
                    let fwd_key = (a, b);
                    let rev_key = (b ^ 1, a ^ 1);
                    if emitted_links.insert(fwd_key.min(rev_key)) {
                        let o1 = if win[0].is_fwd() { '+' } else { '-' };
                        let o2 = if win[1].is_fwd() { '+' } else { '-' };
                        writeln!(
                            self.out,
                            "L\t{}\t{o1}\t{}\t{o2}\t0M",
                            win[0].seg_id(),
                            win[1].seg_id()
                        )?;
                        self.stats.links_emitted += 1;
                    }
                }

                step_buf.clear();
                for v in run {
                    let dir = if v.is_fwd() { '>' } else { '<' };
                    let _ = write!(step_buf, "{dir}{}", v.seg_id());
                }
                writeln!(
                    self.out,
                    "W\t{sample}\t{hap}\t{seq_id}\t{w_start}\t{w_end}\t{step_buf}"
                )?;
                self.stats.walks_emitted += 1;
                run_start = run_end;
            }
        }
        Ok(self.stats)
    }
}

/// Write a reference-free TBA blockset as a native sequence-bearing GFA v1.1 variation graph.
///
/// # Errors
/// Returns an `std::io::Error` if writing to `gfa_out` fails or a block contains non-UTF-8 bytes.
pub fn write_variation_gfa<W: Write>(
    blocks: &[TbaBlock],
    gfa_out: &mut W,
) -> std::io::Result<VariationGfaStats> {
    let mut builder = VariationGfaBuilder::new(gfa_out)?;
    for blk in blocks {
        builder.push_block(blk)?;
    }
    builder.finish()
}

/// Write a reference-free TBA blockset simultaneously to a MAF v1 writer and a native
/// sequence-bearing GFA v1.1 variation graph writer.
///
/// Returns `(total_blocks, multi_species_blocks, variation_gfa_stats)`.
///
/// # Errors
/// Returns an `std::io::Error` if writing to either output fails.
pub fn write_tba_maf_and_variation_gfa<W1: Write, W2: Write>(
    blocks: &[TbaBlock],
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> std::io::Result<(usize, usize, VariationGfaStats)> {
    write_tba_maf_and_variation_gfa_with_seqs(blocks, &[], maf_out, gfa_out)
}

/// Named source genome contigs `(species_name, [(contig_name, dna_bytes)])` for unaligned span splicing.
pub type SourceGenomeContigs = (String, Vec<(String, Vec<u8>)>);

/// Write a reference-free TBA blockset simultaneously to a MAF v1 writer and a native
/// sequence-bearing GFA v1.1 variation graph writer, optionally splicing unaligned/singleton
/// insertion spans (`<= MAX_UNALIGNED_INSERTION_SPAN_BP`) from `source_genomes`.
///
/// # Errors
/// Returns an `std::io::Error` if writing to either output fails.
pub fn write_tba_maf_and_variation_gfa_with_seqs<W1: Write, W2: Write>(
    blocks: &[TbaBlock],
    source_genomes: &[SourceGenomeContigs],
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> std::io::Result<(usize, usize, VariationGfaStats)> {
    maf_out.write_all(b"##maf version=1 scoring=horizon_tba\n")?;
    let mut builder = VariationGfaBuilder::new(gfa_out)?;
    for (sp, ctgs) in source_genomes {
        builder.add_source_contigs(sp, ctgs);
    }
    let mut total_blocks = 0usize;
    let mut multi_species_blocks = 0usize;
    for (idx, blk) in blocks.iter().enumerate() {
        let n_cols = blk.n_cols();
        if blk.rows.len() < 2 || n_cols == 0 {
            continue;
        }
        total_blocks += 1;
        if blk.rows.len() >= 3 {
            multi_species_blocks += 1;
        }
        write_maf_block(idx, blk, maf_out)?;
        builder.push_block(blk)?;
    }
    let stats = builder.finish()?;
    Ok((total_blocks, multi_species_blocks, stats))
}

/// Stream any multi-species MAF v1 reader directly into a native sequence-bearing GFA v1.1
/// variation graph without holding the full MAF blockset in memory.
///
/// # Errors
/// Returns an error if MAF reading or GFA writing fails.
pub fn write_variation_gfa_from_maf_reader<R: BufRead, W: Write>(
    reader: R,
    gfa_out: &mut W,
) -> anyhow::Result<VariationGfaStats> {
    let mut builder = VariationGfaBuilder::new(gfa_out)?;
    for_each_multi_tba_block_reader(reader, |blk| {
        builder.push_block(&blk)?;
        Ok(())
    })?;
    Ok(builder.finish()?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::revcomp;
    use crate::gfa::parse_walk_steps;
    use crate::tba::TbaRow;
    use std::collections::HashMap;

    fn mk_row(species: &str, contig: &str, start: usize, strand: char, seq: &[u8]) -> TbaRow {
        let size = seq.iter().filter(|&&b| b != b'-').count();
        TbaRow {
            species: species.to_owned(),
            contig: contig.to_owned(),
            start,
            size,
            strand,
            src_size: 1000,
            seq: seq.to_vec(),
        }
    }

    /// Reconstruct exact DNA for every `W` walk in a variation GFA and assert `sum(LN) == w_end - w_start`.
    fn reconstruct_walks(gfa: &str) -> HashMap<(String, String, String, usize, usize), Vec<u8>> {
        let mut segs: HashMap<String, Vec<u8>> = HashMap::new();
        let mut out = HashMap::new();
        for line in gfa.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields.first().copied() {
                Some("S") => {
                    let seq = fields[2].as_bytes().to_vec();
                    assert_ne!(seq, b"*", "variation GFA must emit real DNA, not '*'");
                    segs.insert(fields[1].to_owned(), seq);
                }
                Some("W") => {
                    let sample = fields[1].to_owned();
                    let hap = fields[2].to_owned();
                    let ctg = fields[3].to_owned();
                    let w_start: usize = fields[4].parse().expect("w_start");
                    let w_end: usize = fields[5].parse().expect("w_end");
                    let mut seq = Vec::new();
                    for (seg_name, is_fwd) in parse_walk_steps(fields[6]) {
                        let raw = segs.get(seg_name).expect("segment exists");
                        if is_fwd {
                            seq.extend_from_slice(raw);
                        } else {
                            seq.extend_from_slice(&revcomp(raw));
                        }
                    }
                    assert_eq!(
                        seq.len(),
                        w_end - w_start,
                        "walk {sample}#{hap}#{ctg}:{w_start}-{w_end} span mismatch"
                    );
                    out.insert((sample, hap, ctg, w_start, w_end), seq);
                }
                _ => {}
            }
        }
        out
    }

    #[test]
    fn variation_gfa_decomposes_snps_indels_nested_deletions_and_minus_strand() {
        // Block 0:
        //   sp1#1 (+, [10..22)): ACGT A CCGT TT (has SNP A at col 4, nested SNP C at col 6, 2bp ins TT)
        //   sp1#2 (+, [10..22)): ACGT G CTGT TT (has SNP G at col 4, nested SNP T at col 6, 2bp ins TT)
        //   sp2   (-, fwd [50..59)): ungapped block seq = ACGT A ---- TT (length 7 -> fwd_start=50, start=1000-57=943)
        // Wait: ACGT (4) + A (1) + ---- (0) + TT (2) = 7 bp; reverse-complement of forward genomic sequence!
        let sp2_block_seq = b"ACGTA----TT";
        let sp2_ungapped = b"ACGTATT";
        let sp2_fwd_genome = revcomp(sp2_ungapped);
        let sp2_row = TbaRow {
            species: "sp2".to_owned(),
            contig: "chr1".to_owned(),
            start: 1000 - 57, // fwd_interval = [50, 57)
            size: 7,
            strand: '-',
            src_size: 1000,
            seq: sp2_block_seq.to_vec(),
        };
        assert_eq!(sp2_row.fwd_interval(), (50, 57));

        let blocks = vec![TbaBlock {
            score: 50.0,
            rows: vec![
                mk_row("sp1#1", "chr1", 10, '+', b"acgtACCGTtt"),
                mk_row("sp1#2", "chr1", 10, '+', b"ACGTGCTGTTT"),
                sp2_row,
            ],
        }];

        let mut maf_buf = Vec::new();
        let mut gfa_buf = Vec::new();
        let (tot, multi, stats) =
            write_tba_maf_and_variation_gfa(&blocks, &mut maf_buf, &mut gfa_buf).expect("write");
        assert_eq!(tot, 1);
        assert_eq!(multi, 1);
        // Expected compacted segments:
        // col 0..4: ACGT (all 3 agree -> 1 seg)
        // col 4: A (sp1#1, sp2) vs G (sp1#2) -> 2 segs
        // col 5: C (sp1#1, sp1#2; sp2 gap) -> 1 shared seg (compacted!)
        // col 6: C (sp1#1) vs T (sp1#2; sp2 gap) -> 2 segs
        // col 7..9: GT (sp1#1, sp1#2; sp2 gap) -> 1 shared seg (compacted!)
        // col 9..11: TT (all 3 agree -> 1 seg)
        // Total = 8 segments!
        assert_eq!(stats.segments_emitted, 8, "expected 8 compacted segments");
        assert_eq!(stats.walks_emitted, 3);
        assert!(stats.links_emitted >= 8);

        let gfa_str = String::from_utf8(gfa_buf).expect("utf8");
        let walks = reconstruct_walks(&gfa_str);
        assert_eq!(
            walks[&("sp1".to_owned(), "1".to_owned(), "chr1".to_owned(), 10, 21)],
            b"ACGTACCGTTT"
        );
        assert_eq!(
            walks[&("sp1".to_owned(), "2".to_owned(), "chr1".to_owned(), 10, 21)],
            b"ACGTGCTGTTT"
        );
        // Minus-strand walk must reconstruct the EXACT forward genomic sequence!
        assert_eq!(
            walks[&("sp2".to_owned(), "0".to_owned(), "chr1".to_owned(), 50, 57)],
            sp2_fwd_genome
        );

        // Streaming from the emitted MAF must produce byte-identical variation GFA!
        let mut streamed_gfa = Vec::new();
        let streamed_stats =
            write_variation_gfa_from_maf_reader(maf_buf.as_slice(), &mut streamed_gfa)
                .expect("stream from maf");
        assert_eq!(streamed_stats, stats);
        assert_eq!(
            String::from_utf8(streamed_gfa).expect("utf8"),
            gfa_str,
            "streaming MAF -> variation GFA must match direct TbaBlock emission byte-for-byte"
        );
        assert_eq!(std::mem::size_of::<VarWalkVisit>(), 24);
    }

    #[test]
    fn variation_gfa_deduplicates_opposite_strand_canonical_links() {
        // Block 0 (ACGT) immediately adjacent to Block 1 (TTTT) on both sp1 ('+') and sp2 ('-').
        // sp1 walks `>1>2`, while sp2 walks `<2<1` (reverse-complement of the exact same edge).
        // Canonical link deduplication must emit exactly ONE `L` line between segments 1 and 2.
        let blocks = vec![
            TbaBlock {
                score: 20.0,
                rows: vec![
                    mk_row("sp1", "chr1", 0, '+', b"ACGT"),
                    TbaRow {
                        species: "sp2".to_owned(),
                        contig: "chr1".to_owned(),
                        start: 1000 - 14, // fwd [10..14)
                        size: 4,
                        strand: '-',
                        src_size: 1000,
                        seq: b"ACGT".to_vec(),
                    },
                ],
            },
            TbaBlock {
                score: 20.0,
                rows: vec![
                    mk_row("sp1", "chr1", 4, '+', b"TTTT"),
                    TbaRow {
                        species: "sp2".to_owned(),
                        contig: "chr1".to_owned(),
                        start: 1000 - 10, // fwd [6..10), immediately preceding [10..14) on forward sp2!
                        size: 4,
                        strand: '-',
                        src_size: 1000,
                        seq: b"TTTT".to_vec(),
                    },
                ],
            },
        ];
        let mut maf_buf = Vec::new();
        let mut gfa_buf = Vec::new();
        let (_, _, stats) =
            write_tba_maf_and_variation_gfa(&blocks, &mut maf_buf, &mut gfa_buf).expect("write");
        let gfa_str = String::from_utf8(gfa_buf).expect("utf8");
        let l_lines: Vec<&str> = gfa_str.lines().filter(|l| l.starts_with("L\t")).collect();
        assert_eq!(stats.segments_emitted, 2, "gfa:\n{gfa_str}");
        assert_eq!(
            l_lines.len(),
            1,
            "opposite-strand `>1>2` and `<2<1` must deduplicate to 1 canonical L line:\n{gfa_str}"
        );
        assert_eq!(stats.links_emitted, 1);
    }
}
