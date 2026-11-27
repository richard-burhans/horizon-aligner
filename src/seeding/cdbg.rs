//! Colored Compacted de Bruijn Graph (`cDBG`, `TwoPaCo`/`SibeliaZ`-style) for `O(N)`
//! linear-time multi-genome anchor and collinear block extraction.

use crate::seeding::{EisensteinHasher, complement_base};
use std::collections::{BTreeMap, HashMap, HashSet};

/// Occurrence of a colored unitig in one input genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnitigOccurrence {
    pub genome_idx: usize,
    pub contig_idx: usize,
    pub pos: usize,
    pub is_fwd: bool,
}

/// A maximal non-branching colored unitig in a multi-genome compacted de Bruijn graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdbgUnitig {
    pub id: usize,
    pub seq: Vec<u8>,
    /// Sorted distinct genome indices containing this unitig.
    pub genomes: Vec<usize>,
    pub occurrences: Vec<UnitigOccurrence>,
}

/// A multi-genome collinear anchor block extracted directly from shared cDBG unitigs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdbgCollinearBlock {
    /// Unitig IDs in collinear order.
    pub unitig_ids: Vec<usize>,
    /// Participating genomes and their `(genome_idx, contig_idx, start, end, is_fwd)` spans.
    pub spans: Vec<(usize, usize, usize, usize, bool)>,
    pub total_bp: usize,
}

/// Canonical strand classification of a k-mer relative to its reverse complement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum KmerOrient {
    Fwd,
    Rev,
    Palindrome,
}

/// Classify `kmer` against `revcomp(kmer)` (`Fwd`, `Rev`, or `Palindrome` for even-`k` self-rc).
#[inline]
fn classify_kmer_orient(kmer: &[u8]) -> KmerOrient {
    let k = kmer.len();
    for i in 0..k {
        let f = kmer[i].to_ascii_uppercase();
        let r = complement_base(kmer[k - 1 - i]).to_ascii_uppercase();
        if f != r {
            return if f < r {
                KmerOrient::Fwd
            } else {
                KmerOrient::Rev
            };
        }
    }
    KmerOrient::Palindrome
}

#[derive(Debug, Clone, Copy)]
struct KmerOcc {
    genome_idx: usize,
    contig_idx: usize,
    pos: usize,
    orient: KmerOrient,
}

/// Check whether two k-mer occurrence lists `prev` and `next` advance synchronously by `+1`
/// (on same-strand occurrences) / `-1` (on opposite-strand occurrences) across every single
/// occurrence, handling even-`k` palindromes cleanly without spurious unitig breaks.
fn occurrences_extend_synchronously(
    prev: &[KmerOcc],
    prev_anchor_orient: KmerOrient,
    next: &[KmerOcc],
    next_anchor_orient: KmerOrient,
) -> bool {
    if prev.len() != next.len() {
        return false;
    }
    for (p, n) in prev.iter().zip(next.iter()) {
        if p.genome_idx != n.genome_idx || p.contig_idx != n.contig_idx {
            return false;
        }
        let step_fwd = n.pos == p.pos.saturating_add(1);
        let step_rev = p.pos == n.pos.saturating_add(1);
        if !step_fwd && !step_rev {
            return false;
        }
        let prev_same_ok = prev_anchor_orient == KmerOrient::Palindrome
            || p.orient == KmerOrient::Palindrome
            || p.orient == prev_anchor_orient;
        let prev_opp_ok = prev_anchor_orient == KmerOrient::Palindrome
            || p.orient == KmerOrient::Palindrome
            || p.orient != prev_anchor_orient;
        let next_same_ok = next_anchor_orient == KmerOrient::Palindrome
            || n.orient == KmerOrient::Palindrome
            || n.orient == next_anchor_orient;
        let next_opp_ok = next_anchor_orient == KmerOrient::Palindrome
            || n.orient == KmerOrient::Palindrome
            || n.orient != next_anchor_orient;
        let ok = if step_fwd {
            prev_same_ok && next_same_ok
        } else {
            prev_opp_ok && next_opp_ok
        };
        if !ok {
            return false;
        }
    }
    true
}

/// Build a colored compacted de Bruijn graph over `genomes` (`genomes[g]` = slice of contig
/// sequences) for k-mer size `k` in a single `O(sum |G_i|)` pass using canonical Eisenstein hashes.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn build_colored_cdbg(
    genomes: &[&[&[u8]]],
    k: usize,
    max_copy_per_genome: usize,
) -> Vec<CdbgUnitig> {
    if k == 0 || genomes.is_empty() || max_copy_per_genome == 0 {
        return Vec::new();
    }
    let hasher = EisensteinHasher::default();
    let mut kmer_occs: HashMap<u64, Vec<KmerOcc>> = HashMap::new();
    let mut genome_kmer_streams: Vec<Vec<Vec<(usize, u64)>>> = Vec::with_capacity(genomes.len());

    // Single linear O(sum |G_i|) pass over all genomes to index canonical Eisenstein k-mers.
    for (g_idx, contigs) in genomes.iter().enumerate() {
        let mut ctg_streams = Vec::with_capacity(contigs.len());
        for (c_idx, &seq) in contigs.iter().enumerate() {
            let stream = hasher.rolling_kmer_hashes(seq, k);
            for &(pos, h) in &stream {
                let orient = classify_kmer_orient(&seq[pos..pos + k]);
                kmer_occs.entry(h).or_default().push(KmerOcc {
                    genome_idx: g_idx,
                    contig_idx: c_idx,
                    pos,
                    orient,
                });
            }
            ctg_streams.push(stream);
        }
        genome_kmer_streams.push(ctg_streams);
    }

    // Filter out high-copy repeat k-mers exceeding `max_copy_per_genome` in any single genome,
    // and sort occurrences canonically by `(genome_idx, contig_idx, pos)`.
    kmer_occs.retain(|_, occs| {
        occs.sort_by_key(|o| (o.genome_idx, o.contig_idx, o.pos));
        let mut cur_g = usize::MAX;
        let mut count = 0usize;
        for o in occs.iter() {
            if o.genome_idx == cur_g {
                count += 1;
                if count > max_copy_per_genome {
                    return false;
                }
            } else {
                cur_g = o.genome_idx;
                count = 1;
            }
        }
        true
    });

    // Compact maximal non-branching colored runs across all genomes without duplicating unitigs.
    let mut visited: HashSet<u64> = HashSet::with_capacity(kmer_occs.len());
    let mut unitigs: Vec<CdbgUnitig> = Vec::new();

    for (g_idx, contigs) in genomes.iter().enumerate() {
        for (c_idx, &seq) in contigs.iter().enumerate() {
            let stream = &genome_kmer_streams[g_idx][c_idx];
            let mut i = 0usize;
            while i < stream.len() {
                let (start_pos, start_h) = stream[i];
                if visited.contains(&start_h) {
                    i += 1;
                    continue;
                }
                let Some(start_occs) = kmer_occs.get(&start_h) else {
                    i += 1;
                    continue;
                };
                let start_anchor_orient = start_occs
                    .iter()
                    .find(|o| o.genome_idx == g_idx && o.contig_idx == c_idx && o.pos == start_pos)
                    .map_or(KmerOrient::Fwd, |o| o.orient);

                visited.insert(start_h);
                let mut end_i = i;
                let mut cur_occs = start_occs;
                let mut cur_anchor_orient = start_anchor_orient;
                while end_i + 1 < stream.len() {
                    let (next_pos, next_h) = stream[end_i + 1];
                    if next_pos != stream[end_i].0 + 1 || visited.contains(&next_h) {
                        break;
                    }
                    let Some(next_occs) = kmer_occs.get(&next_h) else {
                        break;
                    };
                    let next_anchor_orient = next_occs
                        .iter()
                        .find(|o| {
                            o.genome_idx == g_idx && o.contig_idx == c_idx && o.pos == next_pos
                        })
                        .map_or(KmerOrient::Fwd, |o| o.orient);
                    if !occurrences_extend_synchronously(
                        cur_occs,
                        cur_anchor_orient,
                        next_occs,
                        next_anchor_orient,
                    ) {
                        break;
                    }
                    visited.insert(next_h);
                    cur_occs = next_occs;
                    cur_anchor_orient = next_anchor_orient;
                    end_i += 1;
                }

                let run_len = (end_i - i) + k;
                let unitig_seq = seq[start_pos..start_pos + run_len].to_vec();
                let unitig_occs: Vec<UnitigOccurrence> = start_occs
                    .iter()
                    .zip(cur_occs.iter())
                    .map(|(s_occ, e_occ)| {
                        let same_orient = if end_i > i {
                            e_occ.pos > s_occ.pos
                        } else {
                            start_anchor_orient == KmerOrient::Palindrome
                                || s_occ.orient == KmerOrient::Palindrome
                                || s_occ.orient == start_anchor_orient
                        };
                        let u_start = if same_orient { s_occ.pos } else { e_occ.pos };
                        UnitigOccurrence {
                            genome_idx: s_occ.genome_idx,
                            contig_idx: s_occ.contig_idx,
                            pos: u_start,
                            is_fwd: same_orient,
                        }
                    })
                    .collect();

                let mut distinct_genomes: Vec<usize> =
                    unitig_occs.iter().map(|o| o.genome_idx).collect();
                distinct_genomes.sort_unstable();
                distinct_genomes.dedup();

                unitigs.push(CdbgUnitig {
                    id: unitigs.len(),
                    seq: unitig_seq,
                    genomes: distinct_genomes,
                    occurrences: unitig_occs,
                });
                i = end_i + 1;
            }
        }
    }

    unitigs
}

/// Extract multi-genome collinear blocks supported by at least `min_genomes` from cDBG unitigs.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn extract_cdbg_collinear_blocks(
    unitigs: &[CdbgUnitig],
    min_genomes: usize,
    max_gap_bp: usize,
) -> Vec<CdbgCollinearBlock> {
    if min_genomes == 0 || unitigs.is_empty() {
        return Vec::new();
    }
    // Keep single-copy-per-genome shared unitigs with >= min_genomes support.
    let mut shared: Vec<&CdbgUnitig> = unitigs
        .iter()
        .filter(|u| u.genomes.len() >= min_genomes && u.occurrences.len() == u.genomes.len())
        .collect();
    if shared.is_empty() {
        return Vec::new();
    }

    // Group deterministically by sorted genome-set signature so collinear chaining stays coherent.
    let mut by_sig: BTreeMap<Vec<usize>, Vec<&CdbgUnitig>> = BTreeMap::new();
    for u in shared.drain(..) {
        by_sig.entry(u.genomes.clone()).or_default().push(u);
    }

    let mut blocks = Vec::new();
    for (sig, mut group) in by_sig {
        let anchor_g = sig[0];
        // Sort by anchor genome `(contig_idx, pos, id)`.
        group.sort_by_key(|u| {
            let o = &u.occurrences[0];
            debug_assert_eq!(o.genome_idx, anchor_g);
            (o.contig_idx, o.pos, u.id)
        });

        let mut i = 0usize;
        while i < group.len() {
            let first = group[i];
            let mut chain_ids = vec![first.id];
            let mut total_bp = first.seq.len();
            // Track `(genome_idx, contig_idx, start, end, rel_fwd)` for every participating genome.
            let first_anchor_fwd = first.occurrences[0].is_fwd;
            let mut cur_spans: Vec<(usize, usize, usize, usize, bool)> = first
                .occurrences
                .iter()
                .map(|o| {
                    (
                        o.genome_idx,
                        o.contig_idx,
                        o.pos,
                        o.pos + first.seq.len(),
                        o.is_fwd == first_anchor_fwd,
                    )
                })
                .collect();

            let mut j = i + 1;
            while j < group.len() {
                let nxt = group[j];
                let nxt_anchor_fwd = nxt.occurrences[0].is_fwd;
                let mut can_chain = true;
                let mut updated = cur_spans.clone();
                let mut added_anchor_bp = 0usize;
                for (slot, (sp, no)) in cur_spans.iter().zip(nxt.occurrences.iter()).enumerate() {
                    let rel_fwd = no.is_fwd == nxt_anchor_fwd;
                    if sp.1 != no.contig_idx || sp.4 != rel_fwd {
                        can_chain = false;
                        break;
                    }
                    let n_start = no.pos;
                    let n_end = no.pos + nxt.seq.len();
                    if rel_fwd {
                        // Permit ordered (k-1)-overlapping unitigs separated by a filtered k-mer
                        // or color branch, without double-counting shared overlap bases.
                        if n_start <= sp.2
                            || n_end <= sp.3
                            || n_start.saturating_sub(sp.3) > max_gap_bp
                        {
                            can_chain = false;
                            break;
                        }
                        if slot == 0 {
                            added_anchor_bp = n_end.saturating_sub(sp.3).min(nxt.seq.len());
                        }
                        updated[slot].3 = n_end;
                    } else {
                        if n_start >= sp.2
                            || n_end >= sp.3
                            || sp.2.saturating_sub(n_end) > max_gap_bp
                        {
                            can_chain = false;
                            break;
                        }
                        if slot == 0 {
                            added_anchor_bp = sp.2.saturating_sub(n_start).min(nxt.seq.len());
                        }
                        updated[slot].2 = n_start;
                    }
                }
                if can_chain {
                    chain_ids.push(nxt.id);
                    total_bp += added_anchor_bp;
                    cur_spans = updated;
                    j += 1;
                } else {
                    break;
                }
            }
            blocks.push(CdbgCollinearBlock {
                unitig_ids: chain_ids,
                spans: cur_spans,
                total_bp,
            });
            i = j;
        }
    }

    blocks.sort_by(|a, b| {
        b.total_bp
            .cmp(&a.total_bp)
            .then_with(|| a.spans.cmp(&b.spans))
            .then_with(|| a.unitig_ids.cmp(&b.unitig_ids))
    });
    blocks
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
    fn colored_cdbg_compacts_shared_unitigs_across_4_genomes_including_revcomp_and_snps() {
        let base = synth_dna(0xCDBC_1234, 800);
        let mut g1 = base.clone();
        g1[400] = if g1[400] == b'A' { b'G' } else { b'A' };
        let mut g2 = base.clone();
        g2[400] = g1[400];
        let g3 = revcomp(&base);

        let c0 = [base.as_slice()];
        let c1 = [g1.as_slice()];
        let c2 = [g2.as_slice()];
        let c3 = [g3.as_slice()];
        let genomes: [&[&[u8]]; 4] = [&c0, &c1, &c2, &c3];

        let unitigs = build_colored_cdbg(&genomes, 21, 2);
        let shared_4: Vec<&CdbgUnitig> = unitigs.iter().filter(|u| u.genomes.len() == 4).collect();
        assert!(
            shared_4.len() >= 2,
            "expected >= 2 4-genome shared compacted unitigs flanking SNP at 400, got {}",
            shared_4.len()
        );
        let max_len = shared_4.iter().map(|u| u.seq.len()).max().unwrap_or(0);
        assert!(
            max_len >= 350,
            "expected compacted unitigs >= 350 bp (left/right of SNP), got max {max_len} bp"
        );

        let blocks = extract_cdbg_collinear_blocks(&unitigs, 4, 100);
        assert!(
            !blocks.is_empty(),
            "expected >= 1 4-genome collinear block chaining left & right unitigs across SNP"
        );
        assert!(
            blocks[0].total_bp >= 700,
            "expected chained collinear block >= 700 bp across 800bp locus, got {}",
            blocks[0].total_bp
        );
        assert_eq!(blocks[0].spans.len(), 4);
    }

    #[test]
    fn extract_cdbg_chains_k_minus_1_overlapping_unitigs_across_filtered_repeat_kmer() {
        let left = synth_dna(0x1111_2222, 200);
        let right = synth_dna(0x3333_4444, 200);
        let rep = b"ACGTACGTACGTACGTACGTACGTACGTACG"; // 31 bp repeated twice in genome 0
        let mut g0 = Vec::new();
        g0.extend_from_slice(&left);
        g0.extend_from_slice(rep);
        g0.extend_from_slice(&right);
        g0.extend_from_slice(rep); // pushes `rep` k-mers over max_copy_per_genome=1!
        let mut g1 = Vec::new();
        g1.extend_from_slice(&left);
        g1.extend_from_slice(rep);
        g1.extend_from_slice(&right);

        let unitigs = build_colored_cdbg(&[&[g0.as_slice()], &[g1.as_slice()]], 21, 1);
        let blocks = extract_cdbg_collinear_blocks(&unitigs, 2, 50);
        assert!(
            !blocks.is_empty(),
            "expected collinear block across filtered repeat k-mer"
        );
        assert!(
            blocks[0].unitig_ids.len() >= 2,
            "expected (k-1)-overlapping flanking unitigs to chain into 1 block, got {:?}",
            blocks[0].unitig_ids
        );
        assert!(
            blocks[0].total_bp <= g1.len(),
            "total_bp ({}) must not double-count overlap beyond sequence length ({})",
            blocks[0].total_bp,
            g1.len()
        );
    }
}
