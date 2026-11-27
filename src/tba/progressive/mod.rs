//! Reference-free Progressive Guide-Tree Ancestral Alignment (Figure 1D / Progressive Cactus-style).
//!
//! Walks a binary [`GuideTree`] bottom-up: at each internal node, selects up to `max_outgroups`
//! phylogenetically nearest non-descendant leaf outgroups, runs a small reference-free 2–4 genome
//! subproblem alignment, reconstructs an outgroup-polarized Fitch parsimony ancestral genome
//! (`Anc0`, `Anc1`, ...), threads child clade blocksets upward through the ancestor, and builds a
//! native pure-Rust hierarchical [`HalAlignment`].

use super::merge::{run_tba_tree_merge, tba_merge_two_blocksets};
use super::model::{PairwiseTbaBlock, TbaBlock, TbaRow, canonical_pair_key};
use super::tree::GuideTree;
mod ancestral;
mod hal_build;

use crate::hal::HalAlignment;
use crate::mdl::{AlignConfig, align_contig_pair_runs};
use crate::seeding::{
    ContigIndex, adaptive_top_k_for_fragmentation, prepare_contig_index_owned, select_top_k_pairs,
    shared_seed_count,
};
pub use ancestral::{
    AncestralBranchLengths, reconstruct_ancestral_contigs, reconstruct_ancestral_contigs_transducer,
};
use ancestral::{
    estimate_subproblem_branches, reconstruct_ancestral_contigs_inner,
    reconstruct_ancestral_contigs_transducer_inner,
};
use anyhow::{Context, Result};
use hal_build::build_progressive_hal;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;

/// Configuration for [`run_progressive_ancestral_tba`].
#[derive(Debug, Clone, Copy)]
pub struct ProgressiveAncestralConfig {
    /// Maximum number of nearest non-descendant outgroup leaves per internal subproblem (default: 3).
    pub max_outgroups: usize,
    /// Per-contig top-k pair pruning (`0` = all pairs passing seed threshold).
    pub top_k: usize,
    /// Skip softmasked lowercase bases during seeding.
    pub respect_softmask: bool,
    /// Use branch-length-aware HKY85 + 2-state CTMC indel-polarity ML reconstruction
    /// instead of outgroup-polarized Fitch parsimony.
    pub use_transducer_ml: bool,
}

impl Default for ProgressiveAncestralConfig {
    fn default() -> Self {
        Self {
            // Armstrong et al. (2020) Progressive Cactus default: 3 nearest outgroups.
            max_outgroups: 3,
            top_k: 0,
            respect_softmask: false,
            use_transducer_ml: false,
        }
    }
}

/// Summary statistics for one reconstructed internal ancestral node.
#[derive(Debug, Clone, Default)]
pub struct AncestralNodeStats {
    /// Name of the internal ancestor (`Anc0`, `Anc1`, ...).
    pub name: String,
    /// Left child genome name.
    pub left_child: String,
    /// Right child genome name.
    pub right_child: String,
    /// Outgroup leaves used to polarize indels and substitutions.
    pub outgroups: Vec<String>,
    /// Number of reconstructed ancestral contigs.
    pub n_anc_contigs: usize,
    /// Total reconstructed ancestral bases.
    pub anc_bp: usize,
    /// Bases supported by both child ingroups.
    pub shared_both_children_bp: usize,
    /// Bases present in only one child but supported by >=1 outgroup (lineage-specific deletion).
    pub outgroup_polarized_deletion_bp: usize,
    /// Bases present in only one child with zero outgroup support (lineage-specific insertion).
    pub lineage_insertion_excluded_bp: usize,
    /// Estimated JC69 branch length to `left_child` (populated for ML and branch-length diagnostics).
    pub left_t: Option<f64>,
    /// Estimated JC69 branch length to `right_child`.
    pub right_t: Option<f64>,
    /// Estimated JC69 branch length to the outgroup clade.
    pub outgroup_t: Option<f64>,
}

/// Output of reconstructing a single internal ancestral genome from its subproblem TBA blocks.
#[derive(Debug, Clone)]
pub struct AncestralNodeResult {
    /// Reconstructed ancestral contigs `(contig_name, uppercase_dna)`.
    pub anc_contigs: Vec<(String, Vec<u8>)>,
    /// Exact pairwise blocks linking `Anc_k` -> `left_child`.
    pub anc_to_left: Vec<PairwiseTbaBlock>,
    /// Exact pairwise blocks linking `Anc_k` -> `right_child`.
    pub anc_to_right: Vec<PairwiseTbaBlock>,
    /// Polarization and size telemetry for this node.
    pub stats: AncestralNodeStats,
}

/// Named contig sequence list `(contig_name, dna_bytes)`.
pub type GenomeContigs = Vec<(String, Vec<u8>)>;

/// Reconstructed ancestral genome record `(ancestor_name, contigs, stats)`.
pub type AncestralGenomeEntry = (String, GenomeContigs, AncestralNodeStats);

/// Complete output of a Figure 1D progressive ancestral alignment run.
#[derive(Debug, Clone)]
pub struct ProgressiveAncestralResult {
    /// Reference-free multi-species TBA blocks over all extant leaf species.
    pub merged_blocks: Vec<TbaBlock>,
    /// Reconstructed internal ancestral genomes `(name, contigs, stats)` in post-order.
    pub ancestors: Vec<AncestralGenomeEntry>,
    /// Native pure-Rust hierarchical HAL alignment over all leaves and internal ancestors.
    pub hal: HalAlignment,
    /// Total pairwise subproblem alignments executed (`O(N)` across the tree).
    pub subproblem_pairs_aligned: usize,
}

fn collect_leaf_paths<'a>(
    node: &'a GuideTree,
    prefix: &mut Vec<bool>,
    out: &mut HashMap<&'a str, Vec<bool>>,
) {
    match node {
        GuideTree::Leaf(sp) => {
            out.insert(sp.as_str(), prefix.clone());
        }
        GuideTree::Internal { left, right } => {
            prefix.push(false);
            collect_leaf_paths(left, prefix, out);
            prefix.pop();
            prefix.push(true);
            collect_leaf_paths(right, prefix, out);
            prefix.pop();
        }
    }
}

fn path_distance(a: &[bool], b: &[bool]) -> usize {
    let lcp = a.iter().zip(b.iter()).take_while(|(x, y)| x == y).count();
    (a.len() - lcp) + (b.len() - lcp)
}

/// Select up to `max_outgroups` phylogenetically nearest non-descendant leaf outgroups for an
/// internal clade `node_leaves` within `root`, breaking ties deterministically by species name.
#[must_use]
pub fn select_node_outgroups<S: std::hash::BuildHasher>(
    root: &GuideTree,
    node_leaves: &HashSet<&str, S>,
    max_outgroups: usize,
) -> Vec<String> {
    if max_outgroups == 0 {
        return Vec::new();
    }
    let mut paths = HashMap::new();
    let mut prefix = Vec::new();
    collect_leaf_paths(root, &mut prefix, &mut paths);

    let mut candidates: Vec<(usize, &str)> = Vec::new();
    for (&leaf, leaf_path) in &paths {
        if node_leaves.contains(leaf) {
            continue;
        }
        let min_dist = node_leaves
            .iter()
            .filter_map(|&desc| paths.get(desc).map(|p| path_distance(leaf_path, p)))
            .min()
            .unwrap_or(usize::MAX);
        candidates.push((min_dist, leaf));
    }
    candidates.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(b.1)));
    candidates
        .into_iter()
        .take(max_outgroups)
        .map(|(_, sp)| sp.to_owned())
        .collect()
}

pub(super) const fn fitch_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

pub(super) const fn pick_base_from_mask(mask: u8, fallback: u8) -> u8 {
    if mask & 1 != 0 {
        b'A'
    } else if mask & 2 != 0 {
        b'C'
    } else if mask & 4 != 0 {
        b'G'
    } else if mask & 8 != 0 {
        b'T'
    } else {
        fallback.to_ascii_uppercase()
    }
}

fn align_subproblem_pair(
    sp_a: &str,
    idx_a: &[(String, ContigIndex)],
    sp_b: &str,
    idx_b: &[(String, ContigIndex)],
    cfg: &AlignConfig,
    top_k: usize,
    deep_ancestral: bool,
) -> Result<Vec<PairwiseTbaBlock>> {
    // Internal reconstructed ancestors fragment into multiple contigs and accumulate descendant
    // divergence; when the contig-pair matrix is small, let `align_contig_pair_runs` run its
    // transition (`LASTZ12`) and translated (`AA5`) seeders even when exact 11-of-15 primary
    // seeds drop to zero.
    let small_pair_matrix = idx_a.len().saturating_mul(idx_b.len()) <= 64;
    let min_shared_seeds = if deep_ancestral && small_pair_matrix {
        0usize
    } else if deep_ancestral {
        2usize
    } else {
        8usize
    };
    let pair_cfg = cfg.clone();
    let eff_top_k = adaptive_top_k_for_fragmentation(top_k, idx_a, idx_b);
    let top_mask = (eff_top_k > 0 && min_shared_seeds > 0)
        .then(|| select_top_k_pairs(idx_a, idx_b, eff_top_k, min_shared_seeds));
    let mut contig_pairs = Vec::new();
    for xi in 0..idx_a.len() {
        for yi in 0..idx_b.len() {
            let keep = if min_shared_seeds == 0 {
                true
            } else if let Some(mask) = &top_mask {
                mask[xi][yi]
            } else {
                shared_seed_count(&idx_a[xi].1, &idx_b[yi].1) >= min_shared_seeds
            };
            if keep {
                contig_pairs.push((xi, yi));
            }
        }
    }
    let chunks: Vec<Vec<PairwiseTbaBlock>> = contig_pairs
        .par_iter()
        .map(|&(xi, yi)| {
            let (nx, ix) = &idx_a[xi];
            let (ny, iy) = &idx_b[yi];
            let aln = align_contig_pair_runs(ix, iy, &pair_cfg)?;
            if aln.aligned_runs.is_empty() {
                return Ok(Vec::new());
            }
            Ok(PairwiseTbaBlock::from_run_alignment(
                sp_a,
                nx,
                ix.seq(),
                sp_b,
                ny,
                iy.seq(),
                &aln,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut out = Vec::new();
    for chunk in chunks {
        out.extend(chunk);
    }
    Ok(out)
}

fn build_subproblem_tree(left: &str, right: &str, outgroups: &[String]) -> GuideTree {
    let mut cur = GuideTree::Internal {
        left: Box::new(GuideTree::Leaf(left.to_owned())),
        right: Box::new(GuideTree::Leaf(right.to_owned())),
    };
    for og in outgroups {
        cur = GuideTree::Internal {
            left: Box::new(cur),
            right: Box::new(GuideTree::Leaf(og.clone())),
        };
    }
    cur
}

struct NodeProduct<'a> {
    name: String,
    contigs: std::borrow::Cow<'a, [(String, Vec<u8>)]>,
    clade_blocks: Vec<TbaBlock>,
    clade_species: HashSet<String>,
    newick: String,
}

pub(super) struct HalEdgeRecord {
    pub(super) parent: String,
    pub(super) left_child: String,
    pub(super) right_child: String,
    pub(super) parent_to_left: Vec<PairwiseTbaBlock>,
    pub(super) parent_to_right: Vec<PairwiseTbaBlock>,
}

struct ProgState<'a, S> {
    root: &'a GuideTree,
    leaf_genomes: &'a HashMap<String, Vec<(String, Vec<u8>)>, S>,
    cfg: &'a ProgressiveAncestralConfig,
    align_cfg: AlignConfig,
    anc_counter: usize,
    subproblem_pairs_aligned: usize,
    ancestors: Vec<AncestralGenomeEntry>,
    used_names: HashSet<String>,
    hal_edges: Vec<HalEdgeRecord>,
}

#[allow(clippy::too_many_lines)]
fn run_node_recursive<'a, S: std::hash::BuildHasher>(
    node: &GuideTree,
    st: &mut ProgState<'a, S>,
    is_root: bool,
) -> Result<NodeProduct<'a>> {
    match node {
        GuideTree::Leaf(sp) => {
            let contigs = st
                .leaf_genomes
                .get(sp)
                .with_context(|| format!("missing FASTA contigs for leaf species `{sp}`"))?;
            st.used_names.insert(sp.clone());
            let mut clade_species = HashSet::new();
            clade_species.insert(sp.clone());
            Ok(NodeProduct {
                name: sp.clone(),
                contigs: std::borrow::Cow::Borrowed(contigs.as_slice()),
                clade_blocks: Vec::new(),
                clade_species,
                newick: sp.clone(),
            })
        }
        GuideTree::Internal { left, right } => {
            let left_prod = run_node_recursive(left, st, false)?;
            let right_prod = run_node_recursive(right, st, false)?;

            let mut anc_name = format!("Anc{}", st.anc_counter);
            st.anc_counter += 1;
            while st.leaf_genomes.contains_key(&anc_name) || st.used_names.contains(&anc_name) {
                anc_name = format!("Anc{}", st.anc_counter);
                st.anc_counter += 1;
            }
            st.used_names.insert(anc_name.clone());

            let desc_leaves: HashSet<&str> = node.leaves().into_iter().collect();
            let outgroups = select_node_outgroups(st.root, &desc_leaves, st.cfg.max_outgroups);

            // Index ONLY the 2 children + up to max_outgroups leaves for this subproblem.
            let mut sub_genomes: Vec<(String, Vec<(String, ContigIndex)>)> =
                Vec::with_capacity(2 + outgroups.len());
            let index_one = |name: &str,
                             ctgs: &[(String, Vec<u8>)]|
             -> Result<(String, Vec<(String, ContigIndex)>)> {
                let idxs = ctgs
                    .par_iter()
                    .map(|(cn, seq)| {
                        Ok((
                            cn.clone(),
                            prepare_contig_index_owned(seq.clone(), st.cfg.respect_softmask)?,
                        ))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok((name.to_owned(), idxs))
            };
            sub_genomes.push(index_one(&left_prod.name, &left_prod.contigs)?);
            sub_genomes.push(index_one(&right_prod.name, &right_prod.contigs)?);
            for og in &outgroups {
                let og_ctgs = st
                    .leaf_genomes
                    .get(og)
                    .with_context(|| format!("missing outgroup FASTA `{og}`"))?;
                sub_genomes.push(index_one(og, og_ctgs)?);
            }

            // Align only ingroup-ingroup (left x right) and ingroup-outgroup pairs (skipping
            // unnecessary outgroup x outgroup pairs, following Armstrong et al. 2020).
            let mut sub_pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
            let left_is_anc = left_prod.clade_species.len() > 1;
            let right_is_anc = right_prod.clade_species.len() > 1;
            for i in 0..2.min(sub_genomes.len()) {
                for j in (i + 1)..sub_genomes.len() {
                    let deep_ancestral = (i == 0 && left_is_anc)
                        || (i == 1 && right_is_anc)
                        || (j == 1 && right_is_anc);
                    let blks = align_subproblem_pair(
                        &sub_genomes[i].0,
                        &sub_genomes[i].1,
                        &sub_genomes[j].0,
                        &sub_genomes[j].1,
                        &st.align_cfg,
                        st.cfg.top_k,
                        deep_ancestral,
                    )?;
                    st.subproblem_pairs_aligned += 1;
                    sub_pairs.insert(
                        canonical_pair_key(&sub_genomes[i].0, &sub_genomes[j].0),
                        blks,
                    );
                }
            }
            drop(sub_genomes);

            let sub_tree = build_subproblem_tree(&left_prod.name, &right_prod.name, &outgroups);
            let sub_blocks = run_tba_tree_merge(&sub_tree, &sub_pairs);

            // On non-root internal nodes (`!is_root`), retain uncovered child spans inside the
            // reconstructed ancestor so lineage-specific deletions in one child (or regions missed
            // by the nearest outgroups) remain alignable to distant clades higher up the tree.
            let union_children =
                (!is_root).then_some((&left_prod.contigs[..], &right_prod.contigs[..]));
            let recon = if st.cfg.use_transducer_ml {
                let branches = estimate_subproblem_branches(
                    &left_prod.name,
                    &right_prod.name,
                    &outgroups,
                    &sub_blocks,
                );
                reconstruct_ancestral_contigs_transducer_inner(
                    &anc_name,
                    &left_prod.name,
                    &right_prod.name,
                    &outgroups,
                    &sub_blocks,
                    branches,
                    union_children,
                )
            } else {
                reconstruct_ancestral_contigs_inner(
                    &anc_name,
                    &left_prod.name,
                    &right_prod.name,
                    &outgroups,
                    &sub_blocks,
                    union_children,
                )
            };

            // Thread left clade onto Anc_k, then thread right clade onto the Anc_k-anchored result.
            let left_refs: HashSet<&str> =
                left_prod.clade_species.iter().map(String::as_str).collect();
            let to_left_refs: Vec<&PairwiseTbaBlock> = recon.anc_to_left.iter().collect();
            let left_with_anc =
                tba_merge_two_blocksets(&left_prod.clade_blocks, &[], &left_refs, &to_left_refs);

            let mut acc_species = left_prod.clade_species;
            acc_species.insert(anc_name.clone());
            let acc_refs: HashSet<&str> = acc_species.iter().map(String::as_str).collect();
            let to_right_refs: Vec<&PairwiseTbaBlock> = recon.anc_to_right.iter().collect();
            let merged_clade = tba_merge_two_blocksets(
                &left_with_anc,
                &right_prod.clade_blocks,
                &acc_refs,
                &to_right_refs,
            );
            acc_species.extend(right_prod.clade_species);

            st.ancestors
                .push((anc_name.clone(), recon.anc_contigs.clone(), recon.stats));
            st.hal_edges.push(HalEdgeRecord {
                parent: anc_name.clone(),
                left_child: left_prod.name,
                right_child: right_prod.name,
                parent_to_left: recon.anc_to_left,
                parent_to_right: recon.anc_to_right,
            });

            let newick = format!("({},{}){anc_name}", left_prod.newick, right_prod.newick);
            Ok(NodeProduct {
                name: anc_name,
                contigs: std::borrow::Cow::Owned(recon.anc_contigs),
                clade_blocks: merged_clade,
                clade_species: acc_species,
                newick,
            })
        }
    }
}

/// Execute full bottom-up Figure 1D progressive ancestral alignment along `tree`.
///
/// # Errors
/// Returns an error if any leaf genome in `tree` is missing from `leaf_genomes` or if any
/// subproblem alignment fails.
pub fn run_progressive_ancestral_tba<S: std::hash::BuildHasher>(
    tree: &GuideTree,
    leaf_genomes: &HashMap<String, Vec<(String, Vec<u8>)>, S>,
    cfg: &ProgressiveAncestralConfig,
) -> Result<ProgressiveAncestralResult> {
    let leaf_set: HashSet<String> = tree.leaves().into_iter().map(str::to_owned).collect();

    let mut st = ProgState {
        root: tree,
        leaf_genomes,
        cfg,
        align_cfg: AlignConfig {
            skip_softmasked: cfg.respect_softmask,
            ..AlignConfig::default()
        },
        anc_counter: 0,
        subproblem_pairs_aligned: 0,
        ancestors: Vec::new(),
        used_names: HashSet::new(),
        hal_edges: Vec::new(),
    };

    let root_prod = run_node_recursive(tree, &mut st, true)?;

    // Strip internal Anc* rows from the final extant-species MAF/GFA blockset and drop all-gap columns.
    let mut merged_blocks = Vec::with_capacity(root_prod.clade_blocks.len());
    for blk in root_prod.clade_blocks {
        let leaf_rows: Vec<TbaRow> = blk
            .rows
            .into_iter()
            .filter(|r| leaf_set.contains(&r.species))
            .collect();
        if leaf_rows.len() < 2 {
            continue;
        }
        let n_cols = leaf_rows[0].seq.len();
        let keep_cols: Vec<usize> = (0..n_cols)
            .filter(|&c| leaf_rows.iter().any(|r| r.seq[c] != b'-'))
            .collect();
        if keep_cols.is_empty() {
            continue;
        }
        let compacted: Vec<TbaRow> = leaf_rows
            .into_iter()
            .map(|mut r| {
                r.seq = keep_cols.iter().map(|&c| r.seq[c]).collect();
                r
            })
            .collect();
        merged_blocks.push(TbaBlock {
            score: blk.score,
            rows: compacted,
        });
    }

    let hal = build_progressive_hal(
        &root_prod.name,
        &root_prod.newick,
        leaf_genomes,
        &st.ancestors,
        &st.hal_edges,
    );

    Ok(ProgressiveAncestralResult {
        merged_blocks,
        ancestors: st.ancestors,
        hal,
        subproblem_pairs_aligned: st.subproblem_pairs_aligned,
    })
}

/// Format a tab-separated summary of all reconstructed ancestral nodes (`ancestors_summary.tsv`).
#[must_use]
pub fn format_ancestors_summary_tsv(ancestors: &[AncestralGenomeEntry]) -> String {
    let mut out = String::from(
        "ancestor\tleft_child\tright_child\toutgroups\tn_contigs\tanc_bp\tshared_both_children_bp\toutgroup_polarized_deletion_bp\tlineage_insertion_excluded_bp\tleft_t\tright_t\toutgroup_t\n",
    );
    let fmt_opt = |v: Option<f64>| v.map_or_else(|| "NA".to_owned(), |x| format!("{x:.6}"));
    for (_, _, s) in ancestors {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.name,
            s.left_child,
            s.right_child,
            s.outgroups.join(","),
            s.n_anc_contigs,
            s.anc_bp,
            s.shared_both_children_bp,
            s.outgroup_polarized_deletion_bp,
            s.lineage_insertion_excluded_bp,
            fmt_opt(s.left_t),
            fmt_opt(s.right_t),
            fmt_opt(s.outgroup_t),
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_progressive_hal_preserves_split_reversed_child_runs() {
        // Parent Anc0 has 20 bp. Left child L has 20 bp reversed relative to Anc0; right child R
        // aligns only [0..10) so parent breakpoints split L's reversed run into [0..10) and [10..20).
        // Both pieces of L MUST survive into HAL TopSegments (no monotonic l_used_end drop!).
        let anc_seq = b"ACGTACGTACTTTTGGGGCC".to_vec();
        let l_seq = crate::seeding::revcomp(&anc_seq);
        let r_seq = b"ACGTACGTAC".to_vec();

        let mut all_ctgs: HashMap<String, GenomeContigs> = HashMap::new();
        all_ctgs.insert("Anc0".to_owned(), vec![("a0".to_owned(), anc_seq.clone())]);
        all_ctgs.insert("L".to_owned(), vec![("l0".to_owned(), l_seq.clone())]);
        all_ctgs.insert("R".to_owned(), vec![("r0".to_owned(), r_seq.clone())]);

        let pb_l = PairwiseTbaBlock {
            score: 20.0,
            row_a: TbaRow {
                species: "Anc0".to_owned(),
                contig: "a0".to_owned(),
                start: 0,
                size: 20,
                strand: '+',
                src_size: 20,
                seq: anc_seq.clone(),
            },
            row_b: TbaRow {
                species: "L".to_owned(),
                contig: "l0".to_owned(),
                start: 0,
                size: 20,
                strand: '-',
                src_size: 20,
                seq: anc_seq.clone(),
            },
        };
        let pb_r = PairwiseTbaBlock {
            score: 10.0,
            row_a: TbaRow {
                species: "Anc0".to_owned(),
                contig: "a0".to_owned(),
                start: 0,
                size: 10,
                strand: '+',
                src_size: 20,
                seq: r_seq.clone(),
            },
            row_b: TbaRow {
                species: "R".to_owned(),
                contig: "r0".to_owned(),
                start: 0,
                size: 10,
                strand: '+',
                src_size: 10,
                seq: r_seq,
            },
        };
        let edges = vec![HalEdgeRecord {
            parent: "Anc0".to_owned(),
            left_child: "L".to_owned(),
            right_child: "R".to_owned(),
            parent_to_left: vec![pb_l],
            parent_to_right: vec![pb_r],
        }];
        let hal = build_progressive_hal("Anc0", "(L,R)Anc0", &all_ctgs, &[], &edges);
        let l_genome = hal.genomes.iter().find(|g| g.name == "L").expect("L");
        let mapped_bp: usize = l_genome
            .top_segments
            .iter()
            .filter(|t| t.parent_bottom_idx.is_some())
            .map(|t| t.length)
            .sum();
        assert_eq!(
            mapped_bp, 20,
            "both split halves of reversed child run must survive in HAL (got {mapped_bp}/20)"
        );
    }

    #[test]
    fn build_progressive_hal_never_crosses_contig_boundaries() {
        let mut all_ctgs: HashMap<String, GenomeContigs> = HashMap::new();
        all_ctgs.insert("Anc0".to_owned(), vec![("a0".to_owned(), vec![b'A'; 20]), ("a1".to_owned(), vec![b'C'; 20])]);
        all_ctgs.insert("L".to_owned(), vec![("l0".to_owned(), vec![b'A'; 20]), ("l1".to_owned(), vec![b'C'; 20])]);
        all_ctgs.insert("R".to_owned(), vec![("r0".to_owned(), vec![b'G'; 20])]);
        let mk_pb = |a_ctg: &str, l_ctg: &str, start: usize, base: u8| PairwiseTbaBlock {
            score: 10.0,
            row_a: TbaRow { species: "Anc0".to_owned(), contig: a_ctg.to_owned(), start, size: 10, strand: '+', src_size: 20, seq: vec![base; 10] },
            row_b: TbaRow { species: "L".to_owned(), contig: l_ctg.to_owned(), start, size: 10, strand: '+', src_size: 20, seq: vec![base; 10] },
        };
        let edges = vec![HalEdgeRecord {
            parent: "Anc0".to_owned(),
            left_child: "L".to_owned(),
            right_child: "R".to_owned(),
            parent_to_left: vec![mk_pb("a0", "l0", 0, b'A'), mk_pb("a1", "l1", 10, b'C')],
            parent_to_right: Vec::new(),
        }];
        let hal = build_progressive_hal("Anc0", "(L,R)Anc0", &all_ctgs, &[], &edges);
        for g in &hal.genomes {
            for seq in &g.sequences {
                let s_end = seq.offset + seq.length;
                if !g.top_segments.is_empty() {
                    let sum: usize = g.top_segments.iter().filter(|t| t.start >= seq.offset && t.start < s_end).map(|t| {
                        assert!(t.start + t.length <= s_end, "TopSegment crosses contig end");
                        t.length
                    }).sum();
                    assert_eq!(sum, seq.length, "TopSegment sum mismatch for {}.{}", g.name, seq.name);
                }
                if !g.bottom_segments.is_empty() {
                    let sum: usize = g.bottom_segments.iter().filter(|b| b.start >= seq.offset && b.start < s_end).map(|b| {
                        assert!(b.start + b.length <= s_end, "BottomSegment crosses contig end");
                        b.length
                    }).sum();
                    assert_eq!(sum, seq.length, "BottomSegment sum mismatch for {}.{}", g.name, seq.name);
                }
            }
        }
    }

    #[test]
    fn transducer_ml_excludes_unpolarized_root_insertion_and_zero_child_columns() {
        let mk = |sp: &str, seq: &[u8]| TbaRow {
            species: sp.to_owned(),
            contig: "c1".to_owned(),
            start: 0,
            size: seq.iter().filter(|&&b| b != b'-').count(),
            strand: '+',
            src_size: 100,
            seq: seq.to_vec(),
        };
        // Block 1: unpolarized root (0 outgroups) with equal branches and a 4-bp insertion in L.
        let blk_root = TbaBlock {
            score: 10.0,
            rows: vec![mk("L", b"ACGTAAAAACGT"), mk("R", b"ACGT----ACGT")],
        };
        let branches = AncestralBranchLengths {
            left_t: 0.1,
            right_t: 0.1,
            outgroup_t: 0.25,
            outgroup_stem_t: 0.20,
            outgroup_tip_t: 0.05,
            ins_rate: 0.05,
            del_rate: 0.05,
        };
        let res_root =
            reconstruct_ancestral_contigs_transducer("Anc0", "L", "R", &[], &[blk_root], branches);
        assert_eq!(
            res_root.stats.anc_bp, 8,
            "unpolarized 1-child insertion at root must be excluded"
        );
        assert_eq!(res_root.stats.lineage_insertion_excluded_bp, 4);

        // Block 2: 3 outgroups present where BOTH L and R are gapped -> must NOT fabricate bases.
        let blk_zero_child = TbaBlock {
            score: 10.0,
            rows: vec![
                mk("L", b"ACGT----ACGT"),
                mk("R", b"ACGT----ACGT"),
                mk("O1", b"ACGTTTTTACGT"),
                mk("O2", b"ACGTTTTTACGT"),
                mk("O3", b"ACGTTTTTACGT"),
            ],
        };
        let ogs = vec!["O1".to_owned(), "O2".to_owned(), "O3".to_owned()];
        let res_zc = reconstruct_ancestral_contigs_transducer(
            "Anc1",
            "L",
            "R",
            &ogs,
            &[blk_zero_child],
            branches,
        );
        assert_eq!(
            res_zc.stats.anc_bp, 8,
            "3-outgroup columns with zero ingroup bases must not fabricate ancestral bases"
        );
    }

    #[test]
    fn progressive_ancestral_union_retains_sister_deleted_region_for_distant_clade() {
        // Tree: ((spA, spB), (spC, spD)) with max_outgroups=0 so Anc(spA,spB) has NO outgroup
        // help during reconstruction, yet spA shares a 300 bp middle region with spC that spB
        // completely deleted. Non-root union retention must carry spA's 300 bp span upward so
        // spA aligns to spC at the root!
        let mut state: u64 = 0x1122_3344_5566_7788;
        let bases = *b"ACGT";
        let mut rand_dna = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1);
                    bases[((state >> 32) & 3) as usize]
                })
                .collect()
        };
        let left_flank = rand_dna(400);
        let mid_shared_ac = rand_dna(300);
        let right_flank = rand_dna(400);
        let mut seq_a = Vec::new();
        seq_a.extend_from_slice(&left_flank);
        seq_a.extend_from_slice(&mid_shared_ac);
        seq_a.extend_from_slice(&right_flank);

        let mut seq_b = Vec::new();
        seq_b.extend_from_slice(&left_flank);
        seq_b.extend_from_slice(&right_flank);

        let mut seq_d = Vec::new();
        seq_d.extend_from_slice(&left_flank);
        seq_d.extend_from_slice(&right_flank);

        let tree = crate::tba::parse_newick("((spA,spB),(spC,spD));").expect("tree");
        let mut leaves = HashMap::new();
        leaves.insert("spA".to_owned(), vec![("c1".to_owned(), seq_a.clone())]);
        leaves.insert("spB".to_owned(), vec![("c1".to_owned(), seq_b)]);
        leaves.insert("spC".to_owned(), vec![("c1".to_owned(), seq_a)]);
        leaves.insert("spD".to_owned(), vec![("c1".to_owned(), seq_d)]);

        let cfg = ProgressiveAncestralConfig {
            max_outgroups: 0,
            top_k: 0,
            respect_softmask: false,
            use_transducer_ml: false,
        };
        let res = run_progressive_ancestral_tba(&tree, &leaves, &cfg).expect("progressive");
        // Verify spA and spC align across the middle [400..700) region deleted in spB and spD.
        let mut covered_mid = 0usize;
        for blk in &res.merged_blocks {
            let ra = blk.rows.iter().find(|r| r.species == "spA");
            let rc = blk.rows.iter().find(|r| r.species == "spC");
            if let (Some(a), Some(c)) = (ra, rc) {
                let mut ka = 0usize;
                for (&ba, &bc) in a.seq.iter().zip(&c.seq) {
                    if ba != b'-' {
                        let pos = a.fwd_pos(ka);
                        ka += 1;
                        if bc != b'-' && (400..700).contains(&pos) {
                            covered_mid += 1;
                        }
                    }
                }
            }
        }
        assert!(
            covered_mid >= 250,
            "expected >=250/300 bp of sister-deleted region aligned between spA and spC, got {covered_mid}"
        );
    }
}
