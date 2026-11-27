//! Outgroup-polarized Fitch parsimony and continuous-time HKY85 + reversible indel CTMC
//! ancestral genome reconstruction for progressive TBA.

use super::{AncestralNodeResult, AncestralNodeStats, fitch_mask, pick_base_from_mask};
use crate::coord::{Interval, IntervalSet, Strand};
use crate::formats::MafRow;
use crate::tba::model::{PairwiseTbaBlock, TbaBlock, TbaRow};
use std::collections::{HashMap, HashSet};

struct DraftAncSegment {
    anchor_is_right: bool,
    anchor_contig: String,
    anchor_fwd_start: usize,
    anc_ungapped: Vec<u8>,
    left_meta: Option<(String, usize, char, usize)>,
    left_anc_aln: Vec<u8>,
    left_child_aln: Vec<u8>,
    right_meta: Option<(String, usize, char, usize)>,
    right_anc_aln: Vec<u8>,
    right_child_aln: Vec<u8>,
}

/// Branch lengths `(left_t, right_t, outgroup_t, outgroup_stem_t, outgroup_tip_t)` and 2-state CTMC
/// birth/death rates for [`reconstruct_ancestral_contigs_transducer`].
#[derive(Debug, Clone, Copy)]
pub struct AncestralBranchLengths {
    pub left_t: f64,
    pub right_t: f64,
    pub outgroup_t: f64,
    /// Branch length from the target ancestor to the outgroup clade MRCA (`<= outgroup_t`).
    pub outgroup_stem_t: f64,
    /// Mean tip branch length from the outgroup clade MRCA to each outgroup leaf (`0.0` for `<= 1` outgroup).
    pub outgroup_tip_t: f64,
    /// Per-site insertion rate `λ` (default: `0.05`).
    pub ins_rate: f64,
    /// Per-site deletion rate `μ` (`μ > λ`, default: `0.06`).
    pub del_rate: f64,
}

impl Default for AncestralBranchLengths {
    fn default() -> Self {
        Self {
            left_t: 0.10,
            right_t: 0.10,
            outgroup_t: 0.25,
            outgroup_stem_t: 0.15,
            outgroup_tip_t: 0.10,
            ins_rate: 0.05,
            del_rate: 0.06,
        }
    }
}

/// Two-state continuous-time reversible presence/absence indel kernel over branch length `t`:
/// returns `(p_survive = P(1->1|t), p_insert = P(0->1|t))`. Satisfies exact detailed balance
/// `π_1 (1 - p_surv) = π_0 p_ins` with `π_1 = lam / (lam + mu)`.
fn tkf91_probs(t: f64, lam: f64, mu: f64) -> (f64, f64) {
    let t = t.max(1e-5);
    let mu = mu.max(1e-5);
    let lam = lam.max(1e-5);
    let rate = lam + mu;
    let pi_1 = lam / rate;
    let e_rt = (-rate * t).exp();
    let p_surv = (pi_1 + (1.0 - pi_1) * e_rt).clamp(1e-6, 1.0 - 1e-9);
    let p_ins = (pi_1 * (1.0 - e_rt)).clamp(1e-6, 1.0 - 1e-9);
    (p_surv, p_ins)
}

/// Estimate empirical branch lengths `(left_t, right_t, outgroup_t)` from subproblem TBA blocks
/// using JC69 divergence (falling back to `AncestralBranchLengths::default()` if unaligned).
pub(super) fn estimate_subproblem_branches(
    left_name: &str,
    right_name: &str,
    outgroups: &[String],
    sub_blocks: &[TbaBlock],
) -> AncestralBranchLengths {
    let og_set: HashSet<&str> = outgroups.iter().map(String::as_str).collect();
    let mut pair_lr_cmp = 0usize;
    let mut pair_lr_mm = 0usize;
    let mut left_og_cmp = 0usize;
    let mut left_og_mm = 0usize;
    let mut right_og_cmp = 0usize;
    let mut right_og_mm = 0usize;
    let mut og_og_cmp = 0usize;
    let mut og_og_mm = 0usize;

    let count_pair = |a: &[u8], b: &[u8], cmp: &mut usize, mm: &mut usize| {
        for (&ca, &cb) in a.iter().zip(b) {
            if let (Some(ia), Some(ib)) =
                (crate::bocpd::base_to_idx(ca), crate::bocpd::base_to_idx(cb))
            {
                *cmp += 1;
                if ia != ib {
                    *mm += 1;
                }
            }
        }
    };

    for blk in sub_blocks {
        let lr = blk.rows.iter().find(|r| r.species == left_name);
        let rr = blk.rows.iter().find(|r| r.species == right_name);
        if let (Some(l), Some(r)) = (lr, rr) {
            count_pair(&l.seq, &r.seq, &mut pair_lr_cmp, &mut pair_lr_mm);
        }
        let og_rows: Vec<&TbaRow> = blk
            .rows
            .iter()
            .filter(|r| og_set.contains(r.species.as_str()))
            .collect();
        for (idx, row) in og_rows.iter().enumerate() {
            if let Some(l) = lr {
                count_pair(&l.seq, &row.seq, &mut left_og_cmp, &mut left_og_mm);
            }
            if let Some(r) = rr {
                count_pair(&r.seq, &row.seq, &mut right_og_cmp, &mut right_og_mm);
            }
            for other in &og_rows[(idx + 1)..] {
                count_pair(&row.seq, &other.seq, &mut og_og_cmp, &mut og_og_mm);
            }
        }
    }

    let jc69 = |mm: usize, cmp: usize, fallback: f64| -> f64 {
        if cmp < 10 {
            return fallback;
        }
        #[allow(clippy::cast_precision_loss)]
        let p = (mm as f64 / cmp as f64).clamp(0.0, 0.72);
        let inner = (1.0 - (4.0 / 3.0) * p).max(1e-4);
        (-0.75 * inner.ln()).clamp(0.005, 1.5)
    };

    let dist_lr = jc69(pair_lr_mm, pair_lr_cmp, 0.20);
    let dist_left_og = jc69(left_og_mm, left_og_cmp, dist_lr * 1.5);
    let dist_right_og = jc69(right_og_mm, right_og_cmp, dist_lr * 1.5);

    // Additive 3-taxon branch estimate: t_L = (d_LR + d_LO - d_RO)/2, t_R = (d_LR + d_RO - d_LO)/2.
    let left_t = if left_og_cmp >= 10 && right_og_cmp >= 10 {
        (0.5 * (dist_lr + dist_left_og - dist_right_og)).clamp(0.005, 1.5)
    } else {
        (0.5 * dist_lr).clamp(0.005, 1.5)
    };
    let right_t = if left_og_cmp >= 10 && right_og_cmp >= 10 {
        (0.5 * (dist_lr + dist_right_og - dist_left_og)).clamp(0.005, 1.5)
    } else {
        (0.5 * dist_lr).clamp(0.005, 1.5)
    };
    let outgroup_t = if left_og_cmp >= 10 && right_og_cmp >= 10 {
        (0.5 * (dist_left_og + dist_right_og - dist_lr)).clamp(0.01, 2.0)
    } else {
        0.25
    };
    // Decompose `outgroup_t` into shared stem + tip branches via empirical outgroup-outgroup
    // JC69 distance so correlated sister outgroups peel through their shared clade MRCA
    // without triple-counting, while independent outgroups retain full multi-taxon power.
    let (outgroup_stem_t, outgroup_tip_t) = if outgroups.len() <= 1 {
        (outgroup_t, 0.0)
    } else if og_og_cmp >= 10 {
        let dist_oo = jc69(og_og_mm, og_og_cmp, outgroup_t);
        let tip = (0.5 * dist_oo).clamp(0.005, outgroup_t.max(0.005));
        let stem = (outgroup_t - tip).clamp(0.005, 2.0);
        (stem, tip)
    } else {
        (
            (0.75 * outgroup_t).clamp(0.005, 2.0),
            (0.25 * outgroup_t).clamp(0.005, 1.5),
        )
    };

    AncestralBranchLengths {
        left_t,
        right_t,
        outgroup_t,
        outgroup_stem_t,
        outgroup_tip_t,
        ins_rate: 0.05,
        del_rate: 0.05,
    }
}

/// Choose maximum-likelihood ancestral nucleotide via exact Felsenstein peeling over HKY85 `P(τ)`,
/// marginalizing `>= 2` outgroups through their shared clade MRCA `(p_o_stem, p_o_tip)`.
#[allow(clippy::too_many_arguments)]
fn felsenstein_ml_base(
    lb: Option<u8>,
    rb: Option<u8>,
    ogs: &[u8],
    pi: [f64; 4],
    p_l: &[[f64; 4]; 4],
    p_r: &[[f64; 4]; 4],
    p_o: &[[f64; 4]; 4],
    p_o_stem: &[[f64; 4]; 4],
    p_o_tip: &[[f64; 4]; 4],
    use_clade_peeling: bool,
) -> u8 {
    let bases = *b"ACGT";
    let mut og_idxs = [0usize; 8];
    let mut n_og = 0usize;
    for &ob in ogs {
        if let Some(o_idx) = crate::bocpd::base_to_idx(ob)
            && n_og < og_idxs.len()
        {
            og_idxs[n_og] = o_idx;
            n_og += 1;
        }
    }
    let mut clade_tip_like = [1.0_f64; 4];
    if use_clade_peeling && n_og >= 2 {
        for (g, like) in clade_tip_like.iter_mut().enumerate() {
            let mut prod = 1.0_f64;
            for &o_idx in &og_idxs[..n_og] {
                prod *= p_o_tip[g][o_idx].max(1e-15);
            }
            *like = prod;
        }
    }

    let mut best_b = 0usize;
    let mut best_ll = f64::NEG_INFINITY;
    for anc in 0..4 {
        let mut ll = pi[anc].max(1e-12).ln();
        if let Some(l_idx) = lb.and_then(crate::bocpd::base_to_idx) {
            ll += p_l[anc][l_idx].max(1e-15).ln();
        }
        if let Some(r_idx) = rb.and_then(crate::bocpd::base_to_idx) {
            ll += p_r[anc][r_idx].max(1e-15).ln();
        }
        if n_og == 1 {
            ll += p_o[anc][og_idxs[0]].max(1e-15).ln();
        } else if n_og >= 2 {
            if use_clade_peeling {
                let mut marg = 0.0_f64;
                for g in 0..4 {
                    marg += p_o_stem[anc][g] * clade_tip_like[g];
                }
                ll += marg.max(1e-300).ln();
            } else {
                for &o_idx in &og_idxs[..n_og] {
                    ll += p_o[anc][o_idx].max(1e-15).ln();
                }
            }
        }
        if ll > best_ll {
            best_ll = ll;
            best_b = anc;
        }
    }
    bases[best_b]
}

/// Append uncovered child intervals (`>= min_span_bp`) as carry-forward union segments on non-root
/// internal nodes so child-exclusive regions deleted in the sister lineage (or missed by the
/// nearest outgroups) stay visible to deeper clades higher up the guide tree.
///
/// Spliced before `finalize_draft_ancestral_segments` so `left` uncovered spans interleave in
/// exact syntenic order between aligned blocks on the same chromosome, preserving collinearity.
#[allow(clippy::too_many_lines)]
fn append_child_union_drafts(
    drafts: &mut Vec<DraftAncSegment>,
    left_contigs: &[(String, Vec<u8>)],
    right_contigs: &[(String, Vec<u8>)],
    min_span_bp: usize,
) {
    let mut left_covered: HashMap<&str, IntervalSet> = HashMap::new();
    let mut right_covered: HashMap<&str, IntervalSet> = HashMap::new();
    let mut scratch = Vec::new();
    for d in drafts.iter() {
        if let Some((ref ctg, start, strand_ch, src_size)) = d.left_meta {
            let size = d.left_child_aln.iter().filter(|&&b| b != b'-').count();
            if size > 0 {
                let strand = Strand::from_char(strand_ch).unwrap_or(Strand::Plus);
                let (fs, fe) = Interval::ucsc_to_fwd_interval(start, size, strand, src_size);
                left_covered
                    .entry(ctg.as_str())
                    .or_default()
                    .insert(fs, fe, &mut scratch);
            }
        }
        if let Some((ref ctg, start, strand_ch, src_size)) = d.right_meta {
            let size = d.right_child_aln.iter().filter(|&&b| b != b'-').count();
            if size > 0 {
                let strand = Strand::from_char(strand_ch).unwrap_or(Strand::Plus);
                let (fs, fe) = Interval::ucsc_to_fwd_interval(start, size, strand, src_size);
                right_covered
                    .entry(ctg.as_str())
                    .or_default()
                    .insert(fs, fe, &mut scratch);
            }
        }
    }

    let mut new_drafts: Vec<DraftAncSegment> = Vec::new();
    let mut add_uncovered =
        |contigs: &[(String, Vec<u8>)], covered: &HashMap<&str, IntervalSet>, is_right: bool| {
            for (ctg_name, seq) in contigs {
                if seq.len() < min_span_bp {
                    continue;
                }
                let mut avail = IntervalSet::from_single(0, seq.len());
                if let Some(cov) = covered.get(ctg_name.as_str()) {
                    for (cs, ce) in cov.iter() {
                        avail.subtract(cs, ce, &mut scratch);
                    }
                }
                for (us, ue) in avail.iter() {
                    if ue.saturating_sub(us) < min_span_bp {
                        continue;
                    }
                    // Split uncovered intervals across runs of >= 8 non-ACGT ('N') bases so scaffold
                    // gaps are never copied into reconstructed ancestral contigs.
                    let mut pos = us;
                    while pos < ue {
                        while pos < ue && crate::bocpd::base_to_idx(seq[pos]).is_none() {
                            pos += 1;
                        }
                        let seg_start = pos;
                        let mut last_acgt = pos;
                        let mut n_run = 0usize;
                        while pos < ue {
                            if crate::bocpd::base_to_idx(seq[pos]).is_some() {
                                last_acgt = pos + 1;
                                n_run = 0;
                            } else {
                                n_run += 1;
                                if n_run >= 8 {
                                    break;
                                }
                            }
                            pos += 1;
                        }
                        let seg_end = last_acgt;
                        if seg_end.saturating_sub(seg_start) < min_span_bp {
                            continue;
                        }
                        let span_seq: Vec<u8> = seq[seg_start..seg_end]
                            .iter()
                            .map(u8::to_ascii_uppercase)
                            .collect();
                        let meta = Some((ctg_name.clone(), seg_start, '+', seq.len()));
                        if is_right {
                            new_drafts.push(DraftAncSegment {
                                anchor_is_right: true,
                                anchor_contig: ctg_name.clone(),
                                anchor_fwd_start: seg_start,
                                anc_ungapped: span_seq.clone(),
                                left_meta: None,
                                left_anc_aln: Vec::new(),
                                left_child_aln: Vec::new(),
                                right_meta: meta,
                                right_anc_aln: span_seq.clone(),
                                right_child_aln: span_seq,
                            });
                        } else {
                            new_drafts.push(DraftAncSegment {
                                anchor_is_right: false,
                                anchor_contig: ctg_name.clone(),
                                anchor_fwd_start: seg_start,
                                anc_ungapped: span_seq.clone(),
                                left_meta: meta,
                                left_anc_aln: span_seq.clone(),
                                left_child_aln: span_seq,
                                right_meta: None,
                                right_anc_aln: Vec::new(),
                                right_child_aln: Vec::new(),
                            });
                        }
                    }
                }
            }
        };

    add_uncovered(left_contigs, &left_covered, false);
    add_uncovered(right_contigs, &right_covered, true);
    drafts.extend(new_drafts);
}

/// Minimum uncovered child interval length (bp) carried forward into non-root internal ancestors.
const MIN_UNION_CARRY_BP: usize = 64;

type ChildContigPair<'a> = (&'a [(String, Vec<u8>)], &'a [(String, Vec<u8>)]);

/// Branch-length-aware probabilistic ancestral reconstruction combining
/// Felsenstein peeling over continuous-time HKY85 `P(τ) = exp(τ Q)` with a 2-state CTMC
/// indel-polarity model (bounded to `<= 8` outgroup tips under a 2-level stem/tip hierarchy).
#[must_use]
pub fn reconstruct_ancestral_contigs_transducer(
    anc_name: &str,
    left_name: &str,
    right_name: &str,
    outgroups: &[String],
    subproblem_blocks: &[TbaBlock],
    branches: AncestralBranchLengths,
) -> AncestralNodeResult {
    reconstruct_ancestral_contigs_transducer_inner(
        anc_name,
        left_name,
        right_name,
        outgroups,
        subproblem_blocks,
        branches,
        None,
    )
}

#[allow(clippy::too_many_lines)]
pub(super) fn reconstruct_ancestral_contigs_transducer_inner(
    anc_name: &str,
    left_name: &str,
    right_name: &str,
    outgroups: &[String],
    subproblem_blocks: &[TbaBlock],
    branches: AncestralBranchLengths,
    union_children: Option<ChildContigPair<'_>>,
) -> AncestralNodeResult {
    let og_set: HashSet<&str> = outgroups.iter().map(String::as_str).collect();
    let mut stats = AncestralNodeStats {
        name: anc_name.to_owned(),
        left_child: left_name.to_owned(),
        right_child: right_name.to_owned(),
        outgroups: outgroups.to_vec(),
        left_t: Some(branches.left_t),
        right_t: Some(branches.right_t),
        outgroup_t: Some(branches.outgroup_t),
        ..AncestralNodeStats::default()
    };

    // Estimate empirical base composition across subproblem blocks for HKY85 kernel.
    let mut base_cnt = [1.0_f64; 4];
    for blk in subproblem_blocks {
        for r in &blk.rows {
            for &b in &r.seq {
                if let Some(idx) = crate::bocpd::base_to_idx(b) {
                    base_cnt[idx] += 1.0;
                }
            }
        }
    }
    let kernel = crate::mdl::MarkovKernel::new(base_cnt, 2.0);
    let p_l = kernel.transition_matrix(branches.left_t);
    let p_r = kernel.transition_matrix(branches.right_t);
    let p_o = kernel.transition_matrix(branches.outgroup_t);
    let use_clade_peeling = outgroups.len() >= 2 && branches.outgroup_tip_t > 1e-6;
    let p_o_stem = kernel.transition_matrix(branches.outgroup_stem_t.max(1e-5));
    let p_o_tip = kernel.transition_matrix(branches.outgroup_tip_t.max(1e-5));

    let (surv_l, ins_l) = tkf91_probs(branches.left_t, branches.ins_rate, branches.del_rate);
    let (surv_r, ins_r) = tkf91_probs(branches.right_t, branches.ins_rate, branches.del_rate);
    let (surv_o, ins_o) = tkf91_probs(branches.outgroup_t, branches.ins_rate, branches.del_rate);
    let (surv_o_stem, ins_o_stem) = tkf91_probs(
        branches.outgroup_stem_t.max(1e-5),
        branches.ins_rate,
        branches.del_rate,
    );
    let (surv_o_tip, ins_o_tip) = tkf91_probs(
        branches.outgroup_tip_t.max(1e-5),
        branches.ins_rate,
        branches.del_rate,
    );
    // Use an unbiased 0.5/0.5 column presence prior so unpolarized 1-child insertions at the root
    // are not forced to ancestral presence by the geometric length parameter λ/μ.
    let prior_pres = 0.5_f64;
    let prior_abs = 0.5_f64;

    let mut drafts: Vec<DraftAncSegment> = Vec::new();
    let mut og_col_bases = Vec::new();

    for blk in subproblem_blocks {
        let left_row = blk.rows.iter().find(|r| r.species == left_name);
        let right_row = blk.rows.iter().find(|r| r.species == right_name);
        let og_rows: Vec<&TbaRow> = blk
            .rows
            .iter()
            .filter(|r| og_set.contains(r.species.as_str()))
            .collect();
        if left_row.is_none() && right_row.is_none() {
            continue;
        }
        let n_cols = blk.n_cols();
        let mut anc_ungapped = Vec::with_capacity(n_cols);
        let mut left_anc_aln = Vec::with_capacity(n_cols);
        let mut left_child_aln = Vec::with_capacity(n_cols);
        let mut right_anc_aln = Vec::with_capacity(n_cols);
        let mut right_child_aln = Vec::with_capacity(n_cols);

        for c in 0..n_cols {
            let lb = left_row
                .and_then(|r| r.seq.get(c).copied())
                .filter(|&b| b != b'-');
            let rb = right_row
                .and_then(|r| r.seq.get(c).copied())
                .filter(|&b| b != b'-');
            // Skip columns absent in both ingroup children (cannot anchor an ancestral base).
            if lb.is_none() && rb.is_none() {
                continue;
            }
            og_col_bases.clear();
            let mut og_pres_cnt = 0usize;
            for r in &og_rows {
                if let Some(&ob) = r.seq.get(c)
                    && ob != b'-'
                {
                    og_pres_cnt += 1;
                    og_col_bases.push(ob);
                }
            }
            // Count ALL outgroups missing from the block or gapped at column c as absent.
            let og_abs_cnt = outgroups.len().saturating_sub(og_pres_cnt);

            // 2-state CTMC indel-polarity posterior for ancestral site presence.
            let mut ll_pres = prior_pres.ln();
            let mut ll_abs = prior_abs.ln();
            ll_pres += if lb.is_some() {
                surv_l.ln()
            } else {
                (1.0 - surv_l).max(1e-12).ln()
            };
            ll_abs += if lb.is_some() {
                ins_l.ln()
            } else {
                (1.0 - ins_l).max(1e-12).ln()
            };
            ll_pres += if rb.is_some() {
                surv_r.ln()
            } else {
                (1.0 - surv_r).max(1e-12).ln()
            };
            ll_abs += if rb.is_some() {
                ins_r.ln()
            } else {
                (1.0 - ins_r).max(1e-12).ln()
            };
            if !outgroups.is_empty() {
                if use_clade_peeling {
                    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                    let np = og_pres_cnt as i32;
                    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                    let na = og_abs_cnt as i32;
                    let tip_given_1 = surv_o_tip.powi(np) * (1.0 - surv_o_tip).max(1e-12).powi(na);
                    let tip_given_0 = ins_o_tip.powi(np) * (1.0 - ins_o_tip).max(1e-12).powi(na);
                    let p_og_given_pres =
                        surv_o_stem * tip_given_1 + (1.0 - surv_o_stem).max(1e-12) * tip_given_0;
                    let p_og_given_abs =
                        ins_o_stem * tip_given_1 + (1.0 - ins_o_stem).max(1e-12) * tip_given_0;
                    ll_pres += p_og_given_pres.max(1e-300).ln();
                    ll_abs += p_og_given_abs.max(1e-300).ln();
                } else {
                    #[allow(clippy::cast_precision_loss)]
                    let pres_f = og_pres_cnt as f64;
                    #[allow(clippy::cast_precision_loss)]
                    let abs_f = og_abs_cnt as f64;
                    ll_pres += pres_f * surv_o.ln() + abs_f * (1.0 - surv_o).max(1e-12).ln();
                    ll_abs += pres_f * ins_o.ln() + abs_f * (1.0 - ins_o).max(1e-12).ln();
                }
            }

            // Require strict likelihood dominance (epsilon 1e-9) so unpolarized equal-branch
            // 1-child insertions at the root are not fabricated as ancestral bases.
            let anc_present = (lb.is_some() && rb.is_some()) || ll_pres > ll_abs + 1e-9;
            if anc_present {
                let anc_b = felsenstein_ml_base(
                    lb,
                    rb,
                    &og_col_bases,
                    kernel.pi,
                    &p_l,
                    &p_r,
                    &p_o,
                    &p_o_stem,
                    &p_o_tip,
                    use_clade_peeling,
                );
                anc_ungapped.push(anc_b);
                if lb.is_some() && rb.is_some() {
                    stats.shared_both_children_bp += 1;
                } else {
                    stats.outgroup_polarized_deletion_bp += 1;
                }
                if let Some(l) = lb {
                    left_anc_aln.push(anc_b);
                    left_child_aln.push(l.to_ascii_uppercase());
                } else if left_row.is_some() {
                    left_anc_aln.push(anc_b);
                    left_child_aln.push(b'-');
                }
                if let Some(r) = rb {
                    right_anc_aln.push(anc_b);
                    right_child_aln.push(r.to_ascii_uppercase());
                } else if right_row.is_some() {
                    right_anc_aln.push(anc_b);
                    right_child_aln.push(b'-');
                }
            } else {
                stats.lineage_insertion_excluded_bp +=
                    usize::from(lb.is_some()) + usize::from(rb.is_some());
                if let Some(l) = lb {
                    left_anc_aln.push(b'-');
                    left_child_aln.push(l.to_ascii_uppercase());
                }
                if let Some(r) = rb {
                    right_anc_aln.push(b'-');
                    right_child_aln.push(r.to_ascii_uppercase());
                }
            }
        }

        if anc_ungapped.is_empty() {
            continue;
        }
        let (anchor_is_right, anchor_contig, anchor_fwd_start) = if let Some(lr) = left_row {
            (false, lr.contig.clone(), lr.fwd_interval().0)
        } else if let Some(rr) = right_row {
            (true, rr.contig.clone(), rr.fwd_interval().0)
        } else {
            continue;
        };

        drafts.push(DraftAncSegment {
            anchor_is_right,
            anchor_contig,
            anchor_fwd_start,
            anc_ungapped,
            left_meta: left_row.map(|r| (r.contig.clone(), r.start, r.strand, r.src_size)),
            left_anc_aln,
            left_child_aln,
            right_meta: right_row.map(|r| (r.contig.clone(), r.start, r.strand, r.src_size)),
            right_anc_aln,
            right_child_aln,
        });
    }

    if let Some((left_ctgs, right_ctgs)) = union_children {
        append_child_union_drafts(&mut drafts, left_ctgs, right_ctgs, MIN_UNION_CARRY_BP);
    }

    finalize_draft_ancestral_segments(anc_name, left_name, right_name, drafts, stats)
}

/// Reconstruct an outgroup-polarized ancestral genome (`anc_name`) and exact parent->child
/// pairwise alignments from a 2–4 genome subproblem TBA blockset.
///
/// Inclusion rule (Progressive Cactus indel polarization):
/// - Columns present in **both children** are retained via Fitch parsimony (ties broken by outgroups).
/// - Columns present in **one child + >=1 outgroup** are retained (lineage-specific deletion in the
///   other child).
/// - Columns present in **only one child and 0 outgroups** are excluded from the ancestral genome
///   (lineage-specific insertion).
#[must_use]
pub fn reconstruct_ancestral_contigs(
    anc_name: &str,
    left_name: &str,
    right_name: &str,
    outgroups: &[String],
    subproblem_blocks: &[TbaBlock],
) -> AncestralNodeResult {
    reconstruct_ancestral_contigs_inner(
        anc_name,
        left_name,
        right_name,
        outgroups,
        subproblem_blocks,
        None,
    )
}

#[allow(clippy::too_many_lines)]
pub(super) fn reconstruct_ancestral_contigs_inner(
    anc_name: &str,
    left_name: &str,
    right_name: &str,
    outgroups: &[String],
    subproblem_blocks: &[TbaBlock],
    union_children: Option<ChildContigPair<'_>>,
) -> AncestralNodeResult {
    let og_set: HashSet<&str> = outgroups.iter().map(String::as_str).collect();
    let mut stats = AncestralNodeStats {
        name: anc_name.to_owned(),
        left_child: left_name.to_owned(),
        right_child: right_name.to_owned(),
        outgroups: outgroups.to_vec(),
        ..AncestralNodeStats::default()
    };

    let mut drafts: Vec<DraftAncSegment> = Vec::new();
    for blk in subproblem_blocks {
        let left_row = blk.rows.iter().find(|r| r.species == left_name);
        let right_row = blk.rows.iter().find(|r| r.species == right_name);
        let og_rows: Vec<&TbaRow> = blk
            .rows
            .iter()
            .filter(|r| og_set.contains(r.species.as_str()))
            .collect();

        if left_row.is_none() && right_row.is_none() {
            continue;
        }
        if (left_row.is_none() || right_row.is_none()) && og_rows.is_empty() {
            if let Some(r) = left_row.or(right_row) {
                stats.lineage_insertion_excluded_bp += r.size;
            }
            continue;
        }

        let n_cols = blk.n_cols();
        let mut anc_ungapped = Vec::with_capacity(n_cols);
        let mut left_anc_aln = Vec::with_capacity(n_cols);
        let mut left_child_aln = Vec::with_capacity(n_cols);
        let mut right_anc_aln = Vec::with_capacity(n_cols);
        let mut right_child_aln = Vec::with_capacity(n_cols);

        for c in 0..n_cols {
            let lb = left_row
                .and_then(|r| r.seq.get(c).copied())
                .filter(|&b| b != b'-');
            let rb = right_row
                .and_then(|r| r.seq.get(c).copied())
                .filter(|&b| b != b'-');
            let mut og_mask = 0u8;
            let mut og_has_base = false;
            for r in &og_rows {
                if let Some(&ob) = r.seq.get(c)
                    && ob != b'-'
                {
                    og_has_base = true;
                    og_mask |= fitch_mask(ob);
                }
            }

            match (lb, rb) {
                (Some(l), Some(r)) => {
                    let lm = fitch_mask(l);
                    let rm = fitch_mask(r);
                    let mask = if lm & rm != 0 {
                        lm & rm
                    } else if (lm | rm) & og_mask != 0 {
                        (lm | rm) & og_mask
                    } else {
                        lm | rm
                    };
                    let anc_b = pick_base_from_mask(mask, l);
                    anc_ungapped.push(anc_b);
                    stats.shared_both_children_bp += 1;
                    left_anc_aln.push(anc_b);
                    left_child_aln.push(l.to_ascii_uppercase());
                    right_anc_aln.push(anc_b);
                    right_child_aln.push(r.to_ascii_uppercase());
                }
                (Some(l), None) => {
                    if og_has_base {
                        let lm = fitch_mask(l);
                        let mask = if lm & og_mask != 0 { lm & og_mask } else { lm };
                        let anc_b = pick_base_from_mask(mask, l);
                        anc_ungapped.push(anc_b);
                        stats.outgroup_polarized_deletion_bp += 1;
                        left_anc_aln.push(anc_b);
                        left_child_aln.push(l.to_ascii_uppercase());
                        if right_row.is_some() {
                            right_anc_aln.push(anc_b);
                            right_child_aln.push(b'-');
                        }
                    } else {
                        stats.lineage_insertion_excluded_bp += 1;
                        left_anc_aln.push(b'-');
                        left_child_aln.push(l.to_ascii_uppercase());
                    }
                }
                (None, Some(r)) => {
                    if og_has_base {
                        let rm = fitch_mask(r);
                        let mask = if rm & og_mask != 0 { rm & og_mask } else { rm };
                        let anc_b = pick_base_from_mask(mask, r);
                        anc_ungapped.push(anc_b);
                        stats.outgroup_polarized_deletion_bp += 1;
                        right_anc_aln.push(anc_b);
                        right_child_aln.push(r.to_ascii_uppercase());
                        if left_row.is_some() {
                            left_anc_aln.push(anc_b);
                            left_child_aln.push(b'-');
                        }
                    } else {
                        stats.lineage_insertion_excluded_bp += 1;
                        right_anc_aln.push(b'-');
                        right_child_aln.push(r.to_ascii_uppercase());
                    }
                }
                (None, None) => {}
            }
        }

        if anc_ungapped.is_empty() {
            continue;
        }
        let (anchor_is_right, anchor_contig, anchor_fwd_start) = if let Some(lr) = left_row {
            (false, lr.contig.clone(), lr.fwd_interval().0)
        } else if let Some(rr) = right_row {
            (true, rr.contig.clone(), rr.fwd_interval().0)
        } else {
            continue;
        };

        drafts.push(DraftAncSegment {
            anchor_is_right,
            anchor_contig,
            anchor_fwd_start,
            anc_ungapped,
            left_meta: left_row.map(|r| (r.contig.clone(), r.start, r.strand, r.src_size)),
            left_anc_aln,
            left_child_aln,
            right_meta: right_row.map(|r| (r.contig.clone(), r.start, r.strand, r.src_size)),
            right_anc_aln,
            right_child_aln,
        });
    }

    if let Some((left_ctgs, right_ctgs)) = union_children {
        append_child_union_drafts(&mut drafts, left_ctgs, right_ctgs, MIN_UNION_CARRY_BP);
    }

    finalize_draft_ancestral_segments(anc_name, left_name, right_name, drafts, stats)
}

fn finalize_draft_ancestral_segments(
    anc_name: &str,
    left_name: &str,
    right_name: &str,
    mut drafts: Vec<DraftAncSegment>,
    mut stats: AncestralNodeStats,
) -> AncestralNodeResult {
    // Sort drafts deterministically by synteny anchor so colinear blocks chain into contigs.
    drafts.sort_by(|a, b| {
        a.anchor_is_right
            .cmp(&b.anchor_is_right)
            .then_with(|| a.anchor_contig.cmp(&b.anchor_contig))
            .then_with(|| a.anchor_fwd_start.cmp(&b.anchor_fwd_start))
    });

    // Chain segments sharing the same (anchor_is_right, anchor_contig) into ancestral contigs.
    let mut anc_contigs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut placements: Vec<(usize, usize, DraftAncSegment)> = Vec::with_capacity(drafts.len());
    let mut cur_key: Option<(bool, String)> = None;
    for d in drafts {
        let key = (d.anchor_is_right, d.anchor_contig.clone());
        let ctg_idx = if cur_key.as_ref() == Some(&key) {
            anc_contigs.len() - 1
        } else {
            cur_key = Some(key);
            let idx = anc_contigs.len();
            anc_contigs.push((format!("{anc_name}_ctg{idx}"), Vec::new()));
            idx
        };
        let anc_start = anc_contigs[ctg_idx].1.len();
        anc_contigs[ctg_idx].1.extend_from_slice(&d.anc_ungapped);
        placements.push((ctg_idx, anc_start, d));
    }

    let mut anc_to_left = Vec::new();
    let mut anc_to_right = Vec::new();
    for (ctg_idx, anc_start, d) in placements {
        let (anc_ctg_name, anc_ctg_seq) = &anc_contigs[ctg_idx];
        let anc_src_size = anc_ctg_seq.len();
        let anc_src = format!("{anc_name}.{anc_ctg_name}");

        if let Some((l_ctg, l_start, l_strand, l_src_size)) = d.left_meta
            && d.left_child_aln.iter().any(|&b| b != b'-')
            && d.left_anc_aln.iter().any(|&b| b != b'-')
        {
            let ra = MafRow {
                src: anc_src.clone(),
                start: anc_start,
                strand: '+',
                src_size: anc_src_size,
                seq: String::from_utf8_lossy(&d.left_anc_aln).into_owned(),
            };
            let rb = MafRow {
                src: format!("{left_name}.{l_ctg}"),
                start: l_start,
                strand: l_strand,
                src_size: l_src_size,
                seq: String::from_utf8_lossy(&d.left_child_aln).into_owned(),
            };
            if let Some(pb) = PairwiseTbaBlock::from_maf_rows(anc_name, ra, left_name, rb) {
                anc_to_left.push(pb);
            }
        }

        if let Some((r_ctg, r_start, r_strand, r_src_size)) = d.right_meta
            && d.right_child_aln.iter().any(|&b| b != b'-')
            && d.right_anc_aln.iter().any(|&b| b != b'-')
        {
            let ra = MafRow {
                src: anc_src,
                start: anc_start,
                strand: '+',
                src_size: anc_src_size,
                seq: String::from_utf8_lossy(&d.right_anc_aln).into_owned(),
            };
            let rb = MafRow {
                src: format!("{right_name}.{r_ctg}"),
                start: r_start,
                strand: r_strand,
                src_size: r_src_size,
                seq: String::from_utf8_lossy(&d.right_child_aln).into_owned(),
            };
            if let Some(pb) = PairwiseTbaBlock::from_maf_rows(anc_name, ra, right_name, rb) {
                anc_to_right.push(pb);
            }
        }
    }

    stats.n_anc_contigs = anc_contigs.len();
    stats.anc_bp = anc_contigs.iter().map(|(_, s)| s.len()).sum();

    AncestralNodeResult {
        anc_contigs,
        anc_to_left,
        anc_to_right,
        stats,
    }
}
