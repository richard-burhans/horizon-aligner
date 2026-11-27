//! Linear-scaling partial-order alignment (POA) refinement of shared conspecific tip insertions
//! inside hinge-gap windows.

use super::model::TbaBlock;
use rayon::prelude::*;
use std::collections::HashSet;

const MIN_POA_INS_BP: usize = 8;
const MAX_POA_INS_BP: usize = 2048;
const POA_MIN_BLOCK_BITS: f64 = 12.0;

fn consensus_base_for_col(matrix: &[Vec<u8>], col: usize) -> u8 {
    let mut counts = [0u32; 4];
    for row in matrix {
        match row[col].to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    let bases = *b"ACGT";
    let mut best_i = 0usize;
    let mut best_c = 0u32;
    for (i, &c) in counts.iter().enumerate() {
        if c > best_c {
            best_c = c;
            best_i = i;
        }
    }
    if best_c == 0 { b'N' } else { bases[best_i] }
}

/// Build a compact topological partial-order DAG from the current progressive MSA column matrix.
/// Consecutive columns with identical non-gap row membership collapse into a single `GraphNode`
/// carrying majority-vote consensus bases; edges follow every aligned member's path so indels
/// form clean skip branches. Returns `(nodes, node_col_starts)`.
fn build_poa_dag_from_matrix(matrix: &[Vec<u8>]) -> (Vec<crate::mdl::GraphNode>, Vec<usize>) {
    let width = matrix.first().map_or(0, Vec::len);
    if width == 0 {
        return (Vec::new(), Vec::new());
    }
    let same_mask = |c1: usize, c2: usize| -> bool {
        matrix.iter().all(|r| (r[c1] != b'-') == (r[c2] != b'-'))
    };
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut start = 0usize;
    for c in 1..width {
        if !same_mask(start, c) {
            runs.push((start, c));
            start = c;
        }
    }
    runs.push((start, width));

    let num_nodes = runs.len();
    let mut preds: Vec<Vec<usize>> = vec![Vec::new(); num_nodes];
    for row in matrix {
        let mut prev_node: Option<usize> = None;
        for (ni, &(lo, _)) in runs.iter().enumerate() {
            if row[lo] != b'-' {
                if let Some(pn) = prev_node
                    && !preds[ni].contains(&pn)
                {
                    preds[ni].push(pn);
                }
                prev_node = Some(ni);
            }
        }
    }
    for p in &mut preds {
        p.sort_unstable();
    }

    let mut nodes = Vec::with_capacity(num_nodes);
    let mut col_starts = Vec::with_capacity(num_nodes);
    for (ni, &(lo, hi)) in runs.iter().enumerate() {
        let seq: Vec<u8> = (lo..hi)
            .map(|c| consensus_base_for_col(matrix, c))
            .collect();
        nodes.push(crate::mdl::GraphNode {
            name: format!("poa_{ni}"),
            seq,
            preds: std::mem::take(&mut preds[ni]),
        });
        col_starts.push(lo);
    }
    (nodes, col_starts)
}

/// Thread `seq` into `matrix` along strictly increasing `(msa_col, query_pos)` anchor pairs,
/// preserving every existing row's and `seq`'s ungapped bases byte-for-byte.
fn thread_seq_into_poa_matrix(matrix: &mut Vec<Vec<u8>>, seq: &[u8], pairs: &[(usize, usize)]) {
    let old_width = matrix.first().map_or(0, Vec::len);
    let n_existing = matrix.len();
    let mut next_matrix: Vec<Vec<u8>> = (0..=n_existing)
        .map(|_| Vec::with_capacity(old_width + seq.len()))
        .collect();

    let mut cur_col = 0usize;
    let mut cur_q = 0usize;
    for &(mcol, qpos) in pairs {
        if mcol < cur_col || qpos < cur_q || mcol >= old_width || qpos >= seq.len() {
            continue;
        }
        while cur_col < mcol {
            for r in 0..n_existing {
                next_matrix[r].push(matrix[r][cur_col]);
            }
            next_matrix[n_existing].push(b'-');
            cur_col += 1;
        }
        while cur_q < qpos {
            for row in next_matrix.iter_mut().take(n_existing) {
                row.push(b'-');
            }
            next_matrix[n_existing].push(seq[cur_q]);
            cur_q += 1;
        }
        for r in 0..n_existing {
            next_matrix[r].push(matrix[r][mcol]);
        }
        next_matrix[n_existing].push(seq[qpos]);
        cur_col = mcol + 1;
        cur_q = qpos + 1;
    }

    while cur_col < old_width {
        for r in 0..n_existing {
            next_matrix[r].push(matrix[r][cur_col]);
        }
        next_matrix[n_existing].push(b'-');
        cur_col += 1;
    }
    while cur_q < seq.len() {
        for row in next_matrix.iter_mut().take(n_existing) {
            row.push(b'-');
        }
        next_matrix[n_existing].push(seq[cur_q]);
        cur_q += 1;
    }

    *matrix = next_matrix;
}

/// Align participating non-hinge tip insertions inside a single hinge-gap window `[c0..c1)`
/// using linear-scaling `O(N_tips)` DAG/POA refinement. Returns `true` if the window was
/// compressed and appended into `new_seqs`.
fn try_refine_hinge_gap_window(
    blk: &TbaBlock,
    c0: usize,
    c1: usize,
    hinge_species: &HashSet<&str>,
    params: &crate::mdl::HmmParams,
    new_seqs: &mut [Vec<u8>],
) -> bool {
    let mut eligible: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut bypass: Vec<(usize, Vec<u8>)> = Vec::new();
    for (ri, row) in blk.rows.iter().enumerate() {
        if hinge_species.contains(row.species.as_str()) {
            continue;
        }
        let ug: Vec<u8> = row.seq[c0..c1]
            .iter()
            .copied()
            .filter(|&b| b != b'-')
            .collect();
        if ug.is_empty() {
            continue;
        }
        if ug.len() >= MIN_POA_INS_BP && ug.len() <= MAX_POA_INS_BP {
            eligible.push((ri, ug));
        } else {
            bypass.push((ri, ug));
        }
    }
    if eligible.len() < 2 {
        return false;
    }

    // Deterministic seed order: longest insertion first, tie-broken by species name.
    eligible.sort_by(|a, b| {
        b.1.len()
            .cmp(&a.1.len())
            .then_with(|| blk.rows[a.0].species.cmp(&blk.rows[b.0].species))
    });

    let mut aligned_row_ids: Vec<usize> = vec![eligible[0].0];
    let mut matrix: Vec<Vec<u8>> = vec![eligible[0].1.clone()];
    let mut unaligned: Vec<(usize, Vec<u8>)> = Vec::new();

    for (ri, seq) in eligible.into_iter().skip(1) {
        let (dag_nodes, col_starts) = build_poa_dag_from_matrix(&matrix);
        let aln = crate::mdl::align_seq_to_dag_mdl(&dag_nodes, &seq, params, POA_MIN_BLOCK_BITS);
        let cur_width = matrix.first().map_or(0, Vec::len);
        let min_matches = (seq.len().min(cur_width) / 4).max(MIN_POA_INS_BP / 2);
        if aln.score_bits >= POA_MIN_BLOCK_BITS && aln.kept_pairs.len() >= min_matches {
            let mapped_pairs: Vec<(usize, usize)> = aln
                .kept_pairs
                .into_iter()
                .map(|(nid, npos, qpos)| (col_starts[nid] + npos, qpos))
                .collect();
            thread_seq_into_poa_matrix(&mut matrix, &seq, &mapped_pairs);
            aligned_row_ids.push(ri);
        } else {
            unaligned.push((ri, seq));
        }
    }

    let poa_width = matrix.first().map_or(0, Vec::len);
    let extra_width: usize = unaligned
        .iter()
        .chain(bypass.iter())
        .map(|(_, s)| s.len())
        .sum();
    if aligned_row_ids.len() < 2 || poa_width + extra_width >= c1 - c0 {
        return false;
    }

    for (ri, row_buf) in new_seqs.iter_mut().enumerate() {
        if let Some(pos) = aligned_row_ids.iter().position(|&id| id == ri) {
            row_buf.extend_from_slice(&matrix[pos]);
        } else {
            row_buf.resize(row_buf.len() + poa_width, b'-');
        }
    }
    for (u_ri, u_seq) in unaligned.into_iter().chain(bypass) {
        let ulen = u_seq.len();
        for (ri, row_buf) in new_seqs.iter_mut().enumerate() {
            if ri == u_ri {
                row_buf.extend_from_slice(&u_seq);
            } else {
                row_buf.resize(row_buf.len() + ulen, b'-');
            }
        }
    }
    true
}

/// Refine shared conspecific tip insertions inside hinge-gap windows using linear-scaling
/// `O(N_tips)` partial-order / DAG MDL alignment (`align_seq_to_dag_mdl`).
///
/// For any contiguous column range `[c0..c1)` where every hinge species present in a block has
/// gaps (`-`) and at least two non-hinge tip rows carry non-gap inserted sequence, aligns the
/// participating tip insertions progressively against a compact DAG so homologous insertions
/// collapse into unified multi-tip columns while non-homologous insertions remain unaligned under
/// the MDL gate. Every row's ungapped DNA sequence is preserved byte-for-byte.
pub fn refine_shared_tip_insertions_poa<S: std::hash::BuildHasher>(
    blocks: &mut [TbaBlock],
    hinge_species: &HashSet<String, S>,
) -> usize {
    let params = crate::mdl::HmmParams::default();
    let hinge_refs: HashSet<&str> = hinge_species.iter().map(String::as_str).collect();

    blocks
        .par_iter_mut()
        .map(|blk| {
            let n_rows = blk.rows.len();
            let n_cols = blk.n_cols();
            if n_rows < 3 || n_cols == 0 {
                return 0usize;
            }
            let hinge_rows: Vec<usize> = (0..n_rows)
                .filter(|&ri| hinge_refs.contains(blk.rows[ri].species.as_str()))
                .collect();
            if hinge_rows.is_empty() || n_rows < hinge_rows.len() + 2 {
                return 0usize;
            }

            let mut new_seqs: Vec<Vec<u8>> =
                (0..n_rows).map(|_| Vec::with_capacity(n_cols)).collect();
            let mut windows_refined = 0usize;
            let mut c = 0usize;

            while c < n_cols {
                if !hinge_rows.iter().all(|&ri| blk.rows[ri].seq[c] == b'-') {
                    for (row_buf, row) in new_seqs.iter_mut().zip(&blk.rows) {
                        row_buf.push(row.seq[c]);
                    }
                    c += 1;
                    continue;
                }
                let c0 = c;
                while c < n_cols && hinge_rows.iter().all(|&ri| blk.rows[ri].seq[c] == b'-') {
                    c += 1;
                }
                let c1 = c;

                if try_refine_hinge_gap_window(blk, c0, c1, &hinge_refs, &params, &mut new_seqs) {
                    windows_refined += 1;
                } else {
                    for (row_buf, row) in new_seqs.iter_mut().zip(&blk.rows) {
                        row_buf.extend_from_slice(&row.seq[c0..c1]);
                    }
                }
            }

            if windows_refined > 0 {
                for (row, seq) in blk.rows.iter_mut().zip(new_seqs) {
                    row.seq = seq;
                }
            }
            windows_refined
        })
        .sum()
}
