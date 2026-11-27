//! Probabilistic triplet consistency transformation (`Pecan` / `T-Coffee` / `ProbCons`)
//! over sparse Pair-HMM posterior homology matrices.

use crate::bocpd::{
    NEG_INF, base_probs, bocpd_gc_profile, joint_pair_bits, logsumexp2, logsumexp3,
};
use crate::mdl::HmmParams;

/// Maximum Forward matrix cells (`f64`, 8 bytes/cell = 128 MiB) allowed for a single Pair-HMM
/// posterior computation before tightening the diagonal band.
const MAX_PAIR_HMM_CELLS: usize = 16_000_000;
/// Full-DP threshold (`~12 MiB` for `fm`) above which Forward-Backward switches automatically to
/// a diagonal band around the main diagonal.
const BANDED_THRESHOLD_CELLS: usize = 1_500_000;
/// Default half-band width (columns on each side of the shifted main diagonal) for long loci.
const DEFAULT_HALF_BAND: usize = 256;

/// Sparse row-indexed pairwise posterior homology matrix `P(x_i ~ y_j) >= threshold`.
#[derive(Debug, Clone, Default)]
pub struct SparsePosteriorMatrix {
    pub len_x: usize,
    pub len_y: usize,
    /// `rows[i]` holds `(j, prob)` entries sorted by `j`.
    pub rows: Vec<Vec<(usize, f32)>>,
}

impl SparsePosteriorMatrix {
    /// Compute a sparse posterior homology matrix between `seq_x` and `seq_y` via rolling-row
    /// Pair-HMM Forward-Backward decoding (`P(x_i ~ y_j) >= min_prob`).
    ///
    /// Uses exact full Forward-Backward for small sequences (`(n + 1) * (m + 1) <= 1.5M` cells)
    /// and automatically engages a bounded diagonal band (`<= 128 MiB` peak for `fm`) on longer
    /// loci so 25 kb+ intervals never OOM.
    #[must_use]
    #[allow(clippy::too_many_lines, clippy::similar_names)]
    pub fn from_pair(seq_x: &[u8], seq_y: &[u8], params: &HmmParams, min_prob: f32) -> Self {
        let n = seq_x.len();
        let m = seq_y.len();
        if n == 0 || m == 0 {
            return Self {
                len_x: n,
                len_y: m,
                rows: vec![Vec::new(); n],
            };
        }

        let norm = params.normalized();
        let ln2 = 2.0_f64.ln();
        let go = norm.go_bits * ln2;
        let ge = norm.ge_bits * ln2;
        let mm = norm.mm_bits * ln2;
        let im = norm.im_bits * ln2;

        let gc_x = bocpd_gc_profile(seq_x, 1.0 / 120.0, 80);
        let gc_y = bocpd_gc_profile(seq_y, 1.0 / 120.0, 80);
        let pi_x: Vec<[f64; 4]> = gc_x.into_iter().map(base_probs).collect();
        let pi_y: Vec<[f64; 4]> = gc_y.into_iter().map(base_probs).collect();

        // Determine per-row band width (full DP when small; banded within MAX_PAIR_HMM_CELLS when large).
        let full_cells = (n + 1).saturating_mul(m + 1);
        let half_band = if full_cells <= BANDED_THRESHOLD_CELLS {
            m
        } else {
            let max_width = (MAX_PAIR_HMM_CELLS / (n + 1)).max(1);
            let max_half = max_width.saturating_sub(1) / 2;
            let len_diff = n.abs_diff(m);
            let desired = len_diff.saturating_add(DEFAULT_HALF_BAND).min(m);
            if desired <= max_half {
                desired
            } else if max_half >= len_diff {
                max_half
            } else {
                // Sequence lengths differ by more than our memory budget allows to bridge.
                return Self {
                    len_x: n,
                    len_y: m,
                    rows: vec![Vec::new(); n],
                };
            }
        };

        let band_bounds = |i: usize| -> (usize, usize) {
            if half_band >= m {
                return (0, m);
            }
            // Center diagonal proportionally across [0..=n] -> [0..=m].
            let center = i.saturating_mul(m).checked_div(n).unwrap_or(0);
            let lo = center.saturating_sub(half_band);
            let hi = center.saturating_add(half_band).min(m);
            (lo, hi)
        };

        let row_width = (2 * half_band + 1).min(m + 1);
        let Some(fm_cells) = (n + 1).checked_mul(row_width) else {
            return Self {
                len_x: n,
                len_y: m,
                rows: vec![Vec::new(); n],
            };
        };
        if fm_cells > MAX_PAIR_HMM_CELLS {
            return Self {
                len_x: n,
                len_y: m,
                rows: vec![Vec::new(); n],
            };
        }

        // Store ONLY `fm` across all rows (8 bytes/cell) using compact band offsets.
        let mut fm = vec![NEG_INF; fm_cells];
        let mut row_lo = vec![0usize; n + 1];
        let mut row_hi = vec![0usize; n + 1];
        for i in 0..=n {
            let (lo, hi) = band_bounds(i);
            row_lo[i] = lo;
            row_hi[i] = hi;
        }
        let fm_idx = |i: usize, j: usize, lo: usize| -> usize { i * row_width + (j - lo) };

        // Rolling 1D Forward buffers for Ix and Iy.
        let mut fix_prev = vec![NEG_INF; m + 1];
        let mut fix_cur = vec![NEG_INF; m + 1];
        let mut fiy_prev = vec![NEG_INF; m + 1];
        let mut fiy_cur = vec![NEG_INF; m + 1];

        // Row 0 initialization.
        fm[fm_idx(0, 0, row_lo[0])] = 0.0;
        for j in 1..=row_hi[0] {
            let prev_m = if j > row_lo[0] {
                fm[fm_idx(0, j - 1, row_lo[0])]
            } else {
                NEG_INF
            };
            fiy_prev[j] = logsumexp2(prev_m + go, fiy_prev[j - 1] + ge);
        }

        let mut last_fix = NEG_INF;
        let mut last_fiy = fiy_prev[m];

        for i in 1..=n {
            let lo = row_lo[i];
            let hi = row_hi[i];
            let prev_lo = row_lo[i - 1];
            let prev_hi = row_hi[i - 1];
            let bx = seq_x[i - 1];
            let px = pi_x[i - 1];

            for j in lo..=hi {
                fix_cur[j] = NEG_INF;
                fiy_cur[j] = NEG_INF;
            }
            if lo == 0 {
                let pm0 = if prev_lo == 0 {
                    fm[fm_idx(i - 1, 0, prev_lo)]
                } else {
                    NEG_INF
                };
                fix_cur[0] = logsumexp2(pm0 + go, fix_prev[0] + ge);
            }
            for j in lo.max(1)..=hi {
                let sb = joint_pair_bits(bx, seq_y[j - 1], px, pi_y[j - 1], norm.id_prior) * ln2;
                let (pm_diag, pix_diag, piy_diag) = if j > prev_lo && j - 1 <= prev_hi {
                    (
                        fm[fm_idx(i - 1, j - 1, prev_lo)],
                        fix_prev[j - 1],
                        fiy_prev[j - 1],
                    )
                } else {
                    (NEG_INF, NEG_INF, NEG_INF)
                };
                let cur_m = logsumexp3(pm_diag + mm, pix_diag + im, piy_diag + im) + sb;
                fm[fm_idx(i, j, lo)] = cur_m;

                let (pm_up, pix_up) = if j >= prev_lo && j <= prev_hi {
                    (fm[fm_idx(i - 1, j, prev_lo)], fix_prev[j])
                } else {
                    (NEG_INF, NEG_INF)
                };
                fix_cur[j] = logsumexp2(pm_up + go, pix_up + ge);

                let (pm_left, piy_left) = if j > lo {
                    (fm[fm_idx(i, j - 1, lo)], fiy_cur[j - 1])
                } else {
                    (NEG_INF, NEG_INF)
                };
                fiy_cur[j] = logsumexp2(pm_left + go, piy_left + ge);
            }
            if i == n {
                last_fix = fix_cur[m];
                last_fiy = fiy_cur[m];
            }
            std::mem::swap(&mut fix_prev, &mut fix_cur);
            std::mem::swap(&mut fiy_prev, &mut fiy_cur);
        }

        let last_fm = if m >= row_lo[n] && m <= row_hi[n] {
            fm[fm_idx(n, m, row_lo[n])]
        } else {
            NEG_INF
        };
        let z_tot = logsumexp3(last_fm, last_fix, last_fiy);
        if !z_tot.is_finite() || z_tot <= NEG_INF * 0.5 {
            return Self {
                len_x: n,
                len_y: m,
                rows: vec![Vec::new(); n],
            };
        }

        // Drop Forward rolling vectors before Backward sweep.
        drop(fix_prev);
        drop(fix_cur);
        drop(fiy_prev);
        drop(fiy_cur);

        // Rolling 1D Backward sweep computing sparse posteriors on the fly.
        let mut bm_next = vec![NEG_INF; m + 1];
        let mut bix_next = vec![NEG_INF; m + 1];
        let mut bm_cur = vec![NEG_INF; m + 1];
        let mut bix_cur = vec![NEG_INF; m + 1];
        let mut biy_cur = vec![NEG_INF; m + 1];
        let mut rows = vec![Vec::new(); n];

        for i in (0..=n).rev() {
            let lo = row_lo[i];
            let hi = row_hi[i];
            let (next_lo, next_hi) = if i < n {
                (row_lo[i + 1], row_hi[i + 1])
            } else {
                (m + 1, 0)
            };
            for j in lo..=hi {
                bm_cur[j] = NEG_INF;
                bix_cur[j] = NEG_INF;
                biy_cur[j] = NEG_INF;
            }
            for j in (lo..=hi).rev() {
                if i == n && j == m {
                    bm_cur[m] = 0.0;
                    bix_cur[m] = 0.0;
                    biy_cur[m] = 0.0;
                } else {
                    let diag = if i < n && j < m && j + 1 >= next_lo && j < next_hi {
                        let sb =
                            joint_pair_bits(seq_x[i], seq_y[j], pi_x[i], pi_y[j], norm.id_prior)
                                * ln2;
                        bm_next[j + 1] + sb
                    } else {
                        NEG_INF
                    };
                    let down = if i < n && j >= next_lo && j <= next_hi {
                        bix_next[j]
                    } else {
                        NEG_INF
                    };
                    let right = if j < m && j < hi {
                        biy_cur[j + 1]
                    } else {
                        NEG_INF
                    };
                    bm_cur[j] = logsumexp3(diag + mm, down + go, right + go);
                    bix_cur[j] = logsumexp2(diag + im, down + ge);
                    biy_cur[j] = logsumexp2(diag + im, right + ge);
                }

                if i >= 1 && j >= 1 {
                    let fm_val = fm[fm_idx(i, j, lo)];
                    let log_p = fm_val + bm_cur[j] - z_tot;
                    if log_p.is_finite() {
                        #[allow(clippy::cast_possible_truncation)]
                        let p = log_p.min(0.0).exp() as f32;
                        if p >= min_prob {
                            rows[i - 1].push((j - 1, p));
                        }
                    }
                }
            }
            if i >= 1 {
                rows[i - 1].reverse();
            }
            std::mem::swap(&mut bm_next, &mut bm_cur);
            std::mem::swap(&mut bix_next, &mut bix_cur);
        }

        Self {
            len_x: n,
            len_y: m,
            rows,
        }
    }

    /// Return the transpose `P^T(y_j ~ x_i)`.
    #[must_use]
    pub fn transpose(&self) -> Self {
        let mut rows = vec![Vec::new(); self.len_y];
        for (i, row) in self.rows.iter().enumerate() {
            for &(j, p) in row {
                if j < self.len_y {
                    rows[j].push((i, p));
                }
            }
        }
        Self {
            len_x: self.len_y,
            len_y: self.len_x,
            rows,
        }
    }

    /// Look up posterior `P(x_i ~ y_j)` (0.0 if absent).
    #[must_use]
    pub fn get(&self, i: usize, j: usize) -> f32 {
        let Some(row) = self.rows.get(i) else {
            return 0.0;
        };
        row.binary_search_by_key(&j, |&(col, _)| col)
            .map_or(0.0, |idx| row[idx].1)
    }
}

/// Apply `iterations` rounds of `ProbCons` / `Pecan` probabilistic triplet consistency
/// relaxation across all `N*(N-1)/2` pairwise posterior matrices:
/// `P'(x_i ~ y_j) = (2 * P(x_i ~ y_j) + Σ_{z ≠ x,y} Σ_k P(x_i ~ z_k) P(z_k ~ y_j)) / N`.
///
/// Missing transpose entries (`cur[y][x]` when `cur[x][y]` is present) are auto-populated so
/// callers passing an upper-triangular matrix receive unbiased consistency relaxation.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn triplet_consistency_transform(
    seq_lens: &[usize],
    pairwise: &[Vec<Option<SparsePosteriorMatrix>>],
    iterations: usize,
) -> Vec<Vec<Option<SparsePosteriorMatrix>>> {
    let n = seq_lens.len();
    if n < 3 || iterations == 0 || pairwise.len() < n || pairwise.iter().any(|r| r.len() < n) {
        return pairwise.to_vec();
    }
    let mut cur = pairwise.to_vec();
    // Auto-populate missing symmetric transposes so upper-triangular inputs work accurately.
    for x in 0..n {
        for y in (x + 1)..n {
            match (&cur[x][y], &cur[y][x]) {
                (Some(m), None) => cur[y][x] = Some(m.transpose()),
                (None, Some(m)) => cur[x][y] = Some(m.transpose()),
                _ => {}
            }
        }
    }

    #[allow(clippy::cast_precision_loss)]
    let inv_n = 1.0_f32 / (n as f32);
    let max_len = seq_lens.iter().copied().max().unwrap_or(0);
    let mut acc = vec![0.0_f32; max_len];
    let mut touched: Vec<usize> = Vec::with_capacity(64);

    for _ in 0..iterations {
        let mut next: Vec<Vec<Option<SparsePosteriorMatrix>>> = vec![vec![None; n]; n];
        for x in 0..n {
            for y in (x + 1)..n {
                let len_x = seq_lens[x];
                let len_y = seq_lens[y];
                let mut rows = vec![Vec::new(); len_x];

                for (i, out_row) in rows.iter_mut().enumerate() {
                    if let Some(ref p_xy) = cur[x][y]
                        && let Some(r_xy) = p_xy.rows.get(i)
                    {
                        for &(j, p) in r_xy {
                            if j < len_y {
                                if acc[j] == 0.0 {
                                    touched.push(j);
                                }
                                acc[j] += 2.0 * p;
                            }
                        }
                    }
                    for z in 0..n {
                        if z == x || z == y {
                            continue;
                        }
                        let (Some(p_xz), Some(p_zy)) = (&cur[x][z], &cur[z][y]) else {
                            continue;
                        };
                        let Some(r_xz) = p_xz.rows.get(i) else {
                            continue;
                        };
                        for &(k, pxz) in r_xz {
                            let Some(r_zy) = p_zy.rows.get(k) else {
                                continue;
                            };
                            for &(j, pzy) in r_zy {
                                if j < len_y {
                                    if acc[j] == 0.0 {
                                        touched.push(j);
                                    }
                                    acc[j] += pxz * pzy;
                                }
                            }
                        }
                    }
                    touched.sort_unstable();
                    let mut entries = Vec::with_capacity(touched.len());
                    for &j in &touched {
                        let p = (acc[j] * inv_n).clamp(0.0, 1.0);
                        acc[j] = 0.0;
                        if p >= 0.005 {
                            entries.push((j, p));
                        }
                    }
                    touched.clear();
                    *out_row = entries;
                }
                let mat = SparsePosteriorMatrix { len_x, len_y, rows };
                next[y][x] = Some(mat.transpose());
                next[x][y] = Some(mat);
            }
        }
        cur = next;
    }
    cur
}

#[inline]
fn set_bt2(bt: &mut [u8], idx: usize, code: u8) {
    let shift = (idx & 3) * 2;
    bt[idx >> 2] = (bt[idx >> 2] & !(0b11 << shift)) | ((code & 0b11) << shift);
}

#[inline]
fn get_bt2(bt: &[u8], idx: usize) -> u8 {
    (bt[idx >> 2] >> ((idx & 3) * 2)) & 0b11
}

/// Align two profile blocks (`prof_a` over `seqs_a`, `prof_b` over `seqs_b`) by MEA dynamic
/// programming on average triplet-consistent posterior probabilities, using sparse row scattering.
#[allow(clippy::too_many_lines)]
fn merge_profiles_mea(
    prof_a: &[Vec<u8>],
    ids_a: &[usize],
    prof_b: &[Vec<u8>],
    ids_b: &[usize],
    seq_lens: &[usize],
    posteriors: &[Vec<Option<SparsePosteriorMatrix>>],
) -> Vec<Vec<u8>> {
    let wa = prof_a[0].len();
    let wb = prof_b[0].len();
    // Map each profile column of A to each member sequence's 0-based ungapped coordinate.
    let coords_a: Vec<Vec<Option<usize>>> = prof_a
        .iter()
        .map(|row| {
            let mut pos = 0usize;
            row.iter()
                .map(|&b| {
                    if b == b'-' {
                        None
                    } else {
                        let cur = pos;
                        pos += 1;
                        Some(cur)
                    }
                })
                .collect()
        })
        .collect();
    // Precompute reverse map: for each member `rb` in profile B, map ungapped coordinate `pb` -> column `j`.
    let pos_to_col_b: Vec<Vec<usize>> = prof_b
        .iter()
        .zip(ids_b)
        .map(|(row, &sb)| {
            let mut map = vec![usize::MAX; seq_lens.get(sb).copied().unwrap_or(0)];
            let mut pos = 0usize;
            for (col, &b) in row.iter().enumerate() {
                if b != b'-' {
                    if pos < map.len() {
                        map[pos] = col;
                    }
                    pos += 1;
                }
            }
            map
        })
        .collect();

    let stride = wb + 1;
    // Rolling 1D DP rows (`16 bytes * (wb + 1)` instead of `8 * (wa + 1) * (wb + 1)`) plus
    // 2-bit-packed backtrack codes (`4` cells/byte: 1 = diag, 2 = up, 3 = left).
    // Cuts a 20 kb x 20 kb merge from ~3.6 GB down to ~100 MB with bit-exact identical output.
    let mut dp_prev = vec![0.0_f64; stride];
    let mut dp_cur = vec![0.0_f64; stride];
    let mut bt = vec![0u8; ((wa + 1) * stride).div_ceil(4)];
    for i in 1..=wa {
        set_bt2(&mut bt, i * stride, 2);
    }
    for j in 1..=wb {
        set_bt2(&mut bt, j, 3);
    }

    let mut row_scores = vec![0.0_f64; wb];
    let mut touched_cols: Vec<usize> = Vec::with_capacity(64);

    for i in 1..=wa {
        for (ra, &sa) in ids_a.iter().enumerate() {
            let Some(pa) = coords_a[ra][i - 1] else {
                continue;
            };
            for (rb, &sb) in ids_b.iter().enumerate() {
                let Some(ref mat) = posteriors[sa][sb] else {
                    continue;
                };
                let Some(sparse_row) = mat.rows.get(pa) else {
                    continue;
                };
                let col_map = &pos_to_col_b[rb];
                for &(pb, prob) in sparse_row {
                    if let Some(&col_j) = col_map.get(pb)
                        && col_j < wb
                    {
                        if row_scores[col_j] == 0.0 {
                            touched_cols.push(col_j);
                        }
                        row_scores[col_j] += f64::from(prob);
                    }
                }
            }
        }
        let row_i = i * stride;
        dp_cur[0] = 0.0;
        for j in 1..=wb {
            let col_score = row_scores[j - 1];
            let d = dp_prev[j - 1] + col_score;
            let u = dp_prev[j];
            let l = dp_cur[j - 1];
            if d >= u && d >= l {
                dp_cur[j] = d;
                set_bt2(&mut bt, row_i + j, 1);
            } else if u >= l {
                dp_cur[j] = u;
                set_bt2(&mut bt, row_i + j, 2);
            } else {
                dp_cur[j] = l;
                set_bt2(&mut bt, row_i + j, 3);
            }
        }
        std::mem::swap(&mut dp_prev, &mut dp_cur);
        for &c in &touched_cols {
            row_scores[c] = 0.0;
        }
        touched_cols.clear();
    }

    let mut ops = Vec::with_capacity(wa + wb);
    let mut i = wa;
    let mut j = wb;
    while i > 0 || j > 0 {
        let code = get_bt2(&bt, i * stride + j);
        ops.push(code);
        match code {
            1 => {
                i -= 1;
                j -= 1;
            }
            2 => i -= 1,
            _ => j -= 1,
        }
    }
    ops.reverse();

    let mut merged: Vec<Vec<u8>> = (0..(prof_a.len() + prof_b.len()))
        .map(|_| Vec::with_capacity(ops.len()))
        .collect();
    let mut ca = 0usize;
    let mut cb = 0usize;
    for op in ops {
        match op {
            1 => {
                for (r, row) in prof_a.iter().enumerate() {
                    merged[r].push(row[ca]);
                }
                for (r, row) in prof_b.iter().enumerate() {
                    merged[prof_a.len() + r].push(row[cb]);
                }
                ca += 1;
                cb += 1;
            }
            2 => {
                for (r, row) in prof_a.iter().enumerate() {
                    merged[r].push(row[ca]);
                }
                for r in 0..prof_b.len() {
                    merged[prof_a.len() + r].push(b'-');
                }
                ca += 1;
            }
            _ => {
                for row in merged.iter_mut().take(prof_a.len()) {
                    row.push(b'-');
                }
                for (r, row) in prof_b.iter().enumerate() {
                    merged[prof_a.len() + r].push(row[cb]);
                }
                cb += 1;
            }
        }
    }
    merged
}

/// Build a triplet-consistency MEA multiple sequence alignment across `seqs`, preserving every
/// input sequence's ungapped DNA byte-for-byte.
#[must_use]
pub fn align_locus_probcons_msa(
    seqs: &[&[u8]],
    params: &HmmParams,
    consistency_iters: usize,
) -> Vec<Vec<u8>> {
    let n = seqs.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![seqs[0].to_vec()];
    }

    let lens: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
    let mut pw: Vec<Vec<Option<SparsePosteriorMatrix>>> = vec![vec![None; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let mat = SparsePosteriorMatrix::from_pair(seqs[i], seqs[j], params, 0.01);
            pw[j][i] = Some(mat.transpose());
            pw[i][j] = Some(mat);
        }
    }
    let consistent = triplet_consistency_transform(&lens, &pw, consistency_iters);

    // Precompute N x N pairwise expected-accuracy similarity once before greedy merging.
    let mut pair_sim = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if let Some(ref m) = consistent[i][j] {
                let row_max_sum: f64 = m
                    .rows
                    .iter()
                    .map(|r| r.iter().map(|&(_, p)| f64::from(p)).fold(0.0_f64, f64::max))
                    .sum();
                #[allow(clippy::cast_precision_loss)]
                let denom = lens[i].min(lens[j]).max(1) as f64;
                let s = row_max_sum / denom;
                pair_sim[i][j] = s;
                pair_sim[j][i] = s;
            }
        }
    }

    // Progressive MEA merge over guide order (greedy highest expected-accuracy pair first).
    let mut clusters: Vec<(Vec<usize>, Vec<Vec<u8>>)> = seqs
        .iter()
        .enumerate()
        .map(|(i, &s)| (vec![i], vec![s.to_vec()]))
        .collect();

    while clusters.len() > 1 {
        let mut best_pair = (0usize, 1usize);
        let mut best_sim = f64::NEG_INFINITY;
        for a in 0..clusters.len() {
            for b in (a + 1)..clusters.len() {
                let mut sum_acc = 0.0_f64;
                let mut cnt = 0usize;
                for &ia in &clusters[a].0 {
                    for &ib in &clusters[b].0 {
                        sum_acc += pair_sim[ia][ib];
                        cnt += 1;
                    }
                }
                #[allow(clippy::cast_precision_loss)]
                let sim = if cnt > 0 { sum_acc / (cnt as f64) } else { 0.0 };
                if sim > best_sim {
                    best_sim = sim;
                    best_pair = (a, b);
                }
            }
        }
        let (a, b) = best_pair;
        let (ids_b, prof_b) = clusters.remove(b);
        let (mut ids_a, prof_a) = clusters.remove(a);
        let merged = merge_profiles_mea(&prof_a, &ids_a, &prof_b, &ids_b, &lens, &consistent);
        ids_a.extend(ids_b);
        clusters.push((ids_a, merged));
    }

    let (order, rows) = clusters.pop().unwrap_or_default();
    let mut ordered = vec![Vec::new(); n];
    for (orig_idx, row) in order.into_iter().zip(rows) {
        ordered[orig_idx] = row;
    }
    ordered
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triplet_consistency_boosts_transitively_supported_homolog_pair_even_from_upper_triangle() {
        let n = 3usize;
        let len = 20usize;
        let mk_diag = |p: f32| SparsePosteriorMatrix {
            len_x: len,
            len_y: len,
            rows: (0..len).map(|i| vec![(i, p)]).collect(),
        };
        // Populate ONLY upper triangle to verify auto-symmetrization.
        let mut pw: Vec<Vec<Option<SparsePosteriorMatrix>>> = vec![vec![None; n]; n];
        pw[0][1] = Some(mk_diag(0.90));
        pw[1][2] = Some(mk_diag(0.90));
        pw[0][2] = Some(mk_diag(0.20));

        let relaxed = triplet_consistency_transform(&[len, len, len], &pw, 1);
        let p02 = relaxed[0][2].as_ref().expect("relaxed (0,2)").get(5, 5);
        assert!(
            (p02 - (2.0 * 0.20 + 0.90 * 0.90) / 3.0).abs() < 1e-4,
            "expected ProbCons consistency value ~0.4033, got {p02}"
        );
    }

    #[test]
    fn from_pair_handles_empty_and_banded_long_sequences_without_oom() {
        let empty = SparsePosteriorMatrix::from_pair(b"", b"ACGT", &HmmParams::default(), 0.01);
        assert_eq!(empty.len_x, 0);
        assert_eq!(empty.len_y, 4);

        // 1,600 x 1,600 = 2.56M cells > BANDED_THRESHOLD_CELLS (1.5M) -> engages banded path!
        let mut st = 0x1234_5678_9abc_def0_u64;
        let seq: Vec<u8> = (0..1600)
            .map(|_| {
                st = crate::splitmix64(st);
                b"ACGT"[(st & 3) as usize]
            })
            .collect();
        let mat = SparsePosteriorMatrix::from_pair(&seq, &seq, &HmmParams::default(), 0.1);
        let diag_high = (0..1600).filter(|&i| mat.get(i, i) > 0.8).count();
        assert!(
            diag_high > 1500,
            "banded Pair-HMM must recover high diagonal posteriors on identical 1.6 kb pair (got {diag_high}/1600)"
        );
    }

    #[test]
    fn align_locus_probcons_msa_aligns_3_diverged_homologs_and_preserves_sequences() {
        let s0 = b"ACGTACGTACGTTTTTGGGGCCCCAAAATTTTACGTACGTACGT";
        let s1 = b"ACGTACGTACGTTTTTGGGGAAACCAAAATTTTACGTACGTACGT";
        let s2 = b"ACGTACGTACGTTTTTGGGGCCCAAAATTTTACGTACGTACGT";
        let msa = align_locus_probcons_msa(&[s0, s1, s2], &HmmParams::default(), 2);
        assert_eq!(msa.len(), 3);
        let width = msa[0].len();
        assert!(
            width >= s1.len() && width <= s1.len() + 6,
            "unexpected MSA width {width}"
        );
        for (row, orig) in msa
            .iter()
            .zip([s0.as_slice(), s1.as_slice(), s2.as_slice()])
        {
            assert_eq!(row.len(), width);
            let ug: Vec<u8> = row.iter().copied().filter(|&b| b != b'-').collect();
            assert_eq!(ug, orig, "MSA must preserve ungapped DNA byte-for-byte");
        }
    }
}
