//! Six-frame biochemical amino-acid (5-of-6) spaced seeding.

use super::{
    AnchorPairs, ContigIndex, LcbList, append_two_hit_diagonal_lcbs, base_to_idx,
    cluster_strand_lcbs, codon_to_aa,
};
use std::collections::HashMap;

const fn aa_group(ax: u8) -> u8 {
    const GROUP: [u8; 20] = [0, 1, 2, 2, 3, 0, 4, 5, 4, 5, 5, 6, 7, 6, 4, 8, 8, 5, 3, 3];
    if (ax as usize) < 20 {
        GROUP[ax as usize]
    } else {
        9
    }
}

const fn aa_pair_score(ax: u8, ay: u8) -> i32 {
    if ax == ay {
        return 15;
    }
    let gx = aa_group(ax);
    let gy = aa_group(ay);
    if gx < 9 && gx == gy { 6 } else { -6 }
}

fn index_aa_spaced_4of6_frames(
    seq: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) -> HashMap<u32, Vec<usize>> {
    // Reduced 9-group biochemical alphabet spaced 5-of-6 mask ("111011", 18 bp span).
    const AA_OFFSETS: [usize; 5] = [0, 1, 2, 4, 5];
    let mut map: HashMap<u32, Vec<usize>> = HashMap::new();
    let mut counts: HashMap<u32, usize> = HashMap::new();
    for frame in 0..3 {
        if seq.len() < frame + 18 {
            continue;
        }
        let mut pos = frame;
        while pos + 18 <= seq.len() {
            let mut key = 0_u32;
            let mut valid = true;
            for &c in &AA_OFFSETS {
                let s = &seq[pos + c * 3..pos + c * 3 + 3];
                if skip_softmasked && s.iter().any(u8::is_ascii_lowercase) {
                    valid = false;
                    break;
                }
                let Some(aa) = codon_to_aa(s[0], s[1], s[2]) else {
                    valid = false;
                    break;
                };
                key = key * 9 + u32::from(aa_group(aa));
            }
            if valid {
                let cnt = counts.entry(key).or_insert(0);
                *cnt += 1;
                if *cnt <= max_occ {
                    map.entry(key).or_default().push(pos);
                } else {
                    map.remove(&key);
                }
            }
            pos += 3;
        }
    }
    map
}

/// Find spaced `4-of-6` translated amino-acid (18 bp) hits across all reading frames.
pub(crate) fn find_translated_aa5_hits(
    seq_x: &[u8],
    seq_y: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) -> AnchorPairs {
    let ix = index_aa_spaced_4of6_frames(seq_x, max_occ, skip_softmasked);
    let iy = index_aa_spaced_4of6_frames(seq_y, max_occ, skip_softmasked);
    let mut hits = Vec::new();
    for (k, xs) in &ix {
        if let Some(ys) = iy.get(k) {
            for &px in xs {
                for &py in ys {
                    hits.push((px, py));
                }
            }
        }
    }
    hits
}

/// Extend an ungapped codon-level amino-acid HSP from a 5-codon (`15 bp`) seed.
fn extend_aa_ungapped_hsp(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    min_score: i32,
) -> Option<(usize, usize, usize, usize, i32)> {
    const AA_XDROP: i32 = 30;
    if px + 15 > seq_x.len() || py + 15 > seq_y.len() {
        return None;
    }
    let mut core = 0i32;
    for c in 0..5 {
        let ax = codon_to_aa(
            seq_x[px + c * 3],
            seq_x[px + c * 3 + 1],
            seq_x[px + c * 3 + 2],
        )?;
        let ay = codon_to_aa(
            seq_y[py + c * 3],
            seq_y[py + c * 3 + 1],
            seq_y[py + c * 3 + 2],
        )?;
        core += aa_pair_score(ax, ay);
    }

    let mut cur = 0i32;
    let mut best_left = 0i32;
    let (mut lx, mut ly) = (px, py);
    let (mut ix, mut iy) = (px, py);
    while ix >= 3 && iy >= 3 {
        ix -= 3;
        iy -= 3;
        let (Some(ax), Some(ay)) = (
            codon_to_aa(seq_x[ix], seq_x[ix + 1], seq_x[ix + 2]),
            codon_to_aa(seq_y[iy], seq_y[iy + 1], seq_y[iy + 2]),
        ) else {
            break;
        };
        cur += aa_pair_score(ax, ay);
        if cur > best_left {
            best_left = cur;
            lx = ix;
            ly = iy;
        } else if best_left - cur > AA_XDROP {
            break;
        }
    }

    cur = 0;
    let mut best_right = 0i32;
    let (mut rx, mut ry) = (px + 15, py + 15);
    ix = px + 15;
    iy = py + 15;
    while ix + 3 <= seq_x.len() && iy + 3 <= seq_y.len() {
        let (Some(ax), Some(ay)) = (
            codon_to_aa(seq_x[ix], seq_x[ix + 1], seq_x[ix + 2]),
            codon_to_aa(seq_y[iy], seq_y[iy + 1], seq_y[iy + 2]),
        ) else {
            break;
        };
        cur += aa_pair_score(ax, ay);
        ix += 3;
        iy += 3;
        if cur > best_right {
            best_right = cur;
            rx = ix;
            ry = iy;
        } else if best_right - cur > AA_XDROP {
            break;
        }
    }

    let total = core + best_left + best_right;
    (total >= min_score).then_some((lx, ly, rx, ry, total))
}

pub(crate) fn append_hsp_gated_aa_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    hits: &[(usize, usize)],
    min_aa_score: i32,
) {
    let mut gated = Vec::new();
    let mut strong = Vec::new();
    for &(px, py) in hits {
        if let Some((lx, ly, rx, ry, sc)) =
            extend_aa_ungapped_hsp(seq_x, seq_y, px, py, min_aa_score)
        {
            gated.push((px, py));
            if sc >= min_aa_score + 30 && rx >= lx + 24 && ry >= ly + 24 {
                strong.push(vec![(lx, ly), (rx - 15, ry - 15)]);
            }
        }
    }
    let mut candidates = cluster_strand_lcbs(&gated, 1800, 200, 2, 24, 250.0);
    candidates.extend(strong);
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if !target
            .iter()
            .any(|c| sx <= c[c.len() - 1].0 && ex >= c[0].0)
        {
            target.push(lcb);
        }
    }
}

/// Find 2-hit translated amino-acid `5-mer` (15 bp) anchors across all reading frames.
#[allow(dead_code)]
pub(crate) fn append_translated_aa5_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) {
    let hits = find_translated_aa5_hits(seq_x, seq_y, max_occ, skip_softmasked);
    append_two_hit_diagonal_lcbs(target, &hits, 120, 12);
}

#[inline]
fn hsp_pair_score(a: u8, b: u8) -> i32 {
    let (Some(ia), Some(ib)) = (base_to_idx(a), base_to_idx(b)) else {
        return -15;
    };
    if ia == ib {
        10
    } else if (ia ^ ib) == 2 {
        -3 // transition (A<->G, C<->T)
    } else {
        -12 // transversion
    }
}

/// Extend a 19-bp seed hit `(px, py)` without gaps using transition-aware X-drop;
/// returns `Some((start_x, start_y, end_x, end_y, score))` when `score >= min_score`.
fn extend_ungapped_hsp(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    min_score: i32,
) -> Option<(usize, usize, usize, usize, i32)> {
    const XDROP: i32 = 60;
    if px + 19 > seq_x.len() || py + 19 > seq_y.len() {
        return None;
    }
    let mut core = 0i32;
    for off in 0..19 {
        core += hsp_pair_score(seq_x[px + off], seq_y[py + off]);
    }

    let mut cur = 0i32;
    let mut best_left = 0i32;
    let mut lx = px;
    let mut ly = py;
    let mut ix = px;
    let mut iy = py;
    while ix > 0 && iy > 0 {
        ix -= 1;
        iy -= 1;
        cur += hsp_pair_score(seq_x[ix], seq_y[iy]);
        if cur > best_left {
            best_left = cur;
            lx = ix;
            ly = iy;
        } else if best_left - cur > XDROP {
            break;
        }
    }

    cur = 0;
    let mut best_right = 0i32;
    let mut rx = px + 19;
    let mut ry = py + 19;
    ix = px + 19;
    iy = py + 19;
    while ix < seq_x.len() && iy < seq_y.len() {
        cur += hsp_pair_score(seq_x[ix], seq_y[iy]);
        ix += 1;
        iy += 1;
        if cur > best_right {
            best_right = cur;
            rx = ix;
            ry = iy;
        } else if best_right - cur > XDROP {
            break;
        }
    }

    let total = core + best_left + best_right;
    (total >= min_score).then_some((lx, ly, rx, ry, total))
}

pub(crate) fn append_hsp_gated_lastz12_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    hits: &[(usize, usize)],
    min_hsp_score: i32,
) {
    #[allow(clippy::cast_precision_loss)]
    let space_bits = ((seq_x.len().max(1) as f64) * (seq_y.len().max(1) as f64))
        .log2()
        .max(10.0);
    #[allow(clippy::cast_possible_truncation)]
    let strong_margin = (1.8 * space_bits).round() as i32;
    let mut filtered_hits: Vec<(usize, usize)> = Vec::new();
    let mut strong_hsps: Vec<Vec<(usize, usize)>> = Vec::new();
    for &(px, py) in hits {
        if let Some((lx, ly, rx, ry, sc)) = extend_ungapped_hsp(seq_x, seq_y, px, py, min_hsp_score)
        {
            filtered_hits.push((px, py));
            if sc >= min_hsp_score + strong_margin && rx >= lx + 18 && ry >= ly + 18 {
                strong_hsps.push(vec![(lx, ly), (rx - 12, ry - 12)]);
            }
        }
    }
    let mut candidates = cluster_strand_lcbs(&filtered_hits, 1800, 200, 3, 30, 250.0);
    candidates.extend(cluster_strand_lcbs(&filtered_hits, 150, 25, 2, 16, 150.0));
    candidates.extend(strong_hsps);
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if !target
            .iter()
            .any(|c| sx <= c[c.len() - 1].0 && ex >= c[0].0)
        {
            target.push(lcb);
        }
    }
}

pub(crate) fn find_w12_raw_lcbs(idx_x: &ContigIndex, idx_y: &ContigIndex) -> (LcbList, LcbList) {
    let (hits_fwd, hits_rc) = idx_x.seeds.intersect(&idx_y.seeds);
    let raw_fwd = cluster_strand_lcbs(&hits_fwd, 2500, 250, 7, 70, 220.0);
    let raw_rc = cluster_strand_lcbs(&hits_rc, 2500, 250, 7, 70, 220.0);
    (raw_fwd, raw_rc)
}

pub(crate) fn append_nonoverlapping_lcbs(
    target: &mut LcbList,
    existing: &[Vec<(usize, usize)>],
    hits: &[(usize, usize)],
    min_anchors: usize,
    min_span: usize,
) {
    let mut candidates = cluster_strand_lcbs(hits, 1800, 200, min_anchors, min_span, 250.0);
    // Tight 2-hit diagonal pass (LASTZ-style): 2 nearby hits within 120bp and <=25bp drift.
    candidates.extend(cluster_strand_lcbs(hits, 120, 25, 2, 16, 120.0));
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if !existing
            .iter()
            .chain(target.iter())
            .any(|c| sx <= c[c.len() - 1].0 && ex >= c[0].0)
        {
            target.push(lcb);
        }
    }
}
