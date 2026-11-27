//! Two-pass Fitch parsimony (Fitch 1971) ancestral sequence reconstruction with nibble packing.

pub(super) const FITCH_ALL: u8 = 0b1111;

pub(super) const fn fitch_base_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' => 0b1000,
        _ => FITCH_ALL,
    }
}

#[inline]
pub(super) const fn fitch_step(left: u8, right: u8) -> u8 {
    let inter = left & right;
    if inter != 0 { inter } else { left | right }
}

pub(super) fn add_consensus_base(counts: &mut [u16; 4], base: u8) {
    match base.to_ascii_uppercase() {
        b'A' => counts[0] = counts[0].saturating_add(1),
        b'C' => counts[1] = counts[1].saturating_add(1),
        b'G' => counts[2] = counts[2].saturating_add(1),
        b'T' => counts[3] = counts[3].saturating_add(1),
        _ => {}
    }
}

#[inline]
pub(super) fn pack_fitch_masks_into(masks: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.reserve(masks.len().div_ceil(2));
    let (chunks, rem) = masks.as_chunks::<2>();
    for &[a, b] in chunks {
        dst.push((a & 0x0F) | ((b & 0x0F) << 4));
    }
    if let Some(&last) = rem.first() {
        dst.push(last & 0x0F);
    }
}

#[inline]
pub(super) fn pack_fitch_masks(masks: &[u8]) -> Vec<u8> {
    let mut packed = Vec::with_capacity(masks.len().div_ceil(2));
    pack_fitch_masks_into(masks, &mut packed);
    packed
}

#[inline]
fn unpack_fitch_mask(packed: &[u8], idx: usize) -> u8 {
    let byte = packed[idx >> 1];
    if (idx & 1) == 0 {
        byte & 0x0F
    } else {
        (byte >> 4) & 0x0F
    }
}

/// Resolve root ancestral bases from nibble-packed Fitch state sets (`root_sets`), breaking ties
/// by total aligned leaf counts (`leaf_counts`) and then `ref_dna`.
pub(super) fn fitch_resolve_root(
    ref_dna: &[u8],
    root_sets: &[u8],
    leaf_counts: &[[u16; 4]],
    out: &mut [u8],
) {
    for (pos, ((&rb, counts), dst)) in ref_dna
        .iter()
        .zip(leaf_counts.iter())
        .zip(out.iter_mut())
        .enumerate()
    {
        let mask = unpack_fitch_mask(root_sets, pos);
        let ref_b = rb.to_ascii_uppercase();
        let active = if mask == 0 { FITCH_ALL } else { mask };
        let mut best_b = rb;
        let mut best_cnt = 0u16;
        let mut found = false;
        for (idx, (&b, &cnt)) in b"ACGT".iter().zip(counts.iter()).enumerate() {
            if (active & (1u8 << idx)) == 0 {
                continue;
            }
            if !found || cnt > best_cnt || (cnt == best_cnt && b == ref_b) {
                found = true;
                best_cnt = cnt;
                best_b = b;
            }
        }
        *dst = if found && best_cnt > 0 { best_b } else { rb };
    }
}

/// Top-down Fitch parsimony pass: inherit `parent_dna[i]` whenever it lies in
/// nibble-packed `child_sets[i]` (zero substitutions on the branch), otherwise pick from
/// `child_sets[i]` preferring `ref_dna[i]` on ties.
pub(super) fn fitch_resolve_child(
    ref_dna: &[u8],
    parent_dna: &[u8],
    child_sets: &[u8],
    out: &mut [u8],
) {
    for (pos, ((&rb, &pb), dst)) in ref_dna
        .iter()
        .zip(parent_dna.iter())
        .zip(out.iter_mut())
        .enumerate()
    {
        let mask = unpack_fitch_mask(child_sets, pos);
        let p_mask = fitch_base_mask(pb);
        if p_mask == FITCH_ALL && (mask == FITCH_ALL || mask == 0) {
            *dst = pb;
            continue;
        }
        if p_mask != FITCH_ALL && (mask & p_mask) != 0 {
            *dst = pb;
            continue;
        }
        let ref_b = rb.to_ascii_uppercase();
        let r_mask = fitch_base_mask(ref_b);
        if r_mask != FITCH_ALL && (mask & r_mask) != 0 {
            *dst = ref_b;
            continue;
        }
        let mut chosen = rb;
        for (idx, &b) in b"ACGT".iter().enumerate() {
            if (mask & (1u8 << idx)) != 0 {
                chosen = b;
                break;
            }
        }
        *dst = chosen;
    }
}
