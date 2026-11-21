//! Compact flat sorted spaced-seed index with optional 16-bit prefix bucketing.

use super::AnchorPairs;
use crate::AlignError;

/// Compact flat sorted spaced-seed index (`12 bytes/entry` vs `~100 bytes/entry` for `HashMap<u64, Vec<usize>>`).
///
/// Since canonical seeds (`fwd.min(rc)`) are strand-symmetric, reverse-complement positions are
/// derived directly as `seq_len - span - pos`, halving memory and build time.
#[derive(Debug, Clone)]
pub(crate) struct FlatSeedIndex {
    pub(crate) keys: Vec<u64>,
    pub(crate) pos: Vec<u32>,
    pub(crate) span: usize,
    pub(crate) seq_len: usize,
    pub(super) buckets: Vec<u32>,
    bucket_shift: u8,
}

impl FlatSeedIndex {
    #[allow(clippy::too_many_lines)]
    pub(crate) fn build_with_sdust(
        seq: &[u8],
        mask: &str,
        max_occ: usize,
        skip_softmasked: bool,
        ry_mode: bool,
        sdust: Option<&[bool]>,
    ) -> Result<Self, AlignError> {
        let (span, offsets) = super::parse_mask_offsets(mask, ry_mode)?;
        if seq.len() < span {
            return Ok(Self {
                keys: Vec::new(),
                pos: Vec::new(),
                span,
                seq_len: seq.len(),
                buckets: Vec::new(),
                bucket_shift: 0,
            });
        }
        if seq.len() > u32::MAX as usize {
            return Err(AlignError::InvalidInput(
                "sequence length exceeds 4 Gbp u32 index limit".to_owned(),
            ));
        }

        if let Some(m) = sdust
            && m.len() < seq.len()
        {
            return Err(AlignError::InvalidInput(
                "sdust mask shorter than sequence".to_owned(),
            ));
        }
        let rev_offsets: Vec<usize> = offsets.iter().map(|&o| span - 1 - o).collect();
        let mut pairs: Vec<(u64, u32)> = Vec::with_capacity((seq.len() - span + 1).min(1 << 24));

        for i in 0..=(seq.len() - span) {
            if super::is_masked_window(seq, sdust, i, span) {
                continue;
            }
            let sub = &seq[i..i + span];
            let Some((fwd, rc)) =
                super::encode_spaced_pair(sub, &offsets, &rev_offsets, skip_softmasked, ry_mode)
            else {
                continue;
            };
            #[allow(clippy::cast_possible_truncation)]
            pairs.push((fwd.min(rc), i as u32));
        }

        pairs.sort_unstable();

        let mut keys = Vec::with_capacity(pairs.len());
        let mut pos = Vec::with_capacity(pairs.len());
        let mut i = 0;
        while i < pairs.len() {
            let k = pairs[i].0;
            let mut j = i + 1;
            while j < pairs.len() && pairs[j].0 == k {
                j += 1;
            }
            if j - i <= max_occ {
                for &(key, p) in &pairs[i..j] {
                    keys.push(key);
                    pos.push(p);
                }
            }
            i = j;
        }
        drop(pairs);
        keys.shrink_to_fit();
        pos.shrink_to_fit();

        let (buckets, bucket_shift) = if keys.len() >= 32_768 && u32::try_from(keys.len()).is_ok() {
            let max_key = *keys.last().unwrap_or(&0);
            let max_bits = 64 - max_key.leading_zeros();
            let shift = u8::try_from(max_bits.saturating_sub(16)).unwrap_or(0);
            let mut bkt = vec![0u32; 65_537];
            let mut cur_b = 0usize;
            for (idx, &k) in keys.iter().enumerate() {
                // Fallbacks are unreachable given `shift` and `u32::try_from(keys.len()).is_ok()`,
                // kept to satisfy strict cast lints without `as`.
                let b = usize::try_from(k >> shift).unwrap_or(65_535).min(65_535);
                let idx_u32 = u32::try_from(idx).unwrap_or(u32::MAX);
                while cur_b <= b {
                    bkt[cur_b] = idx_u32;
                    cur_b += 1;
                }
            }
            let total_u32 = u32::try_from(keys.len()).unwrap_or(u32::MAX);
            while cur_b <= 65_536 {
                bkt[cur_b] = total_u32;
                cur_b += 1;
            }
            (bkt, shift)
        } else {
            (Vec::new(), 0)
        };

        Ok(Self {
            keys,
            pos,
            span,
            seq_len: seq.len(),
            buckets,
            bucket_shift,
        })
    }

    pub(crate) fn get(&self, key: u64) -> &[u32] {
        let (lo, hi) = if self.buckets.is_empty() {
            (0, self.keys.len())
        } else {
            let Ok(b) = usize::try_from(key >> self.bucket_shift) else {
                return &[];
            };
            if b >= 65_536 {
                return &[];
            }
            (self.buckets[b] as usize, self.buckets[b + 1] as usize)
        };
        let slice = &self.keys[lo..hi];
        let rel = slice.partition_point(|&k| k < key);
        let start = lo + rel;
        if start >= hi || self.keys[start] != key {
            return &[];
        }
        let mut end = start + 1;
        while end < hi && self.keys[end] == key {
            end += 1;
        }
        &self.pos[start..end]
    }

    pub(crate) fn shared_unique_keys(&self, other: &Self) -> usize {
        let (small, large) = if self.keys.len() <= other.keys.len() {
            (self, other)
        } else {
            (other, self)
        };
        if small.keys.is_empty() || large.keys.is_empty() {
            return 0;
        }
        if small.keys.len().saturating_mul(8) < large.keys.len() {
            let mut count = 0;
            let mut i = 0;
            while i < small.keys.len() {
                let k = small.keys[i];
                if !large.get(k).is_empty() {
                    count += 1;
                }
                while i < small.keys.len() && small.keys[i] == k {
                    i += 1;
                }
            }
            return count;
        }
        let mut i = 0;
        let mut j = 0;
        let mut count = 0;
        while i < self.keys.len() && j < other.keys.len() {
            let ka = self.keys[i];
            let kb = other.keys[j];
            match ka.cmp(&kb) {
                std::cmp::Ordering::Less => {
                    while i < self.keys.len() && self.keys[i] == ka {
                        i += 1;
                    }
                }
                std::cmp::Ordering::Greater => {
                    while j < other.keys.len() && other.keys[j] == kb {
                        j += 1;
                    }
                }
                std::cmp::Ordering::Equal => {
                    count += 1;
                    while i < self.keys.len() && self.keys[i] == ka {
                        i += 1;
                    }
                    while j < other.keys.len() && other.keys[j] == kb {
                        j += 1;
                    }
                }
            }
        }
        count
    }

    pub(crate) fn intersect(&self, other: &Self) -> (AnchorPairs, AnchorPairs) {
        let mut hits_fwd = Vec::new();
        let mut hits_rc = Vec::new();
        let mut i = 0;
        let mut j = 0;
        while i < self.keys.len() && j < other.keys.len() {
            let ka = self.keys[i];
            let kb = other.keys[j];
            match ka.cmp(&kb) {
                std::cmp::Ordering::Less => {
                    while i < self.keys.len() && self.keys[i] == ka {
                        i += 1;
                    }
                }
                std::cmp::Ordering::Greater => {
                    while j < other.keys.len() && other.keys[j] == kb {
                        j += 1;
                    }
                }
                std::cmp::Ordering::Equal => {
                    let i_start = i;
                    while i < self.keys.len() && self.keys[i] == ka {
                        i += 1;
                    }
                    let j_start = j;
                    while j < other.keys.len() && other.keys[j] == kb {
                        j += 1;
                    }
                    for &px in &self.pos[i_start..i] {
                        let x_pos = px as usize;
                        for &py in &other.pos[j_start..j] {
                            let y_pos = py as usize;
                            hits_fwd.push((x_pos, y_pos));
                            hits_rc.push((x_pos, other.seq_len - other.span - y_pos));
                        }
                    }
                }
            }
        }
        (hits_fwd, hits_rc)
    }
}
