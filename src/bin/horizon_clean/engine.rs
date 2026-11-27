use anyhow::Result;
use horizon_aligner::bocpd::{BocpdScratch, bocpd_gc_profile_into};
use horizon_aligner::gfa::{
    GfaPathRef, GfaRecordRef, GfaWalkRef, parse_cigar_overlap_bp, parse_path_sample_fields,
    parse_path_steps, parse_walk_steps, stream_gfa_record_refs,
};
use horizon_aligner::seeding::revcomp;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::ops::ControlFlow;

const MAX_PACKED_SEG_ID: u32 = (1u32 << 31) - 1;
const MIN_JUNCTION_FLANK_BP: usize = 12;

/// Compact 8-byte directed GFA edge `(from, from_fwd) -> (to, to_fwd)` packed into a `u64`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct PackedEdge(u64);

impl PackedEdge {
    #[inline]
    fn new(from: u32, from_fwd: bool, to: u32, to_fwd: bool) -> Self {
        let a = (from << 1) | u32::from(from_fwd);
        let b = (to << 1) | u32::from(to_fwd);
        Self((u64::from(a) << 32) | u64::from(b))
    }

    #[inline]
    fn src_id(self) -> u32 {
        (self.0 >> 33) as u32
    }

    #[inline]
    fn src_fwd(self) -> bool {
        ((self.0 >> 32) & 1) != 0
    }

    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    fn dst_id(self) -> u32 {
        (self.0 as u32) >> 1
    }

    #[inline]
    fn dst_fwd(self) -> bool {
        (self.0 & 1) != 0
    }

    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    fn rev(self) -> Self {
        let a = (self.0 >> 32) as u32;
        let b = self.0 as u32;
        Self((u64::from(b ^ 1) << 32) | u64::from(a ^ 1))
    }

    #[inline]
    fn canon(self) -> Self {
        self.min(self.rev())
    }
}

/// Segment-name interner with a zero-`String` `u32` fast path for canonical decimal GFA segment
/// IDs (standard in variation GFAs) and a `String` fallback for arbitrary alphanumeric names.
#[derive(Default)]
struct SegInterner {
    num_ids: HashMap<u32, u32>,
    str_ids: HashMap<String, u32>,
    next_id: u32,
}

impl SegInterner {
    #[inline]
    fn parse_canonical_u32(name: &str) -> Option<u32> {
        let bytes = name.as_bytes();
        if bytes.is_empty() || bytes.len() > 10 {
            return None;
        }
        if bytes.len() > 1 && bytes[0] == b'0' {
            return None;
        }
        if !bytes.iter().all(u8::is_ascii_digit) {
            return None;
        }
        name.parse::<u32>().ok()
    }

    #[inline]
    fn get(&self, name: &str) -> Option<u32> {
        if let Some(num) = Self::parse_canonical_u32(name) {
            self.num_ids.get(&num).copied()
        } else {
            self.str_ids.get(name).copied()
        }
    }

    fn intern(&mut self, name: &str) -> Result<u32> {
        if let Some(num) = Self::parse_canonical_u32(name) {
            if let Some(&id) = self.num_ids.get(&num) {
                return Ok(id);
            }
            let id = self.alloc_id()?;
            self.num_ids.insert(num, id);
            return Ok(id);
        }
        if let Some(&id) = self.str_ids.get(name) {
            return Ok(id);
        }
        let id = self.alloc_id()?;
        self.str_ids.insert(name.to_owned(), id);
        Ok(id)
    }

    fn alloc_id(&mut self) -> Result<u32> {
        let id = self.next_id;
        anyhow::ensure!(
            id < MAX_PACKED_SEG_ID,
            "segment interner exceeded {MAX_PACKED_SEG_ID} segments"
        );
        self.next_id = id + 1;
        Ok(id)
    }

    #[inline]
    fn len(&self) -> usize {
        self.next_id as usize
    }
}

/// Compact arena slice storing prefix and suffix flanks for segments with `len >= 12`.
#[derive(Clone, Copy)]
struct FlankRef {
    pre_off: u32,
    suf_off: u32,
    pre_len: u16,
    suf_len: u16,
}

#[derive(Default)]
struct FlankStore {
    arena: Vec<u8>,
    map: HashMap<u32, FlankRef>,
}

impl FlankStore {
    fn insert(&mut self, seg_id: u32, seq: &[u8], flank_bp: usize) {
        if seq.len() < MIN_JUNCTION_FLANK_BP {
            return;
        }
        let take = seq.len().min(flank_bp).min(usize::from(u16::MAX));
        let Ok(pre_off) = u32::try_from(self.arena.len()) else {
            return;
        };
        #[allow(clippy::cast_possible_truncation)]
        let take_u16 = take as u16;
        if seq.len() <= take {
            self.arena.extend_from_slice(&seq[..take]);
            self.map.insert(
                seg_id,
                FlankRef {
                    pre_off,
                    suf_off: pre_off,
                    pre_len: take_u16,
                    suf_len: take_u16,
                },
            );
        } else {
            self.arena.extend_from_slice(&seq[..take]);
            let Ok(suf_off) = u32::try_from(self.arena.len()) else {
                return;
            };
            self.arena.extend_from_slice(&seq[seq.len() - take..]);
            self.map.insert(
                seg_id,
                FlankRef {
                    pre_off,
                    suf_off,
                    pre_len: take_u16,
                    suf_len: take_u16,
                },
            );
        }
    }

    fn oriented_3p(&self, seg_id: u32, is_fwd: bool) -> Option<std::borrow::Cow<'_, [u8]>> {
        let f = self.map.get(&seg_id)?;
        if is_fwd {
            let s = f.suf_off as usize;
            Some(std::borrow::Cow::Borrowed(
                &self.arena[s..s + usize::from(f.suf_len)],
            ))
        } else {
            let s = f.pre_off as usize;
            Some(std::borrow::Cow::Owned(revcomp(
                &self.arena[s..s + usize::from(f.pre_len)],
            )))
        }
    }

    fn oriented_5p(&self, seg_id: u32, is_fwd: bool) -> Option<std::borrow::Cow<'_, [u8]>> {
        let f = self.map.get(&seg_id)?;
        if is_fwd {
            let s = f.pre_off as usize;
            Some(std::borrow::Cow::Borrowed(
                &self.arena[s..s + usize::from(f.pre_len)],
            ))
        } else {
            let s = f.suf_off as usize;
            Some(std::borrow::Cow::Owned(revcomp(
                &self.arena[s..s + usize::from(f.suf_len)],
            )))
        }
    }
}

/// Memory-efficient distinct haplotype support tracker per segment ID.
///
/// For `min_haplotypes <= 2` (the standard pangenome case), tracks exact `0 / 1 / >=2` distinct
/// haplotypes in a single flat `Vec<u32>` (`4 B/segment`, zero per-segment `HashSet` allocations).
enum HapTracker {
    Disabled,
    /// `u32::MAX` = 0 haps, `< u32::MAX - 1` = single hap ID, `u32::MAX - 1` = `>= 2` haps.
    UpToTwo(Vec<u32>),
    /// General `k > 2`: tracks exact distinct count up to `min_haps` and frees the small `Vec`
    /// once a segment reaches `min_haps`.
    General {
        min_haps: usize,
        counts: Vec<u16>,
        sets: HashMap<u32, Vec<u32>>,
    },
}

const HAP_NONE: u32 = u32::MAX;
const HAP_MULTI: u32 = u32::MAX - 1;

impl HapTracker {
    fn new(min_haplotypes: usize) -> Self {
        match min_haplotypes {
            0 => Self::Disabled,
            1 | 2 => Self::UpToTwo(Vec::new()),
            k => Self::General {
                min_haps: k,
                counts: Vec::new(),
                sets: HashMap::new(),
            },
        }
    }

    #[inline]
    fn record(&mut self, seg_id: u32, hap_id: u32) {
        let idx = seg_id as usize;
        match self {
            Self::Disabled => {}
            Self::UpToTwo(states) => {
                if idx >= states.len() {
                    states.resize(idx + 1, HAP_NONE);
                }
                let cur = states[idx];
                if cur == HAP_NONE {
                    states[idx] = hap_id.min(HAP_MULTI - 1);
                } else if cur != HAP_MULTI && cur != hap_id.min(HAP_MULTI - 1) {
                    states[idx] = HAP_MULTI;
                }
            }
            Self::General {
                min_haps,
                counts,
                sets,
            } => {
                if idx >= counts.len() {
                    counts.resize(idx + 1, 0);
                }
                if usize::from(counts[idx]) >= *min_haps {
                    return;
                }
                let list = sets.entry(seg_id).or_default();
                if !list.contains(&hap_id) {
                    list.push(hap_id);
                    counts[idx] = counts[idx].saturating_add(1);
                    if usize::from(counts[idx]) >= *min_haps {
                        sets.remove(&seg_id);
                    }
                }
            }
        }
    }

    #[inline]
    fn meets_min(&self, seg_id: u32, min_haplotypes: usize) -> bool {
        let idx = seg_id as usize;
        match self {
            Self::Disabled => true,
            Self::UpToTwo(states) => {
                let st = states.get(idx).copied().unwrap_or(HAP_NONE);
                if min_haplotypes <= 1 {
                    st != HAP_NONE
                } else {
                    st == HAP_MULTI
                }
            }
            Self::General { counts, .. } => {
                usize::from(counts.get(idx).copied().unwrap_or(0)) >= min_haplotypes
            }
        }
    }
}

fn gc_counts(seq: &[u8]) -> (f64, f64) {
    let mut gc = 0.0_f64;
    let mut at = 0.0_f64;
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => gc += 1.0,
            b'A' | b'T' => at += 1.0,
            _ => {}
        }
    }
    (gc, gc + at)
}

fn bernoulli_ll_bits(k: f64, n: f64) -> f64 {
    if n <= 0.0 || k <= 0.0 || k >= n {
        return 0.0;
    }
    let p = (k / n).clamp(1e-6, 1.0 - 1e-6);
    k * p.log2() + (n - k) * (1.0 - p).log2()
}

/// Score a GFA link junction using MDL two-regime changepoint gain (bits) gated by BOCPD GC shift.
/// Positive bits indicate a statistically significant compositional break at the junction.
fn score_junction_cp_bits_with_scratch(
    flank_u: &[u8],
    flank_v: &[u8],
    scratch: &mut BocpdScratch,
    concat: &mut Vec<u8>,
    gc_prof: &mut Vec<f64>,
) -> f64 {
    if flank_u.len() < MIN_JUNCTION_FLANK_BP || flank_v.len() < MIN_JUNCTION_FLANK_BP {
        return 0.0;
    }
    let (ku, nu) = gc_counts(flank_u);
    let (kv, nv) = gc_counts(flank_v);
    let n_tot = nu + nv;
    if nu < 10.0 || nv < 10.0 || n_tot < 24.0 {
        return 0.0;
    }
    let ll_null = bernoulli_ll_bits(ku + kv, n_tot);
    let ll_alt = bernoulli_ll_bits(ku, nu) + bernoulli_ll_bits(kv, nv);
    let mdl_penalty = 1.5 * n_tot.log2();
    let cp_bits = (ll_alt - ll_null) - mdl_penalty;
    if cp_bits <= 0.0 {
        return cp_bits;
    }

    let bocpd_jump_one_way =
        |a: &[u8], b: &[u8], sc: &mut BocpdScratch, cat: &mut Vec<u8>, gp: &mut Vec<f64>| -> f64 {
            cat.clear();
            cat.extend_from_slice(a);
            cat.extend_from_slice(b);
            bocpd_gc_profile_into(cat, 1.0 / 120.0, 80, sc, gp);
            let split = a.len();
            let left_gc = gp[split - 1];
            let right_gc = gp[(split + 15).min(gp.len() - 1)];
            (right_gc - left_gc).abs()
        };
    // Evaluate BOCPD jump symmetrically in both edge orientations (`u++v` and `rc(v)++rc(u)`) so
    // junction severing is strictly invariant to GFA `S`-line input order and edge strand.
    let jump_fwd = bocpd_jump_one_way(flank_u, flank_v, scratch, concat, gc_prof);
    let rev_v = revcomp(flank_v);
    let rev_u = revcomp(flank_u);
    let jump_rev = bocpd_jump_one_way(&rev_v, &rev_u, scratch, concat, gc_prof);
    let bocpd_jump = jump_fwd.max(jump_rev);
    if bocpd_jump < 0.10 { 0.0 } else { cp_bits }
}

pub(crate) struct CleanConfig {
    pub(crate) flank_bp: usize,
    pub(crate) min_cp_bits: f64,
    pub(crate) prune_unsupported: bool,
    pub(crate) min_haplotypes: usize,
    pub(crate) clip_bubbles_bp: usize,
}

pub(crate) struct CleanStats {
    pub(crate) total_links: usize,
    pub(crate) kept_links: usize,
    pub(crate) severed_links: usize,
    pub(crate) dropped_links: usize,
    pub(crate) emitted_bubbles: usize,
    pub(crate) split_walks: usize,
    pub(crate) dropped_low_hap_segs: usize,
    pub(crate) clipped_bubbles: usize,
}

struct CleanPass1 {
    seg_ids: SegInterner,
    seg_lens: Vec<u32>,
    seg_seen_s: Vec<bool>,
    flanks: FlankStore,
    edges: Vec<PackedEdge>,
    /// Canonical edge -> walk support saturated at `2u8`.
    walk_support: HashMap<PackedEdge, u8>,
    link_overlaps: HashMap<PackedEdge, usize>,
    hap_tracker: HapTracker,
    horizon_bubbles: HashSet<u32>,
}

fn record_walk_support(
    steps: &[(&str, bool)],
    hap_id: Option<u32>,
    seg_ids: &mut SegInterner,
    walk_support: &mut HashMap<PackedEdge, u8>,
    hap_tracker: &mut HapTracker,
) -> Result<()> {
    let mut prev: Option<(u32, bool)> = None;
    for &(name, fwd) in steps {
        let id = seg_ids.intern(name)?;
        if let Some(hid) = hap_id {
            hap_tracker.record(id, hid);
        }
        if let Some((u, u_fwd)) = prev {
            let canon = PackedEdge::new(u, u_fwd, id, fwd).canon();
            let entry = walk_support.entry(canon).or_insert(0);
            if *entry < 2 {
                *entry += 1;
            }
        }
        prev = Some((id, fwd));
    }
    Ok(())
}

fn has_horizon_bubble_tag(tags: Option<&str>) -> bool {
    tags.is_some_and(|raw| raw.split('\t').any(|tok| tok == "HB:i:1"))
}

fn intern_hap_id(
    hap_ids: &mut HashMap<(String, String), u32>,
    next_hap_id: &mut u32,
    sample: &str,
    haplotype: &str,
) -> u32 {
    if let Some(&id) = hap_ids.get(&(sample.to_owned(), haplotype.to_owned())) {
        return id;
    }
    let id = *next_hap_id;
    *next_hap_id = next_hap_id.saturating_add(1);
    hap_ids.insert((sample.to_owned(), haplotype.to_owned()), id);
    id
}

fn collect_clean_pass1(
    gfa_path: impl AsRef<std::path::Path>,
    cfg: &CleanConfig,
) -> Result<CleanPass1> {
    let gfa_path = gfa_path.as_ref();
    let mut seg_ids = SegInterner::default();
    let mut seg_lens: Vec<u32> = Vec::new();
    let mut seg_seen_s: Vec<bool> = Vec::new();
    let mut flanks = FlankStore::default();
    let mut edges = Vec::<PackedEdge>::new();
    let mut walk_support = HashMap::<PackedEdge, u8>::new();
    let mut link_overlaps = HashMap::<PackedEdge, usize>::new();
    let mut hap_tracker = HapTracker::new(cfg.min_haplotypes);
    let mut hap_ids = HashMap::<(String, String), u32>::new();
    let mut next_hap_id = 0_u32;
    let mut horizon_bubbles = HashSet::<u32>::new();
    let track_haps = cfg.min_haplotypes > 0;
    let track_bubbles = cfg.clip_bubbles_bp > 0;

    stream_gfa_record_refs(gfa_path, |rec| {
        match rec {
            GfaRecordRef::Segment(s) => {
                let id = seg_ids.intern(s.name)?;
                let idx = id as usize;
                if idx >= seg_lens.len() {
                    seg_lens.resize(idx + 1, 0);
                    seg_seen_s.resize(idx + 1, false);
                }
                let has_seq = s.sequence != b"*";
                let len = s.seq_len();
                seg_lens[idx] = u32::try_from(len).unwrap_or(u32::MAX);
                seg_seen_s[idx] = true;
                if has_seq && len >= MIN_JUNCTION_FLANK_BP {
                    flanks.insert(id, s.sequence, cfg.flank_bp);
                }
                if track_bubbles && has_horizon_bubble_tag(s.tags) {
                    horizon_bubbles.insert(id);
                }
            }
            GfaRecordRef::Link(l) => {
                let u = seg_ids.intern(l.from)?;
                let v = seg_ids.intern(l.to)?;
                let e = PackedEdge::new(u, l.from_orient != b'-', v, l.to_orient != b'-');
                let ov = parse_cigar_overlap_bp(l.overlap);
                if ov > 0 {
                    link_overlaps.entry(e).or_insert(ov);
                    link_overlaps.entry(e.rev()).or_insert(ov);
                }
                edges.push(e);
            }
            GfaRecordRef::Path(p) => {
                let steps = parse_path_steps(p.segment_names);
                let hap_id = track_haps.then(|| {
                    let (sample, hap, _) = parse_path_sample_fields(p.name);
                    intern_hap_id(&mut hap_ids, &mut next_hap_id, sample, hap)
                });
                record_walk_support(
                    &steps,
                    hap_id,
                    &mut seg_ids,
                    &mut walk_support,
                    &mut hap_tracker,
                )?;
            }
            GfaRecordRef::Walk(w) => {
                let steps = parse_walk_steps(w.walk);
                let hap_id = track_haps
                    .then(|| intern_hap_id(&mut hap_ids, &mut next_hap_id, w.sample, w.haplotype));
                record_walk_support(
                    &steps,
                    hap_id,
                    &mut seg_ids,
                    &mut walk_support,
                    &mut hap_tracker,
                )?;
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    let n_segs = seg_ids.len();
    if seg_lens.len() < n_segs {
        seg_lens.resize(n_segs, 0);
        seg_seen_s.resize(n_segs, false);
    }

    Ok(CleanPass1 {
        seg_ids,
        seg_lens,
        seg_seen_s,
        flanks,
        edges,
        walk_support,
        link_overlaps,
        hap_tracker,
        horizon_bubbles,
    })
}

fn select_dropped_segments(pass1: &CleanPass1, cfg: &CleanConfig) -> (Vec<bool>, usize, usize) {
    let n_segs = pass1.seg_ids.len();
    let mut dropped = vec![false; n_segs];
    let mut dropped_low_hap_segs = 0usize;
    let mut clipped_bubbles = 0usize;
    if cfg.min_haplotypes == 0 && cfg.clip_bubbles_bp == 0 {
        return (dropped, 0, 0);
    }
    for idx in 0..n_segs {
        if !pass1.seg_seen_s[idx] {
            continue;
        }
        #[allow(clippy::cast_possible_truncation)]
        let id = idx as u32;
        let len = pass1.seg_lens[idx] as usize;
        let clip_bubble = cfg.clip_bubbles_bp > 0
            && pass1.horizon_bubbles.contains(&id)
            && len >= cfg.clip_bubbles_bp;
        if clip_bubble {
            dropped[idx] = true;
            clipped_bubbles += 1;
            continue;
        }
        if cfg.min_haplotypes > 0 && !pass1.hap_tracker.meets_min(id, cfg.min_haplotypes) {
            dropped[idx] = true;
            dropped_low_hap_segs += 1;
        }
    }
    (dropped, dropped_low_hap_segs, clipped_bubbles)
}

fn score_links_symmetrically(
    pass1: &CleanPass1,
    dropped_segs: &[bool],
    cfg: &CleanConfig,
) -> (HashMap<PackedEdge, f64>, usize, usize, usize) {
    use rayon::prelude::*;

    // Deduplicate surviving canonical edges so BOCPD junction scoring runs in parallel without
    // redundant work across reciprocal orientations.
    let mut dropped_links = 0usize;
    let mut unique_canons: Vec<PackedEdge> = Vec::new();
    let mut seen_canons: HashSet<PackedEdge> = HashSet::new();
    for &edge in &pass1.edges {
        if dropped_segs[edge.src_id() as usize] || dropped_segs[edge.dst_id() as usize] {
            dropped_links += 1;
            continue;
        }
        let canon = edge.canon();
        if seen_canons.insert(canon) {
            unique_canons.push(canon);
        }
    }
    drop(seen_canons);

    let decisions: Vec<(PackedEdge, f64)> = unique_canons
        .into_par_iter()
        .map_init(
            || (BocpdScratch::default(), Vec::new(), Vec::new()),
            |(scratch, concat, gc_prof), canon| {
                let sup = pass1.walk_support.get(&canon).copied().unwrap_or(0);
                let cp = if sup <= 1
                    && let (Some(u3), Some(v5)) = (
                        pass1.flanks.oriented_3p(canon.src_id(), canon.src_fwd()),
                        pass1.flanks.oriented_5p(canon.dst_id(), canon.dst_fwd()),
                    ) {
                    score_junction_cp_bits_with_scratch(&u3, &v5, scratch, concat, gc_prof)
                } else {
                    0.0
                };
                let sever =
                    (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);
                sever.then_some((canon, cp))
            },
        )
        .flatten()
        .collect();

    let mut severed_canons: HashMap<PackedEdge, f64> = HashMap::with_capacity(decisions.len());
    for (canon, cp) in decisions {
        severed_canons.insert(canon, cp);
    }

    let mut severed_info: HashMap<PackedEdge, f64> =
        HashMap::with_capacity(severed_canons.len().saturating_mul(2));
    let mut kept_links = 0usize;
    let mut severed_links = 0usize;
    for &edge in &pass1.edges {
        if dropped_segs[edge.src_id() as usize] || dropped_segs[edge.dst_id() as usize] {
            continue;
        }
        let canon = edge.canon();
        if let Some(&cp_bits) = severed_canons.get(&canon) {
            severed_links += 1;
            severed_info.insert(edge, cp_bits);
            severed_info.insert(edge.rev(), cp_bits);
        } else {
            kept_links += 1;
        }
    }
    (severed_info, kept_links, severed_links, dropped_links)
}

struct Pass2Ctx {
    seg_ids: SegInterner,
    seg_lens: Vec<u32>,
    link_overlaps: HashMap<PackedEdge, usize>,
    dropped_segs: Vec<bool>,
    severed_info: HashMap<PackedEdge, f64>,
}

fn emit_clean_path<W: Write>(p: GfaPathRef<'_>, ctx: &Pass2Ctx, out: &mut W) -> Result<bool> {
    let steps = parse_path_steps(p.segment_names);
    let is_dropped = |name: &str| {
        ctx.seg_ids
            .get(name)
            .is_some_and(|id| ctx.dropped_segs[id as usize])
    };
    let mut chunks: Vec<(usize, Vec<(&str, bool)>)> = Vec::new();
    let mut cur_start_idx = 0usize;
    let mut cur_chunk: Vec<(&str, bool)> = Vec::new();
    for (idx, &step) in steps.iter().enumerate() {
        if is_dropped(step.0) {
            if !cur_chunk.is_empty() {
                chunks.push((cur_start_idx, std::mem::take(&mut cur_chunk)));
            }
            continue;
        }
        if let Some(&prev) = cur_chunk.last() {
            let is_severed = ctx
                .seg_ids
                .get(prev.0)
                .zip(ctx.seg_ids.get(step.0))
                .is_some_and(|(u, v)| {
                    ctx.severed_info
                        .contains_key(&PackedEdge::new(u, prev.1, v, step.1))
                });
            if is_severed {
                chunks.push((cur_start_idx, std::mem::take(&mut cur_chunk)));
                cur_start_idx = idx;
            }
        } else {
            cur_start_idx = idx;
        }
        cur_chunk.push(step);
    }
    if !cur_chunk.is_empty() {
        chunks.push((cur_start_idx, cur_chunk));
    }

    if chunks.is_empty() {
        return Ok(!steps.is_empty());
    }
    if chunks.len() == 1 && chunks[0].1.len() == steps.len() {
        writeln!(out, "P\t{}\t{}\t{}", p.name, p.segment_names, p.overlaps)?;
        return Ok(false);
    }
    let orig_ovs: Vec<&str> = if p.overlaps == "*" || p.overlaps.is_empty() {
        Vec::new()
    } else {
        p.overlaps.split(',').collect()
    };
    for (ci, (chunk_start, chunk)) in chunks.iter().enumerate() {
        let segs: Vec<String> = chunk
            .iter()
            .map(|(s, fwd)| format!("{s}{}", if *fwd { '+' } else { '-' }))
            .collect();
        let ovs = if chunk.len() > 1 {
            if orig_ovs.len() >= chunk_start + chunk.len() - 1 {
                orig_ovs[*chunk_start..chunk_start + chunk.len() - 1].join(",")
            } else {
                vec!["0M"; chunk.len() - 1].join(",")
            }
        } else {
            "*".to_string()
        };
        writeln!(
            out,
            "P\t{}_part{}\t{}\t{}",
            p.name,
            ci + 1,
            segs.join(","),
            ovs
        )?;
    }
    Ok(true)
}

fn emit_clean_walk<W: Write>(w: GfaWalkRef<'_>, ctx: &Pass2Ctx, out: &mut W) -> Result<bool> {
    type WalkChunk<'a> = (usize, usize, Vec<(&'a str, bool)>);
    let emit_walk = |out: &mut dyn Write, cs: usize, ce: usize, walk_str: &str| -> Result<()> {
        if let Some(tags) = w.tags {
            writeln!(
                out,
                "W\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                w.sample, w.haplotype, w.seq_id, cs, ce, walk_str, tags
            )?;
        } else {
            writeln!(
                out,
                "W\t{}\t{}\t{}\t{}\t{}\t{}",
                w.sample, w.haplotype, w.seq_id, cs, ce, walk_str
            )?;
        }
        Ok(())
    };
    let steps = parse_walk_steps(w.walk);
    if steps.is_empty() {
        emit_walk(out, w.start, w.end, w.walk)?;
        return Ok(false);
    }

    let mut chunks: Vec<WalkChunk<'_>> = Vec::new();
    let mut cur_start = w.start;
    let mut cur_end = w.start;
    let mut cur_steps: Vec<(&str, bool)> = Vec::new();
    let mut prev_step_end = w.start;
    let mut prev_id: Option<u32> = None;

    for (idx, &step) in steps.iter().enumerate() {
        let step_id = ctx.seg_ids.get(step.0);
        let step_len = step_id.map_or(0, |id| ctx.seg_lens[id as usize] as usize);
        let edge_opt = if idx > 0 {
            prev_id
                .zip(step_id)
                .map(|(u, v)| PackedEdge::new(u, steps[idx - 1].1, v, step.1))
        } else {
            None
        };
        let ov_bp = edge_opt
            .and_then(|e| ctx.link_overlaps.get(&e).copied())
            .unwrap_or(0);
        let step_start = if idx == 0 {
            w.start
        } else {
            prev_step_end.saturating_sub(ov_bp)
        };
        let step_end = step_start.saturating_add(step_len);
        let dropped = step_id.is_some_and(|id| ctx.dropped_segs[id as usize]);

        if dropped {
            if !cur_steps.is_empty() {
                chunks.push((cur_start, cur_end, std::mem::take(&mut cur_steps)));
            }
            prev_step_end = step_end;
            prev_id = step_id;
            continue;
        }

        let is_severed =
            !cur_steps.is_empty() && edge_opt.is_some_and(|e| ctx.severed_info.contains_key(&e));
        if cur_steps.is_empty() {
            cur_start = step_start;
            cur_end = step_end;
            cur_steps.push(step);
        } else if is_severed {
            chunks.push((cur_start, cur_end, std::mem::take(&mut cur_steps)));
            cur_start = step_start;
            cur_end = step_end;
            cur_steps.push(step);
        } else {
            cur_steps.push(step);
            cur_end = step_end;
        }
        prev_step_end = step_end;
        prev_id = step_id;
    }
    if !cur_steps.is_empty() {
        chunks.push((cur_start, cur_end, cur_steps));
    }

    if chunks.is_empty() {
        return Ok(true);
    }
    if chunks.len() == 1 && chunks[0].2.len() == steps.len() {
        emit_walk(out, w.start, w.end, w.walk)?;
        return Ok(false);
    }
    for (cs, ce, chunk) in chunks {
        let mut walk_str = String::new();
        for (s, fwd) in chunk {
            walk_str.push(if fwd { '>' } else { '<' });
            walk_str.push_str(s);
        }
        emit_walk(out, cs, ce, &walk_str)?;
    }
    Ok(true)
}

pub(crate) fn clean_gfa<W: Write>(
    gfa_path: impl AsRef<std::path::Path>,
    cfg: &CleanConfig,
    out: &mut W,
) -> Result<CleanStats> {
    let gfa_path = gfa_path.as_ref();
    let pass1 = collect_clean_pass1(gfa_path, cfg)?;
    let total_links = pass1.edges.len();
    let (dropped_segs, dropped_low_hap_segs, clipped_bubbles) =
        select_dropped_segments(&pass1, cfg);
    let (severed_info, kept_links, severed_links, dropped_links) =
        score_links_symmetrically(&pass1, &dropped_segs, cfg);

    // Explicitly drop Pass 1 flank arena, walk support, edge list, and haplotype tracker before
    // streaming Pass 2 so peak memory remains minimal.
    let ctx = Pass2Ctx {
        seg_ids: pass1.seg_ids,
        seg_lens: pass1.seg_lens,
        link_overlaps: pass1.link_overlaps,
        dropped_segs,
        severed_info,
    };

    writeln!(out, "H\tVN:Z:1.1")?;
    let mut emitted_bubbles = 0usize;
    let mut split_walks = 0usize;

    stream_gfa_record_refs(gfa_path, |rec| {
        match rec {
            GfaRecordRef::Segment(s) => {
                if ctx
                    .seg_ids
                    .get(s.name)
                    .is_some_and(|id| ctx.dropped_segs[id as usize])
                {
                    return Ok(ControlFlow::Continue(()));
                }
                let seq_str = String::from_utf8_lossy(s.sequence);
                if let Some(tags) = s.tags {
                    writeln!(out, "S\t{}\t{}\t{}", s.name, seq_str, tags)?;
                } else {
                    writeln!(out, "S\t{}\t{}", s.name, seq_str)?;
                }
            }
            GfaRecordRef::Link(l) => {
                let u_opt = ctx.seg_ids.get(l.from);
                let v_opt = ctx.seg_ids.get(l.to);
                if u_opt.is_some_and(|u| ctx.dropped_segs[u as usize])
                    || v_opt.is_some_and(|v| ctx.dropped_segs[v as usize])
                {
                    return Ok(ControlFlow::Continue(()));
                }
                let edge_cp = u_opt.zip(v_opt).and_then(|(u, v)| {
                    let edge = PackedEdge::new(u, l.from_orient != b'-', v, l.to_orient != b'-');
                    ctx.severed_info.get(&edge).copied()
                });
                if let Some(cp_bits) = edge_cp {
                    emitted_bubbles += 1;
                    writeln!(
                        out,
                        "#HB\tfrom={}\tto={}\tid={}__{}__hb_{}\tHB:i:1\tCP:f:{:.2}",
                        l.from, l.to, l.from, l.to, emitted_bubbles, cp_bits
                    )?;
                } else {
                    writeln!(
                        out,
                        "L\t{}\t{}\t{}\t{}\t{}",
                        l.from, l.from_orient as char, l.to, l.to_orient as char, l.overlap
                    )?;
                }
            }
            GfaRecordRef::Path(p) => {
                if emit_clean_path(p, &ctx, out)? {
                    split_walks += 1;
                }
            }
            GfaRecordRef::Walk(w) => {
                if emit_clean_walk(w, &ctx, out)? {
                    split_walks += 1;
                }
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    Ok(CleanStats {
        total_links,
        kept_links,
        severed_links,
        dropped_links,
        emitted_bubbles,
        split_walks,
        dropped_low_hap_segs,
        clipped_bubbles,
    })
}
