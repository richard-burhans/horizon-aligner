//! Multi-sample VCFv4.2 export from Tarjan biedged snarls (`LV=0`) and optional
//! nested reference-anchored local bubbles (`LV=1`) projected onto a reference walk.

use anyhow::{Result, bail};
use horizon_aligner::dna::complement_base;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as _;
use std::io::Write;

use crate::resolve::SegSeqArena;
use crate::snarl::{PackedEdge, PackedGraph, PackedStep};
use crate::vcf_nested::{genotype_chrom_nested_bubbles, plan_nested_bubbles_by_chrom};

/// Summary statistics for VCF records emitted from Tarjan snarls and optional nested bubbles.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub(crate) struct VcfStats {
    pub(crate) sites_emitted: usize,
    pub(crate) nested_sites: usize,
    pub(crate) snps: usize,
    pub(crate) mnvs: usize,
    pub(crate) indels: usize,
    pub(crate) multiallelic: usize,
    pub(crate) samples: usize,
}

/// A resolved GFA walk/path captured for snarl-to-VCF projection.
#[derive(Debug, Clone)]
pub(crate) struct CapturedWalk {
    pub(crate) sample: String,
    pub(crate) haplotype: String,
    pub(crate) seq_id: String,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) steps: Vec<PackedStep>,
}

pub(super) type SnarlKey = (PackedStep, PackedStep);
pub(super) type HapId = u32;
type Lv0Job = (usize, usize, usize, SnarlKey, bool);

#[inline]
pub(super) fn canonical_snarl_key(u: PackedStep, v: PackedStep) -> (SnarlKey, bool) {
    let fwd = (u, v);
    let rev = (v.rev(), u.rev());
    if fwd <= rev {
        (fwd, false)
    } else {
        (rev, true)
    }
}

/// Normalize a DNA byte to the strict VCFv4.2 allele alphabet (`A`, `C`, `G`, `T`, `N`).
///
/// Draft genome assemblies occasionally contain IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`,
/// `M`, `B`, `D`, `H`, `V`, etc.), which VCFv4.2 forbids in `REF`/`ALT` alleles and `bcftools norm`
/// rejects (`Non-ACGTN reference/alternate allele`). Mapping any non-`ACGT` byte to `N` matches
/// `htslib`/`bcftools` reference normalization and ensures distinct IUPAC traversals that collapse
/// to the same sanitized sequence are deduplicated cleanly before `AC`/`AN`/`AT` emission.
#[inline]
const fn sanitize_vcf_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' => b'T',
        _ => b'N',
    }
}

/// O(1) zero-allocation lookup for the last base of oriented segment `step`.
#[inline]
fn last_oriented_base(seg_seqs: &SegSeqArena, step: PackedStep) -> u8 {
    let raw = seg_seqs.get(step.id()).unwrap_or_default();
    let b = if step.fwd() {
        raw.last().copied()
    } else {
        raw.first().copied().map(complement_base)
    };
    b.map_or(b'N', sanitize_vcf_base)
}

fn render_traversal(
    graph: &PackedGraph,
    u: PackedStep,
    interior: &[PackedStep],
    v: PackedStep,
) -> String {
    let mut s = String::with_capacity((interior.len() + 2).saturating_mul(6));
    for &step in std::iter::once(&u)
        .chain(interior.iter())
        .chain(std::iter::once(&v))
    {
        s.push(if step.fwd() { '>' } else { '<' });
        graph.push_name(step.id(), &mut s);
    }
    s
}

/// Append oriented interior DNA bases (trimming non-0M overlaps) into a reusable buffer `out`.
/// Returns `false` if any interior segment lacks an explicit DNA sequence (`*`).
fn extract_interior_seq_into(
    graph: &PackedGraph,
    seg_seqs: &SegSeqArena,
    u: PackedStep,
    interior: &[PackedStep],
    v: PackedStep,
    out: &mut Vec<u8>,
) -> bool {
    out.clear();
    let mut prev = u;
    for &step in interior {
        let Some(raw) = seg_seqs.get(step.id()) else {
            return false;
        };
        let ov = graph
            .link_overlaps
            .get(&PackedEdge::directed(prev, step))
            .copied()
            .unwrap_or(0) as usize;
        if ov < raw.len() {
            if step.fwd() {
                out.extend(raw[ov..].iter().copied().map(sanitize_vcf_base));
            } else {
                // Reverse-complement slice after skipping the first `ov` revcomp bases
                // (which correspond to the last `ov` forward bases of `raw`).
                let keep = raw.len() - ov;
                out.extend(
                    raw[..keep]
                        .iter()
                        .rev()
                        .map(|&b| sanitize_vcf_base(complement_base(b))),
                );
            }
        }
        prev = step;
    }
    // Trim closing link overlap into anchor `v` on non-0M GFA v1.0 graphs so alleles do not
    // duplicate suffix bases already belonging to closing anchor `v`.
    if !interior.is_empty() {
        let close_ov = graph
            .link_overlaps
            .get(&PackedEdge::directed(prev, v))
            .copied()
            .unwrap_or(0) as usize;
        if close_ov > 0 {
            let new_len = out.len().saturating_sub(close_ov);
            out.truncate(new_len);
        }
    }
    true
}

/// Compact haplotype traversal state inside a snarl or nested bubble.
///
/// Zero heap allocation for deletions (`Empty`) and SNPs/single-node alleles (`One`).
pub(super) enum TraversalEntry {
    Empty,
    One(PackedStep),
    Many(Box<[PackedStep]>),
    /// Multiple distinct visits (e.g., tandem repeat or split walk): emit missing GT `.`.
    Conflict,
}

impl TraversalEntry {
    #[inline]
    fn from_slice(raw: &[PackedStep], rev: bool) -> Self {
        match raw {
            [] => Self::Empty,
            [single] => Self::One(if rev { single.rev() } else { *single }),
            many => {
                let boxed: Box<[PackedStep]> = if rev {
                    many.iter().rev().map(|&s| s.rev()).collect()
                } else {
                    many.into()
                };
                Self::Many(boxed)
            }
        }
    }

    /// Populate `buf` with oriented interior steps (or return a direct slice when zero-copy).
    #[inline]
    fn orient_into<'a>(
        &'a self,
        rev: bool,
        buf: &'a mut Vec<PackedStep>,
    ) -> Option<&'a [PackedStep]> {
        match self {
            Self::Conflict => None,
            Self::Empty => Some(&[]),
            Self::One(s) => {
                if rev {
                    buf.clear();
                    buf.push(s.rev());
                    Some(buf.as_slice())
                } else {
                    Some(std::slice::from_ref(s))
                }
            }
            Self::Many(slice) => {
                if rev {
                    buf.clear();
                    buf.extend(slice.iter().rev().map(|&s| s.rev()));
                    Some(buf.as_slice())
                } else {
                    Some(slice.as_ref())
                }
            }
        }
    }
}

#[inline]
pub(super) fn insert_traversal(
    map: &mut HashMap<HapId, TraversalEntry>,
    hap_id: HapId,
    raw: &[PackedStep],
    rev: bool,
) {
    match map.entry(hap_id) {
        std::collections::hash_map::Entry::Vacant(vac) => {
            vac.insert(TraversalEntry::from_slice(raw, rev));
        }
        std::collections::hash_map::Entry::Occupied(mut occ) => {
            occ.insert(TraversalEntry::Conflict);
        }
    }
}

pub(super) struct VcfBuildCtx<'a> {
    pub(super) graph: &'a PackedGraph,
    pub(super) seg_seqs: &'a SegSeqArena,
    pub(super) samples: &'a [(String, Vec<HapId>)],
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum VariantClass {
    Snp,
    Mnv,
    Indel,
}

/// Single-allocation packed VCF record (`<REF><ALT_CSV><INFO\tGT\tGTS>`), replacing
/// ~80 per-site heap `String`/`Vec` allocations with one exact-fit `Box<str>`.
pub(super) struct VcfRecord {
    pos: usize,
    /// 1-based position of the snarl's reference interior before any 1-bp left-anchor shift.
    interior_pos_1: usize,
    level: u8,
    variant_class: VariantClass,
    is_multiallelic: bool,
    ref_len: u32,
    alt_len: u32,
    payload: Box<str>,
}

impl VcfRecord {
    #[inline]
    fn ref_allele(&self) -> &str {
        &self.payload[..self.ref_len as usize]
    }

    #[inline]
    fn alt_csv(&self) -> &str {
        let s = self.ref_len as usize;
        &self.payload[s..s + self.alt_len as usize]
    }

    #[inline]
    fn info_and_gts(&self) -> &str {
        &self.payload[(self.ref_len as usize) + (self.alt_len as usize)..]
    }
}

/// Minimum reference allele length (bp) to annotate an inverted snarl allele with `SVTYPE=INV`.
const MIN_SV_INV_BP: usize = 50;

/// Return `true` when an ALT allele represents a structural inversion (`>= MIN_SV_INV_BP` bp) of
/// `ref_seq`: either a deduplicated reverse-oriented graph traversal covering `>= 50%` of
/// `max(|REF|, |ALT|)` or an indel-tolerant 9-mer reverse-complement match (`>= 70%` RC k-mer
/// containment, `< 35%` forward containment).
fn is_inversion_allele(
    seg_seqs: &SegSeqArena,
    ref_seq: &[u8],
    alt_seq: &[u8],
    ref_interior: &[PackedStep],
    alt_interior: &[PackedStep],
) -> bool {
    const K: usize = 9;
    let rlen = ref_seq.len();
    let alen = alt_seq.len();
    if rlen < MIN_SV_INV_BP || alen < MIN_SV_INV_BP {
        return false;
    }
    let max_len = rlen.max(alen);
    if !ref_interior.is_empty() && !alt_interior.is_empty() {
        let mut ref_counts: HashMap<PackedStep, usize> = HashMap::new();
        for &r in ref_interior {
            *ref_counts.entry(r).or_insert(0) += 1;
        }
        let mut inverted_shared_bp = 0usize;
        for &a in alt_interior {
            if let Some(cnt) = ref_counts.get_mut(&a.rev())
                && *cnt > 0
            {
                *cnt -= 1;
                inverted_shared_bp = inverted_shared_bp.saturating_add(seg_seqs.seq_len(a.id()));
            }
        }
        if inverted_shared_bp >= MIN_SV_INV_BP && inverted_shared_bp.saturating_mul(2) >= max_len {
            return true;
        }
    }
    if rlen.abs_diff(alen) > (rlen / 5).max(50) || rlen < K || alen < K {
        return false;
    }
    let enc = |b: u8| -> Option<u32> {
        match b.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    };
    let mask = (1u32 << (2 * K)) - 1;
    let mut ref_kmers: HashSet<u32> = HashSet::with_capacity(rlen);
    let (mut code, mut valid) = (0u32, 0usize);
    for &b in ref_seq {
        if let Some(v) = enc(b) {
            code = ((code << 2) | v) & mask;
            valid += 1;
            if valid >= K {
                ref_kmers.insert(code);
            }
        } else {
            code = 0;
            valid = 0;
        }
    }
    if ref_kmers.len() < 16 {
        return false;
    }
    let (mut fwd_code, mut rc_code, mut a_valid) = (0u32, 0u32, 0usize);
    let (mut total, mut fwd_hits, mut rc_hits) = (0usize, 0usize, 0usize);
    let shift = 2 * (K - 1);
    for &b in alt_seq {
        if let Some(v) = enc(b) {
            fwd_code = ((fwd_code << 2) | v) & mask;
            rc_code = (rc_code >> 2) | ((3 - v) << shift);
            a_valid += 1;
            if a_valid >= K {
                total += 1;
                if ref_kmers.contains(&fwd_code) {
                    fwd_hits += 1;
                }
                if ref_kmers.contains(&rc_code) {
                    rc_hits += 1;
                }
            }
        } else {
            fwd_code = 0;
            rc_code = 0;
            a_valid = 0;
        }
    }
    total >= 16
        && rc_hits.saturating_mul(100) >= total.saturating_mul(70)
        && fwd_hits.saturating_mul(100) < total.saturating_mul(35)
}

/// Build a single polymorphic VCF record for anchor pair `(u, v)` on the current reference contig.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn build_vcf_record(
    ctx: &VcfBuildCtx<'_>,
    u: PackedStep,
    v: PackedStep,
    ref_interior: &[PackedStep],
    ref_interior_start_0: usize,
    ref_is_rev: bool,
    hap_map: &HashMap<HapId, TraversalEntry>,
    level: u8,
    scratch_int: &mut Vec<PackedStep>,
    scratch_seq: &mut Vec<u8>,
) -> Option<VcfRecord> {
    if !ctx.seg_seqs.has_seq(u.id()) || !ctx.seg_seqs.has_seq(v.id()) {
        return None;
    }
    if !extract_interior_seq_into(ctx.graph, ctx.seg_seqs, u, ref_interior, v, scratch_seq) {
        return None;
    }
    let ref_seq = scratch_seq.clone();
    let mut allele_seqs: Vec<Vec<u8>> = vec![ref_seq];
    let mut allele_travs: Vec<String> = Vec::new();
    let mut allele_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();
    allele_to_idx.insert(allele_seqs[0].clone(), 0);
    let mut has_inversion = false;

    let mut per_sample_calls: Vec<Vec<Option<usize>>> = Vec::with_capacity(ctx.samples.len());
    for (_, haps) in ctx.samples {
        let mut hap_calls = Vec::with_capacity(haps.len());
        for &hap_id in haps {
            let Some(entry) = hap_map.get(&hap_id) else {
                hap_calls.push(None);
                continue;
            };
            let Some(oriented_int) = entry.orient_into(ref_is_rev, scratch_int) else {
                hap_calls.push(None);
                continue;
            };
            if extract_interior_seq_into(ctx.graph, ctx.seg_seqs, u, oriented_int, v, scratch_seq) {
                let aidx = if let Some(&idx) = allele_to_idx.get(scratch_seq.as_slice()) {
                    idx
                } else {
                    if allele_travs.is_empty() {
                        // Lazily render reference traversal only once a site proves polymorphic.
                        allele_travs.push(render_traversal(ctx.graph, u, ref_interior, v));
                    }
                    if !has_inversion
                        && is_inversion_allele(
                            ctx.seg_seqs,
                            &allele_seqs[0],
                            scratch_seq,
                            ref_interior,
                            oriented_int,
                        )
                    {
                        has_inversion = true;
                    }
                    let idx = allele_seqs.len();
                    let trav = render_traversal(ctx.graph, u, oriented_int, v);
                    let owned = scratch_seq.clone();
                    allele_to_idx.insert(owned.clone(), idx);
                    allele_seqs.push(owned);
                    allele_travs.push(trav);
                    idx
                };
                hap_calls.push(Some(aidx));
            } else {
                hap_calls.push(None);
            }
        }
        per_sample_calls.push(hap_calls);
    }

    if allele_seqs.len() <= 1 {
        return None;
    }

    // Apply VCFv4.2 1-bp left-anchor base rule whenever REF or any ALT is empty.
    let needs_anchor = allele_seqs.iter().any(Vec::is_empty);
    let vcf_pos = if needs_anchor {
        let anchor_base = last_oriented_base(ctx.seg_seqs, u);
        for a in &mut allele_seqs {
            a.insert(0, anchor_base);
        }
        ref_interior_start_0.max(1)
    } else {
        ref_interior_start_0.saturating_add(1)
    };

    let mut ac = vec![0usize; allele_seqs.len() - 1];
    let mut an = 0usize;
    let mut ns = 0usize;
    for calls in &per_sample_calls {
        let mut any_called = false;
        for call in calls {
            if let Some(aidx) = *call {
                any_called = true;
                an += 1;
                if aidx > 0 {
                    ac[aidx - 1] += 1;
                }
            }
        }
        if any_called {
            ns += 1;
        }
    }

    let ref_b_len = allele_seqs[0].len();
    let variant_class = if ref_b_len == 1 && allele_seqs[1..].iter().all(|a| a.len() == 1) {
        VariantClass::Snp
    } else if ref_b_len > 1 && allele_seqs[1..].iter().all(|a| a.len() == ref_b_len) {
        VariantClass::Mnv
    } else {
        VariantClass::Indel
    };
    let is_multiallelic = allele_seqs.len() > 2;

    let mut payload = String::with_capacity(
        allele_seqs[0].len()
            + allele_seqs[1..].iter().map(Vec::len).sum::<usize>()
            + allele_travs.iter().map(String::len).sum::<usize>()
            + ctx.samples.len().saturating_mul(3)
            + 48,
    );
    payload.push_str(&String::from_utf8_lossy(&allele_seqs[0]));
    let ref_len = u32::try_from(payload.len()).ok()?;
    for (idx, alt) in allele_seqs[1..].iter().enumerate() {
        if idx > 0 {
            payload.push(',');
        }
        payload.push_str(&String::from_utf8_lossy(alt));
    }
    let alt_len = u32::try_from(payload.len() - ref_len as usize).ok()?;

    payload.push_str("AT=");
    for (idx, trav) in allele_travs.iter().enumerate() {
        if idx > 0 {
            payload.push(',');
        }
        payload.push_str(trav);
    }
    payload.push_str(";AC=");
    for (idx, count) in ac.iter().enumerate() {
        if idx > 0 {
            payload.push(',');
        }
        let _ = write!(payload, "{count}");
    }
    if has_inversion {
        let svlen = allele_seqs[0].len();
        let end_pos = vcf_pos.saturating_add(svlen.saturating_sub(1));
        let _ = write!(
            payload,
            ";AN={an};NS={ns};SVTYPE=INV;END={end_pos};SVLEN={svlen};LV={level}\tGT"
        );
    } else {
        let _ = write!(payload, ";AN={an};NS={ns};LV={level}\tGT");
    }
    for calls in &per_sample_calls {
        payload.push('\t');
        for (h_idx, call) in calls.iter().enumerate() {
            if h_idx > 0 {
                payload.push('|');
            }
            if let Some(aidx) = *call {
                let _ = write!(payload, "{aidx}");
            } else {
                payload.push('.');
            }
        }
    }

    Some(VcfRecord {
        pos: vcf_pos,
        interior_pos_1: ref_interior_start_0.saturating_add(1),
        level,
        variant_class,
        is_multiallelic,
        ref_len,
        alt_len,
        payload: payload.into_boxed_str(),
    })
}

fn write_vcf_header<W: Write>(
    out: &mut W,
    contigs: &BTreeMap<String, usize>,
    samples: &[(String, Vec<HapId>)],
) -> Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##source=horizon_resolve")?;
    for (contig, len) in contigs {
        writeln!(out, "##contig=<ID={contig},length={len}>")?;
    }
    writeln!(
        out,
        "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele traversal as path in graph\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in snarl tree (0=top-level, 1=nested)\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant (e.g., INV)\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant on reference\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for (sample, _) in samples {
        write!(out, "\t{sample}")?;
    }
    writeln!(out)?;
    Ok(())
}

/// Project Tarjan snarls (`LV=0`) and optional nested reference-anchored bubbles (`LV=1`) onto
/// `ref_sample` (or the first observed sample when `None`) and stream a multi-sample VCFv4.2
/// one reference chromosome at a time with Rayon-parallel allele/genotype construction.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn export_snarl_vcf<W: Write>(
    graph: &PackedGraph,
    seg_seqs: &SegSeqArena,
    walks: &[CapturedWalk],
    ref_sample: Option<&str>,
    nested_snarls: bool,
    max_nested_steps: usize,
    out: &mut W,
) -> Result<VcfStats> {
    if walks.is_empty() {
        if let Some(r) = ref_sample {
            bail!("--vcf-ref sample `{r}` not found in GFA walks/paths");
        }
        writeln!(out, "##fileformat=VCFv4.2")?;
        writeln!(out, "##source=horizon_resolve")?;
        writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
        return Ok(VcfStats::default());
    }

    let ref_name: &str = if let Some(req) = ref_sample {
        if !walks.iter().any(|w| w.sample == req) {
            bail!("--vcf-ref sample `{req}` not found in GFA walks/paths");
        }
        req
    } else {
        walks[0].sample.as_str()
    };

    // Collect deterministic sorted samples and intern every `(sample, haplotype)` to a `u32` HapId.
    let mut sample_haps_map: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for w in walks {
        sample_haps_map
            .entry(w.sample.clone())
            .or_default()
            .insert(w.haplotype.clone());
    }
    let mut hap_to_id: HashMap<(&str, &str), HapId> = HashMap::new();
    let mut samples_raw: Vec<(String, Vec<String>)> = Vec::with_capacity(sample_haps_map.len());
    for (s, hs) in sample_haps_map {
        samples_raw.push((s, hs.into_iter().collect()));
    }
    let mut samples: Vec<(String, Vec<HapId>)> = Vec::with_capacity(samples_raw.len());
    let mut next_hap_id: HapId = 0;
    for (s, hs) in &samples_raw {
        let mut ids = Vec::with_capacity(hs.len());
        for h in hs {
            let id = next_hap_id;
            next_hap_id = next_hap_id.saturating_add(1);
            hap_to_id.insert((s.as_str(), h.as_str()), id);
            ids.push(id);
        }
        samples.push((s.clone(), ids));
    }
    let walk_hap_ids: Vec<HapId> = walks
        .iter()
        .map(|w| hap_to_id[&(w.sample.as_str(), w.haplotype.as_str())])
        .collect();

    let primary_ref_hap: String = samples_raw
        .iter()
        .find(|(s, _)| s == ref_name)
        .and_then(|(_, hs)| hs.first().cloned())
        .unwrap_or_else(|| "0".to_owned());

    let ctx = VcfBuildCtx {
        graph,
        seg_seqs,
        samples: &samples,
    };

    // Restrict top-level VCF snarl boundary anchors to segments actually visited by the primary
    // reference walk AND strictly Tarjan bridges (`is_bridge_seg`), falling back to `is_anchor`
    // only on reference walks that contain `< 2` bridges (e.g. cyclic contigs). This prevents
    // `>= 200 bp` shared non-reference insertions and `>= 200 bp` bypassed reference deletions from
    // splitting snarl keys and being silently dropped.
    let mut vcf_is_anchor = vec![false; graph.num_segs()];
    for w in walks {
        if w.sample != ref_name || w.haplotype != primary_ref_hap {
            continue;
        }
        let bridge_cnt = w
            .steps
            .iter()
            .filter(|step| {
                let idx = step.id() as usize;
                graph.is_anchor.get(idx).copied().unwrap_or(false)
                    && graph.is_bridge_seg.get(idx).copied().unwrap_or(false)
            })
            .count();
        let use_strict_bridge = bridge_cnt >= 2;
        for &step in &w.steps {
            let idx = step.id() as usize;
            let ok = if use_strict_bridge {
                graph.is_anchor.get(idx).copied().unwrap_or(false)
                    && graph.is_bridge_seg.get(idx).copied().unwrap_or(false)
            } else {
                graph.is_anchor.get(idx).copied().unwrap_or(false)
            };
            if ok && idx < vcf_is_anchor.len() {
                vcf_is_anchor[idx] = true;
            }
        }
    }

    // Index every walk's traversal across consecutive Tarjan anchor pairs under canonical snarl key.
    let mut snarl_traversals: HashMap<SnarlKey, HashMap<HapId, TraversalEntry>> = HashMap::new();
    for (w, &hap_id) in walks.iter().zip(&walk_hap_ids) {
        let anchors: Vec<usize> = w
            .steps
            .iter()
            .enumerate()
            .filter_map(|(idx, &step)| {
                vcf_is_anchor
                    .get(step.id() as usize)
                    .copied()
                    .unwrap_or(false)
                    .then_some(idx)
            })
            .collect();
        for pair in anchors.windows(2) {
            let u = w.steps[pair[0]];
            let v = w.steps[pair[1]];
            let raw_interior = &w.steps[(pair[0] + 1)..pair[1]];
            let (canon_key, is_rev) = canonical_snarl_key(u, v);
            insert_traversal(
                snarl_traversals.entry(canon_key).or_default(),
                hap_id,
                raw_interior,
                is_rev,
            );
        }
    }

    // Collect primary reference walks, sorted contigs, and LV=0 snarl jobs.
    let mut contigs: BTreeMap<String, usize> = BTreeMap::new();
    let mut ref_walk_indices: Vec<usize> = Vec::new();
    let mut ref_step_ends: Vec<Vec<usize>> = Vec::new();
    for (w_idx, w) in walks.iter().enumerate() {
        if w.sample != ref_name || w.haplotype != primary_ref_hap {
            continue;
        }
        let mut step_ends: Vec<usize> = Vec::with_capacity(w.steps.len());
        let mut cursor = w.start;
        for (idx, &step) in w.steps.iter().enumerate() {
            let slen = seg_seqs.seq_len(step.id()).max(graph.seg_len(step.id()));
            let ov = if idx > 0 {
                graph
                    .link_overlaps
                    .get(&PackedEdge::directed(w.steps[idx - 1], step))
                    .copied()
                    .unwrap_or(0) as usize
            } else {
                0
            };
            cursor = cursor.saturating_add(slen.saturating_sub(ov));
            step_ends.push(cursor);
        }
        let contig_len = w.end.max(cursor);
        let entry = contigs.entry(w.seq_id.clone()).or_insert(0);
        *entry = (*entry).max(contig_len);
        ref_walk_indices.push(w_idx);
        ref_step_ends.push(step_ends);
    }

    write_vcf_header(out, &contigs, &samples)?;

    // Map each reference walk to its alphabetically sorted contig index.
    let chrom_names: Vec<String> = contigs.keys().cloned().collect();
    let chrom_to_idx: HashMap<&str, usize> = chrom_names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    let ref_walk_chrom_ids: Vec<usize> = ref_walk_indices
        .iter()
        .map(|&w_idx| chrom_to_idx[walks[w_idx].seq_id.as_str()])
        .collect();

    let mut lv0_keys: HashSet<SnarlKey> = HashSet::new();
    let mut lv0_jobs_by_chrom: Vec<Vec<Lv0Job>> = vec![Vec::new(); chrom_names.len()];
    for (rw_i, &w_idx) in ref_walk_indices.iter().enumerate() {
        let w = &walks[w_idx];
        let c_idx = ref_walk_chrom_ids[rw_i];
        let anchors: Vec<usize> = w
            .steps
            .iter()
            .enumerate()
            .filter_map(|(idx, &step)| {
                vcf_is_anchor
                    .get(step.id() as usize)
                    .copied()
                    .unwrap_or(false)
                    .then_some(idx)
            })
            .collect();
        for pair in anchors.windows(2) {
            let (canon_key, ref_is_rev) = canonical_snarl_key(w.steps[pair[0]], w.steps[pair[1]]);
            lv0_keys.insert(canon_key);
            if snarl_traversals.contains_key(&canon_key) {
                lv0_jobs_by_chrom[c_idx].push((rw_i, pair[0], pair[1], canon_key, ref_is_rev));
            }
        }
    }

    let mut lv0_records_by_chrom: Vec<Vec<VcfRecord>> = lv0_jobs_by_chrom
        .into_iter()
        .map(|jobs| {
            jobs.par_iter()
                .map_init(
                    || (Vec::new(), Vec::new()),
                    |(scratch_int, scratch_seq), &(rw_i, a_start, a_end, canon_key, ref_is_rev)| {
                        let hap_map = snarl_traversals.get(&canon_key)?;
                        let rw = &walks[ref_walk_indices[rw_i]];
                        build_vcf_record(
                            &ctx,
                            rw.steps[a_start],
                            rw.steps[a_end],
                            &rw.steps[(a_start + 1)..a_end],
                            ref_step_ends[rw_i][a_start],
                            ref_is_rev,
                            hap_map,
                            0,
                            scratch_int,
                            scratch_seq,
                        )
                    },
                )
                .flatten()
                .collect()
        })
        .collect();

    // Free top-level snarl map before discovering and streaming per-chromosome LV=1 bubbles.
    drop(snarl_traversals);

    let max_steps = max_nested_steps.max(2);
    let mut nested_plan = nested_snarls.then(|| {
        plan_nested_bubbles_by_chrom(
            &ctx,
            walks,
            &ref_walk_indices,
            &ref_walk_chrom_ids,
            chrom_names.len(),
            &lv0_keys,
            max_steps,
        )
    });
    drop(lv0_keys);

    let mut stats = VcfStats {
        samples: samples.len(),
        ..VcfStats::default()
    };
    let mut snarl_idx = 0usize;

    for (c_idx, chrom_name) in chrom_names.iter().enumerate() {
        let mut chrom_records = std::mem::take(&mut lv0_records_by_chrom[c_idx]);
        if let Some(plan) = nested_plan.as_mut() {
            let bubbles = std::mem::take(&mut plan.bubbles_by_chrom[c_idx]);
            let touching = std::mem::take(&mut plan.walks_by_chrom[c_idx]);
            if !bubbles.is_empty() {
                let lv1 = genotype_chrom_nested_bubbles(
                    &ctx,
                    walks,
                    &walk_hap_ids,
                    &ref_walk_indices,
                    &ref_walk_chrom_ids,
                    &ref_step_ends,
                    &plan.seg_ref_loc,
                    c_idx,
                    &bubbles,
                    &touching,
                    max_steps,
                );
                chrom_records.extend(lv1);
            }
        }
        chrom_records.par_sort_unstable_by(|a, b| {
            a.pos
                .cmp(&b.pos)
                .then_with(|| a.level.cmp(&b.level))
                .then_with(|| b.ref_len.cmp(&a.ref_len))
                .then_with(|| a.ref_allele().cmp(b.ref_allele()))
                .then_with(|| a.alt_csv().split(',').cmp(b.alt_csv().split(',')))
        });

        // Pre-index every native LV=0 record's span `(pos, end_pos, assigned_snarl_id)` so a
        // left-anchored LV=1 indel/insertion that shifted 1 bp left of its enclosing LV=0 parent
        // still links via `PS=` instead of promoting into an overlapping LV=0 record.
        let base_snarl_idx = snarl_idx;
        let native_lv0: Vec<(usize, usize, usize)> = chrom_records
            .iter()
            .enumerate()
            .filter(|(_, r)| r.level == 0)
            .map(|(i, r)| {
                (
                    r.pos,
                    r.pos.saturating_add(r.ref_len as usize),
                    base_snarl_idx + i + 1,
                )
            })
            .collect();
        let mut next_lv0_i = 0usize;
        let mut active_lv0: Option<(usize, usize)> = None;

        for rec in &chrom_records {
            snarl_idx += 1;
            stats.sites_emitted += 1;
            match rec.variant_class {
                VariantClass::Snp => stats.snps += 1,
                VariantClass::Mnv => stats.mnvs += 1,
                VariantClass::Indel => stats.indels += 1,
            }
            if rec.is_multiallelic {
                stats.multiallelic += 1;
            }
            let rec_end = rec.pos.saturating_add(rec.ref_len as usize);
            let parent_id = if rec.level == 0 {
                active_lv0 = Some((snarl_idx, rec_end));
                None
            } else {
                while next_lv0_i < native_lv0.len() && native_lv0[next_lv0_i].1 <= rec.pos {
                    next_lv0_i += 1;
                }
                let active_match =
                    active_lv0.and_then(|(pid, end_pos)| (rec.pos < end_pos).then_some(pid));
                active_match.or_else(|| {
                    native_lv0
                        .get(next_lv0_i)
                        .and_then(|&(lv0_pos, lv0_end, lv0_pid)| {
                            let overlaps_span = lv0_pos < rec_end;
                            let interior_inside =
                                rec.interior_pos_1 >= lv0_pos && rec.interior_pos_1 < lv0_end;
                            (overlaps_span || interior_inside).then_some(lv0_pid)
                        })
                })
            };
            if let Some(pid) = parent_id
                && let Some((info, gt_part)) = rec.info_and_gts().split_once('\t')
            {
                stats.nested_sites += 1;
                writeln!(
                    out,
                    "{chrom_name}\t{}\tsnarl_{snarl_idx}\t{}\t{}\t.\tPASS\t{info};PS=snarl_{pid}\t{gt_part}",
                    rec.pos,
                    rec.ref_allele(),
                    rec.alt_csv(),
                )?;
            } else if rec.level > 0 {
                // Promote orphaned nested bubble (whose enclosing LV=0 snarl was uncallable/unemitted)
                // to LV=0 and update `active_lv0` so overlapping inner orphan bubbles attach via PS=
                // instead of emitting conflicting overlapping LV=0 records.
                active_lv0 = Some((snarl_idx, rec_end));
                let promoted = rec.info_and_gts().replacen(";LV=1\t", ";LV=0\t", 1);
                writeln!(
                    out,
                    "{chrom_name}\t{}\tsnarl_{snarl_idx}\t{}\t{}\t.\tPASS\t{promoted}",
                    rec.pos,
                    rec.ref_allele(),
                    rec.alt_csv(),
                )?;
            } else {
                writeln!(
                    out,
                    "{chrom_name}\t{}\tsnarl_{snarl_idx}\t{}\t{}\t.\tPASS\t{}",
                    rec.pos,
                    rec.ref_allele(),
                    rec.alt_csv(),
                    rec.info_and_gts(),
                )?;
            }
        }
    }

    Ok(stats)
}
