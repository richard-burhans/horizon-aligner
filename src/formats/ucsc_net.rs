//! Hierarchical UCSC `.net` builder, serializer, and combined chain/net file writer.

use super::ucsc::{UcscChain, UcscChainBlock, build_ucsc_chains, write_ucsc_chains};
use crate::coord::IntervalSet;
use crate::mdl::AlignedRun;
use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::io::Write;

/// Minimum empty target gap (bp) emitted as a leaf `gap` line in `.net` (matches `chainNet -minSpace=25`).
const MIN_NET_GAP_BP: usize = 25;
/// Maximum nesting depth of `fill` levels inside a `.net` hierarchy.
const MAX_NET_DEPTH: usize = 15;
/// Maximum query distance (bp) from parent gap to classify a nested fill as `syn` or `inv`.
const MAX_SYNTENIC_QFAR_BP: usize = 200_000;

/// One `fill` node in a UCSC `.net` hierarchy (coordinates are always `+` strand).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UcscNetFill {
    /// 0-based forward target start.
    pub t_start: usize,
    /// Target span length (bp).
    pub t_size: usize,
    /// Query contig name.
    pub q_name: String,
    /// Query alignment orientation (`'+'` or `'-'`).
    pub q_strand: char,
    /// 0-based **forward** (`+`) query start.
    pub q_start: usize,
    /// Query span length (bp).
    pub q_size: usize,
    /// Originating [`UcscChain::id`].
    pub chain_id: u64,
    /// Score of the chain slice used in this fill.
    pub score: i64,
    /// Aligned ungapped bases inside this fill.
    pub ali: usize,
    /// Overlap (bp) between this fill's query span and its parent gap's query span.
    pub q_over: usize,
    /// Distance (bp) from the parent gap's query span (`-1` when query contigs differ).
    pub q_far: i64,
    /// Query bases in this fill already covered by higher-scoring fills.
    pub q_dup: usize,
    /// UCSC `netSyntenic` class (`"top"`, `"syn"`, `"inv"`, `"dup"`, or `"nonSyn"`).
    pub fill_type: &'static str,
    /// Nested gaps inside this fill (in ascending `t_start` order).
    pub gaps: Vec<UcscNetGap>,
}

/// One `gap` node inside a [`UcscNetFill`] (coordinates are always `+` strand).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UcscNetGap {
    /// 0-based forward target gap start.
    pub t_start: usize,
    /// Target gap length (bp).
    pub t_size: usize,
    /// Query contig name of enclosing fill.
    pub q_name: String,
    /// Query orientation of enclosing fill.
    pub q_strand: char,
    /// 0-based **forward** (`+`) query gap start.
    pub q_start: usize,
    /// Query gap length (bp).
    pub q_size: usize,
    /// Nested child fills placed inside this target gap.
    pub fills: Vec<UcscNetFill>,
}

/// Complete UCSC `.net` hierarchy for a single target contig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UcscTargetNet {
    /// Target contig name.
    pub t_name: String,
    /// Total target contig length (bp).
    pub t_size: usize,
    /// Disjoint top-level fills sorted by `t_start`.
    pub fills: Vec<UcscNetFill>,
}

/// Slice a chain's ungapped blocks strictly into `[t_lo, t_hi)`.
fn slice_chain_to_target(ch: &UcscChain, t_lo: usize, t_hi: usize) -> Vec<UcscChainBlock> {
    let mut out = Vec::new();
    for &b in &ch.blocks {
        let b_end = b.t_start + b.size;
        if b_end <= t_lo {
            continue;
        }
        if b.t_start >= t_hi {
            break;
        }
        let s = b.t_start.max(t_lo);
        let e = b_end.min(t_hi);
        if e > s {
            let off = s - b.t_start;
            out.push(UcscChainBlock {
                t_start: s,
                q_start: b.q_start + off,
                size: e - s,
            });
        }
    }
    out
}

#[derive(Clone)]
struct GapNode {
    t_start: usize,
    t_end: usize,
    q_name: String,
    q_strand: char,
    q_fwd_start: usize,
    q_fwd_end: usize,
    used_t: IntervalSet,
    fills: Vec<FillNode>,
}

#[derive(Clone)]
struct FillNode {
    t_start: usize,
    t_end: usize,
    q_name: String,
    q_strand: char,
    q_fwd_start: usize,
    q_fwd_end: usize,
    chain_id: u64,
    score: i64,
    ali: usize,
    q_over: usize,
    q_far: i64,
    q_dup: usize,
    fill_type: &'static str,
    gaps: Vec<GapNode>,
}

struct ParentGapCtx<'a> {
    name: &'a str,
    strand: char,
    fwd_start: usize,
    fwd_end: usize,
}

#[derive(Default)]
struct NetScratch {
    used_query: HashMap<String, IntervalSet>,
    merge_keys: Vec<usize>,
}

fn interval_overlap_len(a0: usize, a1: usize, b0: usize, b1: usize) -> usize {
    a1.min(b1).saturating_sub(a0.max(b0))
}

fn covered_bp_in_set(set: &IntervalSet, lo: usize, hi: usize) -> usize {
    if hi <= lo {
        return 0;
    }
    let mut free = Vec::new();
    set.free_subintervals(lo, hi, &mut free);
    let free_bp: usize = free.iter().map(|&(a, b)| b - a).sum();
    (hi - lo).saturating_sub(free_bp)
}

fn build_fill_node(
    ch: &UcscChain,
    blocks: &[UcscChainBlock],
    parent: Option<&ParentGapCtx<'_>>,
    scratch: &mut NetScratch,
) -> FillNode {
    let first = blocks[0];
    let last = blocks[blocks.len() - 1];
    let t_start = first.t_start;
    let t_end = last.t_start + last.size;
    let q_s_lo = first.q_start;
    let q_s_hi = last.q_start + last.size;
    let (q_fwd_start, q_fwd_end) = ch.strand_to_fwd(q_s_lo, q_s_hi);
    let q_size = q_fwd_end.saturating_sub(q_fwd_start);
    let ali: usize = blocks.iter().map(|b| b.size).sum();
    let full_ali = ch.aligned_bp().max(1);
    let score = if ali == full_ali {
        ch.score
    } else {
        let num = i128::from(ch.score) * (ali as i128);
        i64::try_from((num / (full_ali as i128)).max(1)).unwrap_or(ch.score)
    };

    let q_set = scratch.used_query.entry(ch.q_name.clone()).or_default();
    let q_dup = covered_bp_in_set(q_set, q_fwd_start, q_fwd_end);
    q_set.insert(q_fwd_start, q_fwd_end, &mut scratch.merge_keys);

    let (q_over, q_far, fill_type) = match parent {
        None => (0, -1, "top"),
        Some(pg) => {
            if pg.name == ch.q_name {
                let ov = interval_overlap_len(q_fwd_start, q_fwd_end, pg.fwd_start, pg.fwd_end);
                let dist = if ov > 0 {
                    0usize
                } else if q_fwd_end <= pg.fwd_start {
                    pg.fwd_start - q_fwd_end
                } else {
                    q_fwd_start.saturating_sub(pg.fwd_end)
                };
                let q_far = i64::try_from(dist).unwrap_or(i64::MAX);
                let ftype = if q_size > 0
                    && q_dup.saturating_mul(2) >= q_size
                    && dist > MAX_SYNTENIC_QFAR_BP
                {
                    "dup"
                } else if dist <= MAX_SYNTENIC_QFAR_BP {
                    if ch.q_strand == pg.strand {
                        "syn"
                    } else {
                        "inv"
                    }
                } else if q_size > 0 && q_dup.saturating_mul(2) >= q_size {
                    "dup"
                } else {
                    "nonSyn"
                };
                (ov, q_far, ftype)
            } else {
                (0, -1, "nonSyn")
            }
        }
    };

    let mut gaps = Vec::new();
    for pair in blocks.windows(2) {
        let b0 = pair[0];
        let b1 = pair[1];
        let gt0 = b0.t_start + b0.size;
        let gt1 = b1.t_start;
        if gt1 > gt0 {
            let gq0 = b0.q_start + b0.size;
            let gq1 = b1.q_start;
            let (qf0, qf1) = ch.strand_to_fwd(gq0, gq1);
            gaps.push(GapNode {
                t_start: gt0,
                t_end: gt1,
                q_name: ch.q_name.clone(),
                q_strand: ch.q_strand,
                q_fwd_start: qf0,
                q_fwd_end: qf1,
                used_t: IntervalSet::new(),
                fills: Vec::new(),
            });
        }
    }

    FillNode {
        t_start,
        t_end,
        q_name: ch.q_name.clone(),
        q_strand: ch.q_strand,
        q_fwd_start,
        q_fwd_end,
        chain_id: ch.id,
        score,
        ali,
        q_over,
        q_far,
        q_dup,
        fill_type,
        gaps,
    }
}

fn insert_chain_into_level(
    ch: &UcscChain,
    span: (usize, usize),
    depth: usize,
    used_t: &mut IntervalSet,
    fills: &mut Vec<FillNode>,
    parent: Option<&ParentGapCtx<'_>>,
    scratch: &mut NetScratch,
) {
    let (t_lo, t_hi) = span;
    if depth > MAX_NET_DEPTH || t_hi <= t_lo {
        return;
    }
    let mut free = Vec::new();
    used_t.free_subintervals(t_lo, t_hi, &mut free);

    // 1. Place unclaimed slices at this level
    let mut new_fills = Vec::new();
    for (fa, fb) in free {
        let blocks = slice_chain_to_target(ch, fa, fb);
        if blocks.is_empty() {
            continue;
        }
        let node = build_fill_node(ch, &blocks, parent, scratch);
        used_t.insert(node.t_start, node.t_end, &mut scratch.merge_keys);
        new_fills.push(node);
    }

    // 2. Recursively place any remainder into gaps of already-existing fills at this level
    if depth + 2 <= MAX_NET_DEPTH {
        for existing in fills.iter_mut() {
            if existing.t_end <= t_lo || existing.t_start >= t_hi {
                continue;
            }
            for gap in &mut existing.gaps {
                let gl = gap.t_start.max(t_lo);
                let gh = gap.t_end.min(t_hi);
                if gh > gl {
                    let ctx = ParentGapCtx {
                        name: &gap.q_name,
                        strand: gap.q_strand,
                        fwd_start: gap.q_fwd_start,
                        fwd_end: gap.q_fwd_end,
                    };
                    insert_chain_into_level(
                        ch,
                        (gl, gh),
                        depth + 2,
                        &mut gap.used_t,
                        &mut gap.fills,
                        Some(&ctx),
                        scratch,
                    );
                }
            }
        }
    }

    fills.extend(new_fills);
}

fn finalize_fill(mut node: FillNode) -> UcscNetFill {
    node.gaps.sort_unstable_by_key(|g| g.t_start);
    let gaps = node
        .gaps
        .into_iter()
        .filter(|g| !g.fills.is_empty() || g.t_end.saturating_sub(g.t_start) >= MIN_NET_GAP_BP)
        .map(|mut g| {
            g.fills.sort_unstable_by_key(|f| f.t_start);
            UcscNetGap {
                t_start: g.t_start,
                t_size: g.t_end - g.t_start,
                q_name: g.q_name,
                q_strand: g.q_strand,
                q_start: g.q_fwd_start,
                q_size: g.q_fwd_end.saturating_sub(g.q_fwd_start),
                fills: g.fills.into_iter().map(finalize_fill).collect(),
            }
        })
        .collect();
    UcscNetFill {
        t_start: node.t_start,
        t_size: node.t_end - node.t_start,
        q_name: node.q_name,
        q_strand: node.q_strand,
        q_start: node.q_fwd_start,
        q_size: node.q_fwd_end.saturating_sub(node.q_fwd_start),
        chain_id: node.chain_id,
        score: node.score,
        ali: node.ali,
        q_over: node.q_over,
        q_far: node.q_far,
        q_dup: node.q_dup,
        fill_type: node.fill_type,
        gaps,
    }
}

/// Build hierarchical UCSC `.net` trees from score-ranked [`UcscChain`]s.
///
/// `target_order` supplies `(t_name, t_size)` in output order; any target contig present in
/// `chains` with at least one fill is emitted.
#[must_use]
pub fn build_ucsc_nets(chains: &[UcscChain], target_order: &[(&str, usize)]) -> Vec<UcscTargetNet> {
    let mut by_target: HashMap<&str, Vec<&UcscChain>> = HashMap::new();
    for ch in chains {
        by_target.entry(ch.t_name.as_str()).or_default().push(ch);
    }
    for list in by_target.values_mut() {
        list.sort_unstable_by(|a, b| b.score.cmp(&a.score).then_with(|| a.id.cmp(&b.id)));
    }

    let mut ordered_targets: Vec<(&str, usize)> = Vec::with_capacity(target_order.len());
    let mut seen = HashSet::new();
    for &(name, size) in target_order {
        if by_target.contains_key(name) && seen.insert(name) {
            ordered_targets.push((name, size));
        }
    }
    let mut extras: Vec<(&str, usize)> = by_target
        .iter()
        .filter(|&(&name, _)| !seen.contains(name))
        .map(|(&name, list)| (name, list[0].t_size))
        .collect();
    extras.sort_unstable_by_key(|&(name, _)| name);
    ordered_targets.extend(extras);

    let mut scratch = NetScratch::default();
    let mut out = Vec::new();

    for (t_name, t_size) in ordered_targets {
        let Some(t_chains) = by_target.get(t_name) else {
            continue;
        };
        let mut used_t = IntervalSet::new();
        let mut fills: Vec<FillNode> = Vec::new();
        for &ch in t_chains {
            insert_chain_into_level(
                ch,
                (ch.t_start, ch.t_end),
                1,
                &mut used_t,
                &mut fills,
                None,
                &mut scratch,
            );
        }
        fills.sort_unstable_by_key(|f| f.t_start);
        if !fills.is_empty() {
            out.push(UcscTargetNet {
                t_name: t_name.to_owned(),
                t_size,
                fills: fills.into_iter().map(finalize_fill).collect(),
            });
        }
    }
    out
}

fn write_fill_rec<W: Write>(
    fill: &UcscNetFill,
    depth: usize,
    writer: &mut W,
) -> std::io::Result<()> {
    let indent = " ".repeat(depth);
    let mut line = String::with_capacity(128);
    let _ = write!(
        line,
        "{indent}fill {} {} {} {} {} {} id {} score {} ali {}",
        fill.t_start,
        fill.t_size,
        fill.q_name,
        fill.q_strand,
        fill.q_start,
        fill.q_size,
        fill.chain_id,
        fill.score,
        fill.ali,
    );
    if depth > 1 {
        let _ = write!(
            line,
            " qOver {} qFar {} qDup {}",
            fill.q_over, fill.q_far, fill.q_dup
        );
    } else if fill.q_dup > 0 {
        let _ = write!(line, " qDup {}", fill.q_dup);
    }
    let _ = writeln!(line, " type {}", fill.fill_type);
    writer.write_all(line.as_bytes())?;

    let gap_indent = " ".repeat(depth + 1);
    for gap in &fill.gaps {
        writeln!(
            writer,
            "{gap_indent}gap {} {} {} {} {} {}",
            gap.t_start, gap.t_size, gap.q_name, gap.q_strand, gap.q_start, gap.q_size
        )?;
        for child in &gap.fills {
            write_fill_rec(child, depth + 2, writer)?;
        }
    }
    Ok(())
}

/// Serialize target nets in standard UCSC `.net` format to `writer`.
///
/// # Errors
/// Returns an `std::io::Error` if writing to `writer` fails.
pub fn write_ucsc_nets<W: Write>(nets: &[UcscTargetNet], writer: &mut W) -> std::io::Result<()> {
    writeln!(writer, "# net generated by horizon_wga (UCSC net format)")?;
    for tn in nets {
        writeln!(writer, "net {} {}", tn.t_name, tn.t_size)?;
        for fill in &tn.fills {
            write_fill_rec(fill, 1, writer)?;
        }
    }
    Ok(())
}

/// Format target nets as a UCSC `.net` string.
#[must_use]
pub fn format_ucsc_nets(nets: &[UcscTargetNet]) -> String {
    let mut buf = Vec::new();
    let _ = write_ucsc_nets(nets, &mut buf);
    String::from_utf8(buf).unwrap_or_default()
}

/// Build UCSC chains (and optionally hierarchical nets) from aligned contig pairs and write them
/// to `chain_out` and/or `net_out` (each supporting `.gz` and `-` for stdout).
///
/// # Errors
/// Returns an error if creating or writing either output file fails.
#[allow(clippy::too_many_arguments)]
pub fn write_ucsc_chain_and_net_files<T, Q>(
    aligned_pairs: &[(usize, usize, &[AlignedRun])],
    names_x: &[&str],
    seqs_x: &[T],
    names_y: &[&str],
    seqs_y: &[Q],
    chain_out: Option<&std::path::Path>,
    net_out: Option<&std::path::Path>,
) -> Result<(usize, usize)>
where
    T: AsRef<[u8]>,
    Q: AsRef<[u8]>,
{
    let chains = build_ucsc_chains(aligned_pairs, names_x, seqs_x, names_y, seqs_y);
    if let Some(p) = chain_out {
        let mut w = crate::io_util::open_text_writer(Some(p.to_string_lossy().as_ref()))
            .context("creating UCSC chain output")?;
        write_ucsc_chains(&chains, &mut w)?;
        w.finish()?;
    }
    let mut net_fills = 0usize;
    if let Some(p) = net_out {
        let tgt_order: Vec<(&str, usize)> = names_x
            .iter()
            .zip(seqs_x.iter())
            .map(|(&n, s)| (n, s.as_ref().len()))
            .collect();
        let nets = build_ucsc_nets(&chains, &tgt_order);
        net_fills = nets.iter().map(|n| n.fills.len()).sum();
        let mut w = crate::io_util::open_text_writer(Some(p.to_string_lossy().as_ref()))
            .context("creating UCSC net output")?;
        write_ucsc_nets(&nets, &mut w)?;
        w.finish()?;
    }
    Ok((chains.len(), net_fills))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::revcomp;
    use crate::formats::ucsc::{format_ucsc_chains, parse_ucsc_chains};

    fn synth_dna(seed: u64, len: usize) -> Vec<u8> {
        let bases = *b"ACGT";
        let mut s = seed;
        (0..len)
            .map(|_| {
                s = crate::splitmix64(s.wrapping_add(1));
                bases[usize::try_from(s & 3).unwrap_or(0)]
            })
            .collect()
    }

    #[test]
    fn ucsc_chain_reverse_strand_coords_and_roundtrip() {
        let seq_x = synth_dna(0x1010, 400);
        let rc_x = revcomp(&seq_x);
        // Embed revcomp(seq_x[50..150]) at y[200..300] and revcomp(seq_x[180..260]) at y[90..170]
        let mut seq_y = synth_dna(0x2020, 500);
        seq_y[200..300].copy_from_slice(&rc_x[250..350]); // revcomp of x[50..150]
        seq_y[90..170].copy_from_slice(&rc_x[140..220]); // revcomp of x[180..260]

        let runs = vec![
            AlignedRun {
                x_start: 50,
                y_start: 299,
                len: 100,
                is_rc: true,
            },
            AlignedRun {
                x_start: 180,
                y_start: 169,
                len: 80,
                is_rc: true,
            },
        ];
        let chains = build_ucsc_chains(
            &[(0, 0, runs.as_slice())],
            &["chrT"],
            &[seq_x.as_slice()],
            &["chrQ"],
            &[seq_y.as_slice()],
        );
        assert_eq!(chains.len(), 1);
        let ch = &chains[0];
        assert_eq!(ch.q_strand, '-');
        assert_eq!((ch.t_start, ch.t_end), (50, 260));
        // In '-' UCSC chain coords: qStart = 500 - 300 = 200, qEnd = 500 - 90 = 410
        assert_eq!((ch.q_start, ch.q_end), (200, 410));
        assert_eq!(ch.q_fwd_interval(), (90, 300));
        assert_eq!(ch.blocks.len(), 2);

        let text = format_ucsc_chains(&chains);
        let parsed = parse_ucsc_chains(&text).expect("round-trip parse");
        assert_eq!(parsed, chains);

        // Net must emit forward (+) query coordinates [90, 300) => qStart=90 qSize=210
        let nets = build_ucsc_nets(&chains, &[("chrT", seq_x.len())]);
        assert_eq!(nets.len(), 1);
        let fill = &nets[0].fills[0];
        assert_eq!(fill.q_strand, '-');
        assert_eq!((fill.q_start, fill.q_size), (90, 210));
        assert_eq!(fill.ali, 180);
        assert_eq!(fill.gaps.len(), 1);
        assert_eq!((fill.gaps[0].t_start, fill.gaps[0].t_size), (150, 30));
        assert_eq!((fill.gaps[0].q_start, fill.gaps[0].q_size), (170, 30));
    }

    #[test]
    fn ucsc_net_nests_inversion_inside_syntenic_gap() {
        // Target: [left 300][inv 200][right 300] = 800 bp
        // Query:  [left 300][revcomp(inv) 200][right 300] = 800 bp
        let left = synth_dna(0xAAAA, 300);
        let mid = synth_dna(0xBBBB, 200);
        let right = synth_dna(0xCCCC, 300);
        let mut seq_x = Vec::new();
        seq_x.extend_from_slice(&left);
        seq_x.extend_from_slice(&mid);
        seq_x.extend_from_slice(&right);
        let mut seq_y = Vec::new();
        seq_y.extend_from_slice(&left);
        seq_y.extend_from_slice(&revcomp(&mid));
        seq_y.extend_from_slice(&right);

        let runs = vec![
            AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 300,
                is_rc: false,
            },
            AlignedRun {
                x_start: 300,
                y_start: 499,
                len: 200,
                is_rc: true,
            },
            AlignedRun {
                x_start: 500,
                y_start: 500,
                len: 300,
                is_rc: false,
            },
        ];
        let chains = build_ucsc_chains(
            &[(0, 0, runs.as_slice())],
            &["chrT"],
            &[seq_x.as_slice()],
            &["chrQ"],
            &[seq_y.as_slice()],
        );
        assert_eq!(
            chains.len(),
            2,
            "expected 1 outer '+' chain and 1 inner '-' chain"
        );
        assert_eq!(chains[0].q_strand, '+');
        assert_eq!(chains[1].q_strand, '-');

        let nets = build_ucsc_nets(&chains, &[("chrT", seq_x.len())]);
        assert_eq!(nets.len(), 1);
        let top = &nets[0].fills[0];
        assert_eq!(top.fill_type, "top");
        assert_eq!(top.gaps.len(), 1);
        assert_eq!(top.gaps[0].fills.len(), 1);
        let nested = &top.gaps[0].fills[0];
        assert_eq!(nested.fill_type, "inv");
        assert_eq!(nested.q_strand, '-');
        assert_eq!((nested.t_start, nested.t_size), (300, 200));
        assert_eq!((nested.q_start, nested.q_size), (300, 200));

        let net_text = format_ucsc_nets(&nets);
        assert!(net_text.contains("net chrT 800\n"));
        assert!(net_text.contains(" fill 0 800 chrQ + 0 800 id 1 "));
        assert!(net_text.contains("  gap 300 200 chrQ + 300 200\n"));
        assert!(net_text.contains("   fill 300 200 chrQ - 300 200 id 2 "));
        assert!(net_text.contains(" type inv\n"));
    }
}
