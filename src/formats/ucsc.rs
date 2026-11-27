//! Native UCSC `.chain` and `.net` builders, serializers, and chain parser.
//!
//! Coordinates follow the official UCSC Genome Browser specifications:
//! - In **`.chain`**, `tStrand` is always `+`, and when `qStrand == '-'`, `qStart`/`qEnd`
//!   are 0-based half-open reverse-complement coordinates (`qSize - qFwdEnd`, `qSize - qFwdStart`).
//! - In **`.net`**, `fill` and `gap` query intervals (`qStart`, `qSize`) are **always** in
//!   0-based forward (`+`) strand coordinates, even when `qStrand == '-'`.

use crate::mdl::AlignedRun;
use crate::seeding::complement_base;
use anyhow::{Context, Result, bail};
use std::io::Write;

pub use super::ucsc_net::{
    UcscNetFill, UcscNetGap, UcscTargetNet, build_ucsc_nets, format_ucsc_nets,
    write_ucsc_chain_and_net_files, write_ucsc_nets,
};

/// Maximum target or query gap (bp) bridged inside a single UCSC chain.
const MAX_CHAIN_GAP_BP: usize = 500_000;
/// Bounded predecessor lookback for collinear DP chaining.
const CHAIN_DP_LOOKBACK: usize = 64;
/// Minimum aligned bases required to emit a UCSC chain.
const MIN_CHAIN_ALI_BP: usize = 12;

/// One ungapped aligned block inside a [`UcscChain`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UcscChainBlock {
    /// 0-based forward start on the target sequence.
    pub t_start: usize,
    /// 0-based start on the query sequence in `q_strand` coordinates.
    pub q_start: usize,
    /// Ungapped block length in base pairs (`> 0`).
    pub size: usize,
}

/// One UCSC pairwise alignment chain (`chain` header + ungapped blocks).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UcscChain {
    /// Unique 1-based chain identifier (assigned descending by score).
    pub id: u64,
    /// Chain score (`>= 1`).
    pub score: i64,
    /// Target contig name.
    pub t_name: String,
    /// Total target contig length (bp).
    pub t_size: usize,
    /// 0-based half-open target span `[t_start, t_end)`.
    pub t_start: usize,
    /// 0-based half-open target end.
    pub t_end: usize,
    /// Query contig name.
    pub q_name: String,
    /// Total query contig length (bp).
    pub q_size: usize,
    /// Query orientation (`'+'` or `'-'`).
    pub q_strand: char,
    /// 0-based half-open query start in `q_strand` coordinates.
    pub q_start: usize,
    /// 0-based half-open query end in `q_strand` coordinates.
    pub q_end: usize,
    /// Strictly increasing, non-overlapping ungapped blocks (`dt >= 0`, `dq >= 0`, `dt + dq > 0`).
    pub blocks: Vec<UcscChainBlock>,
}

impl UcscChain {
    /// Total aligned ungapped bases across all blocks in this chain.
    #[must_use]
    pub fn aligned_bp(&self) -> usize {
        self.blocks.iter().map(|b| b.size).sum()
    }

    /// Convert a `[q_s, q_e)` interval in `q_strand` coordinates to forward (`+`) coordinates.
    #[must_use]
    pub fn strand_to_fwd(&self, q_s: usize, q_e: usize) -> (usize, usize) {
        if self.q_strand == '-' {
            (
                self.q_size.saturating_sub(q_e),
                self.q_size.saturating_sub(q_s),
            )
        } else {
            (q_s, q_e)
        }
    }

    /// Forward (`+`) strand query interval `[q_fwd_start, q_fwd_end)` spanned by this chain.
    #[must_use]
    pub fn q_fwd_interval(&self) -> (usize, usize) {
        self.strand_to_fwd(self.q_start, self.q_end)
    }
}

/// Piecewise-linear gap penalty modeled on UCSC `medium.lineargap` so large syntenic indels
/// remain chainable without falsely linking unrelated regions across megabases.
#[inline]
fn linear_gap_cost(dt: usize, dq: usize) -> i64 {
    fn one_dim(g: usize) -> i64 {
        if g == 0 {
            return 0;
        }
        let g = i64::try_from(g).unwrap_or(i64::MAX / 8);
        if g <= 100 {
            300 + g * 5
        } else if g <= 1_000 {
            800 + (g - 100) * 2
        } else if g <= 10_000 {
            2_600 + (g - 1_000) / 3
        } else {
            5_600 + (g - 10_000) / 20
        }
    }
    let single = one_dim(dt).saturating_add(one_dim(dq));
    // Slight double-sided gap surcharge when both target and query skip unaligned bases
    let both_bonus = if dt > 0 && dq > 0 { 250 } else { 0 };
    single.saturating_add(both_bonus)
}

#[derive(Clone, Copy)]
struct RawBlock {
    t_start: usize,
    q_start: usize,
    size: usize,
    score: i64,
}

impl RawBlock {
    #[inline]
    fn t_end(self) -> usize {
        self.t_start + self.size
    }
    #[inline]
    fn q_end(self) -> usize {
        self.q_start + self.size
    }
}

fn score_raw_block(run: AlignedRun, seq_x: &[u8], seq_y: &[u8]) -> Option<RawBlock> {
    let t_size = seq_x.len();
    let q_size = seq_y.len();
    if run.len == 0 || run.x_start >= t_size {
        return None;
    }
    let avail_x = t_size - run.x_start;
    let (len, q_strand_start) = if run.is_rc {
        if run.y_start >= q_size {
            return None;
        }
        let l = run.len.min(avail_x).min(run.y_start + 1);
        let q_rc_start = q_size - (run.y_start + 1);
        (l, q_rc_start)
    } else {
        if run.y_start >= q_size {
            return None;
        }
        let l = run.len.min(avail_x).min(q_size - run.y_start);
        (l, run.y_start)
    };
    if len == 0 {
        return None;
    }
    let mut score = 0i64;
    let sx = &seq_x[run.x_start..run.x_start + len];
    if run.is_rc {
        let y_lo = run.y_start + 1 - len;
        for (&bx, &by) in sx.iter().zip(seq_y[y_lo..=run.y_start].iter().rev()) {
            let ux = bx.to_ascii_uppercase();
            let uy = complement_base(by).to_ascii_uppercase();
            if matches!(ux, b'A' | b'C' | b'G' | b'T') && matches!(uy, b'A' | b'C' | b'G' | b'T') {
                score += if ux == uy { 100 } else { -100 };
            }
        }
    } else {
        let sy = &seq_y[run.y_start..run.y_start + len];
        for (&bx, &by) in sx.iter().zip(sy) {
            let ux = bx.to_ascii_uppercase();
            let uy = by.to_ascii_uppercase();
            if matches!(ux, b'A' | b'C' | b'G' | b'T') && matches!(uy, b'A' | b'C' | b'G' | b'T') {
                score += if ux == uy { 100 } else { -100 };
            }
        }
    }
    Some(RawBlock {
        t_start: run.x_start,
        q_start: q_strand_start,
        size: len,
        score: score.max(1),
    })
}

fn merge_adjacent_blocks(path: &[RawBlock]) -> (Vec<UcscChainBlock>, usize) {
    let mut merged: Vec<UcscChainBlock> = Vec::with_capacity(path.len());
    let mut total_ali = 0usize;
    for &b in path {
        total_ali += b.size;
        if let Some(last) = merged.last_mut()
            && last.t_start + last.size == b.t_start
            && last.q_start + last.size == b.q_start
        {
            last.size += b.size;
            continue;
        }
        merged.push(UcscChainBlock {
            t_start: b.t_start,
            q_start: b.q_start,
            size: b.size,
        });
    }
    (merged, total_ali)
}

fn chain_strand_blocks(
    mut blocks: Vec<RawBlock>,
    t_name: &str,
    t_size: usize,
    q_name: &str,
    q_size: usize,
    q_strand: char,
    out: &mut Vec<UcscChain>,
) {
    if blocks.is_empty() {
        return;
    }
    blocks.sort_unstable_by(|a, b| {
        a.t_start
            .cmp(&b.t_start)
            .then_with(|| a.q_start.cmp(&b.q_start))
            .then_with(|| b.size.cmp(&a.size))
    });
    blocks.dedup_by(|cur, prev| cur.t_start == prev.t_start && cur.q_start == prev.q_start);

    let n = blocks.len();
    let mut used = vec![false; n];
    let mut active: Vec<usize> = (0..n).collect();

    while !active.is_empty() {
        let m = active.len();
        let mut dp = vec![0i64; m];
        let mut prev = vec![usize::MAX; m];
        let mut best_idx = 0usize;
        let mut best_score = i64::MIN;

        for i in 0..m {
            let bi = blocks[active[i]];
            let mut s_best = bi.score;
            let mut p_best = usize::MAX;
            let lo = i.saturating_sub(CHAIN_DP_LOOKBACK);
            for j in (lo..i).rev() {
                let bj = blocks[active[j]];
                if bj.t_end() > bi.t_start {
                    continue;
                }
                let dt = bi.t_start - bj.t_end();
                if dt > MAX_CHAIN_GAP_BP {
                    break;
                }
                if bj.q_end() > bi.q_start {
                    continue;
                }
                let dq = bi.q_start - bj.q_end();
                if dq > MAX_CHAIN_GAP_BP {
                    continue;
                }
                let cand = dp[j]
                    .saturating_sub(linear_gap_cost(dt, dq))
                    .saturating_add(bi.score);
                if cand > s_best {
                    s_best = cand;
                    p_best = j;
                }
            }
            dp[i] = s_best;
            prev[i] = p_best;
            if s_best > best_score {
                best_score = s_best;
                best_idx = i;
            }
        }

        let mut path = Vec::new();
        let mut cur = best_idx;
        while cur != usize::MAX {
            let orig = active[cur];
            used[orig] = true;
            path.push(blocks[orig]);
            cur = prev[cur];
        }
        path.reverse();

        let (merged, total_ali) = merge_adjacent_blocks(&path);
        if total_ali >= MIN_CHAIN_ALI_BP && best_score > 0 && !merged.is_empty() {
            let first = merged[0];
            let last = merged[merged.len() - 1];
            out.push(UcscChain {
                id: 0,
                score: best_score.max(1),
                t_name: t_name.to_owned(),
                t_size,
                t_start: first.t_start,
                t_end: last.t_start + last.size,
                q_name: q_name.to_owned(),
                q_size,
                q_strand,
                q_start: first.q_start,
                q_end: last.q_start + last.size,
                blocks: merged,
            });
        }
        active.retain(|&idx| !used[idx]);
    }
}

/// Build score-ranked UCSC [`UcscChain`] records across all aligned `(target, query)` contig pairs.
///
/// Chains are sorted globally descending by score and assigned deterministic 1-based IDs.
#[must_use]
pub fn build_ucsc_chains<T, Q>(
    aligned_pairs: &[(usize, usize, &[AlignedRun])],
    names_x: &[&str],
    seqs_x: &[T],
    names_y: &[&str],
    seqs_y: &[Q],
) -> Vec<UcscChain>
where
    T: AsRef<[u8]>,
    Q: AsRef<[u8]>,
{
    let mut chains = Vec::new();
    for &(x_i, y_i, runs) in aligned_pairs {
        if runs.is_empty() {
            continue;
        }
        let (Some(&nx), Some(sx), Some(&ny), Some(sy)) = (
            names_x.get(x_i),
            seqs_x.get(x_i).map(AsRef::as_ref),
            names_y.get(y_i),
            seqs_y.get(y_i).map(AsRef::as_ref),
        ) else {
            continue;
        };
        let mut fwd = Vec::new();
        let mut rev = Vec::new();
        for &r in runs {
            if let Some(rb) = score_raw_block(r, sx, sy) {
                if r.is_rc {
                    rev.push(rb);
                } else {
                    fwd.push(rb);
                }
            }
        }
        chain_strand_blocks(fwd, nx, sx.len(), ny, sy.len(), '+', &mut chains);
        chain_strand_blocks(rev, nx, sx.len(), ny, sy.len(), '-', &mut chains);
    }

    chains.sort_unstable_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| b.aligned_bp().cmp(&a.aligned_bp()))
            .then_with(|| a.t_name.cmp(&b.t_name))
            .then_with(|| a.t_start.cmp(&b.t_start))
            .then_with(|| a.q_name.cmp(&b.q_name))
            .then_with(|| a.q_start.cmp(&b.q_start))
    });
    for (idx, ch) in chains.iter_mut().enumerate() {
        ch.id = (idx as u64) + 1;
    }
    chains
}

/// Serialize chains in standard UCSC `.chain` format to `writer`.
///
/// # Errors
/// Returns an `std::io::Error` if writing to `writer` fails.
pub fn write_ucsc_chains<W: Write>(chains: &[UcscChain], writer: &mut W) -> std::io::Result<()> {
    for ch in chains {
        writeln!(
            writer,
            "chain {} {} {} + {} {} {} {} {} {} {} {}",
            ch.score,
            ch.t_name,
            ch.t_size,
            ch.t_start,
            ch.t_end,
            ch.q_name,
            ch.q_size,
            ch.q_strand,
            ch.q_start,
            ch.q_end,
            ch.id
        )?;
        for (i, blk) in ch.blocks.iter().enumerate() {
            if i + 1 < ch.blocks.len() {
                let next = ch.blocks[i + 1];
                let dt = next.t_start.saturating_sub(blk.t_start + blk.size);
                let dq = next.q_start.saturating_sub(blk.q_start + blk.size);
                writeln!(writer, "{}\t{dt}\t{dq}", blk.size)?;
            } else {
                writeln!(writer, "{}", blk.size)?;
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Format chains as a UCSC `.chain` string.
#[must_use]
pub fn format_ucsc_chains(chains: &[UcscChain]) -> String {
    let mut buf = Vec::new();
    let _ = write_ucsc_chains(chains, &mut buf);
    String::from_utf8(buf).unwrap_or_default()
}

/// Parse UCSC `.chain` text into validated [`UcscChain`] records.
///
/// # Errors
/// Returns an error if a header or data line is malformed or violates UCSC coordinate invariants.
pub fn parse_ucsc_chains(text: &str) -> Result<Vec<UcscChain>> {
    let mut out = Vec::new();
    let mut cur: Option<(UcscChain, usize, usize)> = None;

    for (line_idx, raw) in text.lines().enumerate() {
        let line_no = line_idx + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            if let Some((ch, t_cur, q_cur)) = cur.take() {
                anyhow::ensure!(
                    t_cur == ch.t_end && q_cur == ch.q_end,
                    "UCSC chain {} span mismatch at line {line_no}",
                    ch.id
                );
                out.push(ch);
            }
            continue;
        }
        if let Some(rest) = line.strip_prefix("chain ") {
            if let Some((ch, t_cur, q_cur)) = cur.take() {
                anyhow::ensure!(t_cur == ch.t_end && q_cur == ch.q_end);
                out.push(ch);
            }
            let f: Vec<&str> = rest.split_whitespace().collect();
            anyhow::ensure!(
                f.len() == 12,
                "malformed chain header at line {line_no}: expected 12 fields after 'chain'"
            );
            let score: i64 = f[0].parse().context("chain score")?;
            let t_size: usize = f[2].parse().context("chain tSize")?;
            anyhow::ensure!(f[3] == "+", "chain tStrand must be '+' at line {line_no}");
            let t_start: usize = f[4].parse().context("chain tStart")?;
            let t_end: usize = f[5].parse().context("chain tEnd")?;
            let q_size: usize = f[7].parse().context("chain qSize")?;
            let q_strand = match f[8] {
                "+" => '+',
                "-" => '-',
                other => bail!("invalid chain qStrand `{other}` at line {line_no}"),
            };
            let q_start: usize = f[9].parse().context("chain qStart")?;
            let q_end: usize = f[10].parse().context("chain qEnd")?;
            let id: u64 = f[11].parse().context("chain id")?;
            anyhow::ensure!(
                t_start < t_end && t_end <= t_size && q_start < q_end && q_end <= q_size,
                "out-of-bounds chain coordinates at line {line_no}"
            );
            let ch = UcscChain {
                id,
                score,
                t_name: f[1].to_owned(),
                t_size,
                t_start,
                t_end,
                q_name: f[6].to_owned(),
                q_size,
                q_strand,
                q_start,
                q_end,
                blocks: Vec::new(),
            };
            cur = Some((ch, t_start, q_start));
        } else if let Some((ch, t_cur, q_cur)) = cur.as_mut() {
            let f: Vec<&str> = line.split_whitespace().collect();
            anyhow::ensure!(
                f.len() == 1 || f.len() == 3,
                "invalid chain block line at {line_no}"
            );
            let size: usize = f[0].parse().context("chain block size")?;
            anyhow::ensure!(size > 0, "zero-length chain block at line {line_no}");
            ch.blocks.push(UcscChainBlock {
                t_start: *t_cur,
                q_start: *q_cur,
                size,
            });
            *t_cur += size;
            *q_cur += size;
            if f.len() == 3 {
                let dt: usize = f[1].parse().context("chain dt")?;
                let dq: usize = f[2].parse().context("chain dq")?;
                *t_cur += dt;
                *q_cur += dq;
            } else if let Some((finished, tc, qc)) = cur.take() {
                anyhow::ensure!(
                    tc == finished.t_end && qc == finished.q_end,
                    "UCSC chain {} terminal coordinate mismatch at line {line_no}",
                    finished.id
                );
                out.push(finished);
            }
        } else {
            bail!("unexpected chain data line outside chain block at line {line_no}");
        }
    }
    if let Some((ch, t_cur, q_cur)) = cur.take() {
        anyhow::ensure!(t_cur == ch.t_end && q_cur == ch.q_end);
        out.push(ch);
    }
    Ok(out)
}
