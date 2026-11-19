//! Benchmark binary running KP-14 + banded MDL/BOCPD on real Plasmodium Chr1 FASTAs.

use anyhow::{Context, Result};
use horizon_aligner::{MASK_PRIMARY_11_15, banded_mdl_gap, canonical_spaced_seeds, lis_chain};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn read_first_contig(path: &Path) -> Result<Vec<u8>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut seq = Vec::with_capacity(700_000);
    let mut started = false;
    for line in text.lines() {
        if line.starts_with('>') {
            if started {
                break;
            }
            started = true;
        } else if started {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    Ok(seq)
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        anyhow::bail!("usage: bench_chr1 <target.fa> <query.fa>");
    }
    let (fa_x, fa_y) = (PathBuf::from(&args[1]), PathBuf::from(&args[2]));

    let t0 = Instant::now();
    let seq_x = read_first_contig(&fa_x)?;
    let seq_y_raw = read_first_contig(&fa_y)?;

    let ix = canonical_spaced_seeds(&seq_x, MASK_PRIMARY_11_15, 1)?;
    let iy_fwd = canonical_spaced_seeds(&seq_y_raw, MASK_PRIMARY_11_15, 1)?;
    let mut hits_fwd = Vec::new();
    for (k, px) in &ix {
        if let Some(py) = iy_fwd.get(k) {
            hits_fwd.push((px[0], py[0]));
        }
    }
    let chain_fwd = lis_chain(&hits_fwd);

    let seq_y_rc = revcomp(&seq_y_raw);
    let iy_rc = canonical_spaced_seeds(&seq_y_rc, MASK_PRIMARY_11_15, 1)?;
    let mut hits_rc = Vec::new();
    for (k, px) in &ix {
        if let Some(py) = iy_rc.get(k) {
            hits_rc.push((px[0], py[0]));
        }
    }
    let chain_rc = lis_chain(&hits_rc);

    let (anchors, seq_y) = if chain_rc.len() > chain_fwd.len() {
        (chain_rc, seq_y_rc)
    } else {
        (chain_fwd, seq_y_raw)
    };

    let span = MASK_PRIMARY_11_15.len();
    let mut macro_bubbles = 0usize;
    let mut dp_gaps = 0usize;
    let mut kept_gap_cols = 0usize;
    let mut internal_bubbles = 0usize;

    let mut gap_spans = Vec::with_capacity(anchors.len() + 1);
    let mut prev_end = (0usize, 0usize);
    for &(x, y) in &anchors {
        gap_spans.push((prev_end, (x, y)));
        prev_end = (x + span, y + span);
    }
    gap_spans.push((prev_end, (seq_x.len(), seq_y.len())));

    for ((x0, y0), (x1, y1)) in gap_spans {
        if x1 <= x0 || y1 <= y0 {
            continue;
        }
        let lx = x1 - x0;
        let ly = y1 - y0;
        if lx > 800 || ly > 800 {
            macro_bubbles += 1;
            continue;
        }
        if lx >= 25 && ly >= 25 {
            dp_gaps += 1;
            let res = banded_mdl_gap(&seq_x[x0..x1], &seq_y[y0..y1], 48, 10.0);
            kept_gap_cols += res.kept_pairs.len();
            internal_bubbles += res.bubbles;
        }
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "len_x={} len_y={} anchors={} macro_bubbles={} dp_gaps={} kept_gap_cols={} internal_bubbles={} elapsed_s={:.4}",
        seq_x.len(),
        seq_y.len(),
        anchors.len(),
        macro_bubbles,
        dp_gaps,
        kept_gap_cols,
        internal_bubbles,
        elapsed
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_and_read_first_contig_work() {
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
        let fa_path = std::env::temp_dir().join("bench_chr1_test.fa");
        fs::write(&fa_path, ">c1\nACGT\nTTAA\n>c2\nGGGG\n").unwrap();
        let seq = read_first_contig(&fa_path).unwrap();
        let _ = fs::remove_file(&fa_path);
        assert_eq!(seq, b"ACGTTTAA");
    }
}
