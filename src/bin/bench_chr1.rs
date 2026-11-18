//! Benchmark binary running KP-14 + banded MDL/BOCPD on real Plasmodium Chr1 FASTAs.

use anyhow::{Context, Result};
use horizon_aligner::{MASK_KP14, banded_mdl_gap, canonical_spaced_seeds, lis_chain};
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
    let (fa_x, fa_y) = if args.len() >= 3 {
        (PathBuf::from(&args[1]), PathBuf::from(&args[2]))
    } else {
        let base = "/home/agent/.gemini/antigravity-cli/brain/c77b5ff5-88ad-43af-973a-c37dbadc186d/scratch/plasmodium_fastas";
        (
            Path::new(base).join("GCA_032713075.1__Plasmodium_falciparum.fa"),
            Path::new(base).join("GCA_900097025.1__Plasmodium_reichenowi.fa"),
        )
    };

    let t0 = Instant::now();
    let seq_x = read_first_contig(&fa_x)?;
    let seq_y_raw = read_first_contig(&fa_y)?;

    let ix = canonical_spaced_seeds(&seq_x, MASK_KP14, 1)?;
    let iy_fwd = canonical_spaced_seeds(&seq_y_raw, MASK_KP14, 1)?;
    let mut hits_fwd = Vec::new();
    for (k, px) in &ix {
        if let Some(py) = iy_fwd.get(k) {
            hits_fwd.push((px[0], py[0]));
        }
    }
    let chain_fwd = lis_chain(&hits_fwd);

    let seq_y_rc = revcomp(&seq_y_raw);
    let iy_rc = canonical_spaced_seeds(&seq_y_rc, MASK_KP14, 1)?;
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

    let span = MASK_KP14.len();
    let mut macro_bubbles = 0usize;
    let mut dp_gaps = 0usize;
    let mut kept_gap_cols = 0usize;
    let mut internal_bubbles = 0usize;

    let mut pts = Vec::with_capacity(anchors.len() + 2);
    pts.push((0usize, 0usize));
    for &(x, y) in &anchors {
        pts.push((x + span, y + span));
    }
    pts.push((seq_x.len(), seq_y.len()));

    for w in pts.windows(2) {
        let (x0, y0) = w[0];
        let (x1, y1) = w[1];
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
