//! Dump live Rust BOCPD and Pair-HMM Forward-Backward posterior traces as JSON
//! for Figure 4.1 (`scripts/figures/plot_pair_hmm_mdl.py`, R9-FIG4-LIVE).

use horizon_aligner::{
    GapPosteriorTrace, HmmParams, banded_mdl_gap_posterior_trace, bocpd_gc_profile,
};

fn print_trace_json(label: &str, t: &GapPosteriorTrace, trailing_comma: bool) {
    println!("  \"{label}\": {{");
    println!("    \"n\": {},", t.n);
    println!("    \"m\": {},", t.m);
    println!("    \"band_w\": {},", t.band_w);
    println!("    \"score_bits\": {:.6},", t.score_bits);
    println!("    \"bubbles\": {},", t.bubbles);
    println!("    \"max_block_bits\": {:.6},", t.max_block_bits);
    println!("    \"kept_count\": {},", t.kept_pairs.len());
    print!("    \"posterior_cells\": [");
    for (idx, &(i, j, p)) in t.posterior_cells.iter().enumerate() {
        if idx > 0 {
            print!(",");
        }
        print!("[{i},{j},{p:.6}]");
    }
    println!("],");
    print!("    \"mea_cols\": [");
    for (idx, &(xi, yj, bits, post, abst)) in t.mea_cols.iter().enumerate() {
        if idx > 0 {
            print!(",");
        }
        let xs = xi.map_or_else(|| "null".to_owned(), |v| v.to_string());
        let ys = yj.map_or_else(|| "null".to_owned(), |v| v.to_string());
        print!("[{xs},{ys},{bits:.6},{post:.6},{abst}]");
    }
    println!("]");
    if trailing_comma {
        println!("  }},");
    } else {
        println!("  }}");
    }
}

fn main() {
    // Panel A: 130 bp composition shift (30 bp warm-up + 44 bp 20% GC + 33 bp 58% GC + 23 bp 24% GC)
    let mut seq_a = Vec::with_capacity(130);
    for i in 0..30 {
        seq_a.push(if i % 5 == 0 { b'G' } else { b'A' });
    }
    for i in 0..44 {
        seq_a.push(if i % 5 == 2 { b'G' } else { b'A' });
    }
    for i in 0..33 {
        seq_a.push(if (i * 7) % 12 < 7 { b'G' } else { b'A' });
    }
    for i in 0..23 {
        seq_a.push(if i % 4 == 1 { b'G' } else { b'A' });
    }
    let gc_a = bocpd_gc_profile(&seq_a, 1.0 / 120.0, 80);

    let params = HmmParams::default();

    // Shared 48x48 live sequence pair (band_w=8) for Panels C (2D banded posterior heatmap)
    // and D (1D MEA posterior profile & Horizon Bubble abstention):
    // 16 bp conserved flank + 16 bp divergent middle + 16 bp conserved flank.
    let pair_x = b"GCTAGCATCGAATTGGCCCCCAAAACCCCAAACGATGCATTCGGAATT";
    let pair_y = b"GCTAGCATCGAATTGGTTTTGGGGTTTTGGGGCGATGCATTCGGAATT";
    let trace_cd = banded_mdl_gap_posterior_trace(pair_x, pair_y, 8, 10.0, &params, 0.15);

    println!("{{");
    print!("  \"panel_a_gc\": [");
    for (idx, &v) in gc_a[30..].iter().enumerate() {
        if idx > 0 {
            print!(",");
        }
        print!("{v:.6}");
    }
    println!("],");
    print_trace_json("panel_cd", &trace_cd, false);
    println!("}}");
}
