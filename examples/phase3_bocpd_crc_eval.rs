//! Phase 3 S2 BOCPD & CRC exact library evaluator (invoked by `scripts/validation/run_phase3_*.py`).
//!
//! Subcommands:
//! - `bocpd-steady <hazard> <max_run> <true_gc> <burnin_bp> <seqs.fa>`:
//!   Evaluates `bocpd_gc_profile_into` steady-state RMSE over `[burnin_bp..]` for each FASTA record.
//! - `bocpd-switch <hazard> <max_run> <switch_pos> <gc_before> <gc_after> <seqs.fa>`:
//!   Evaluates midpoint-crossing lag (`pos - switch_pos`) and steady-state RMSE on each record.
//! - `crc-report <seq_x.fa> <seq_y.fa> <n_cal> <alpha> <base_floor>`:
//!   Runs `calibrate_null_crc_report` on real biological sequences and prints JSON diagnostics.
//! - `crc-fallthrough <n_cal> <base_floor>`:
//!   Audits `calibrate_crc_threshold` across `alpha` values including `< 1/(n_cal+1)` fallthrough.

use anyhow::{Context, Result, bail};
use horizon_aligner::bocpd::{BocpdScratch, bocpd_gc_profile_into};
use horizon_aligner::mdl::{HmmParams, calibrate_crc_threshold, calibrate_null_crc_report};
use horizon_aligner::read_fasta_contigs;
use std::env;
use std::path::Path;

fn parse_f64(s: &str, name: &str) -> Result<f64> {
    s.parse::<f64>().with_context(|| format!("parsing {name}"))
}

fn parse_usize(s: &str, name: &str) -> Result<usize> {
    s.parse::<usize>()
        .with_context(|| format!("parsing {name}"))
}

fn run_bocpd_steady(args: &[String]) -> Result<()> {
    anyhow::ensure!(
        args.len() == 5,
        "usage: bocpd-steady <hazard> <max_run> <true_gc> <burnin_bp> <seqs.fa>"
    );
    let hazard = parse_f64(&args[0], "hazard")?;
    let max_run = parse_usize(&args[1], "max_run")?;
    let true_gc = parse_f64(&args[2], "true_gc")?;
    let burnin = parse_usize(&args[3], "burnin_bp")?;
    let contigs = read_fasta_contigs(Path::new(&args[4]), 0)?;

    let mut scratch = BocpdScratch::default();
    let mut gc_est = Vec::new();
    let mut rows = Vec::with_capacity(contigs.len());
    for (id, seq) in &contigs {
        bocpd_gc_profile_into(seq, hazard, max_run, &mut scratch, &mut gc_est);
        let slice = &gc_est[burnin.min(gc_est.len())..];
        let mut sq_sum = 0.0_f64;
        let mut mean_sum = 0.0_f64;
        for &g in slice {
            let d = g - true_gc;
            sq_sum += d * d;
            mean_sum += g;
        }
        #[allow(clippy::cast_precision_loss)]
        let denom = slice.len().max(1) as f64;
        let rmse = (sq_sum / denom).sqrt();
        let mean_gc = mean_sum / denom;
        rows.push(format!(
            "{{\"id\":\"{id}\",\"len\":{},\"rmse\":{rmse:.6},\"mean_gc\":{mean_gc:.6}}}",
            seq.len()
        ));
    }
    println!("[{}]", rows.join(","));
    Ok(())
}

fn run_bocpd_switch(args: &[String]) -> Result<()> {
    anyhow::ensure!(
        args.len() == 6,
        "usage: bocpd-switch <hazard> <max_run> <switch_pos> <gc_before> <gc_after> <seqs.fa>"
    );
    let hazard = parse_f64(&args[0], "hazard")?;
    let max_run = parse_usize(&args[1], "max_run")?;
    let switch_pos = parse_usize(&args[2], "switch_pos")?;
    let gc_before = parse_f64(&args[3], "gc_before")?;
    let gc_after = parse_f64(&args[4], "gc_after")?;
    let contigs = read_fasta_contigs(Path::new(&args[5]), 0)?;

    let mid = f64::midpoint(gc_before, gc_after);
    let rising = gc_after > gc_before;
    let mut scratch = BocpdScratch::default();
    let mut gc_est = Vec::new();
    let mut rows = Vec::with_capacity(contigs.len());
    for (id, seq) in &contigs {
        bocpd_gc_profile_into(seq, hazard, max_run, &mut scratch, &mut gc_est);
        let lag = gc_est
            .iter()
            .enumerate()
            .skip(switch_pos)
            .find_map(|(pos, &g)| {
                let crossed = if rising { g >= mid } else { g <= mid };
                crossed.then_some(pos - switch_pos)
            });
        // Post-switch steady-state RMSE (100 bp after switch to end)
        let settle = (switch_pos + 100).min(gc_est.len());
        let tail = &gc_est[settle..];
        let sq: f64 = tail.iter().map(|&g| (g - gc_after).powi(2)).sum();
        #[allow(clippy::cast_precision_loss)]
        let post_rmse = (sq / tail.len().max(1) as f64).sqrt();
        let lag_json = lag.map_or_else(|| "null".to_owned(), |v| v.to_string());
        rows.push(format!(
            "{{\"id\":\"{id}\",\"lag_bp\":{lag_json},\"post_rmse\":{post_rmse:.6}}}"
        ));
    }
    println!("[{}]", rows.join(","));
    Ok(())
}

fn run_crc_report(args: &[String]) -> Result<()> {
    anyhow::ensure!(
        args.len() == 5,
        "usage: crc-report <seq_x.fa> <seq_y.fa> <n_cal> <alpha> <base_floor>"
    );
    let cx = read_fasta_contigs(Path::new(&args[0]), 0)?;
    let cy = read_fasta_contigs(Path::new(&args[1]), 0)?;
    let n_cal = parse_usize(&args[2], "n_cal")?;
    let alpha = parse_f64(&args[3], "alpha")?;
    let base_floor = parse_f64(&args[4], "base_floor")?;
    let sx = &cx.first().context("empty seq_x")?.1;
    let sy = &cy.first().context("empty seq_y")?.1;
    let params = HmmParams::default();
    let rep = calibrate_null_crc_report(sx, sy, &[], 64, n_cal, alpha, base_floor, &params, 0.15);
    println!(
        "{{\"threshold_bits\":{:.6},\"base_floor_bits\":{:.6},\"n_cal_bio\":{},\"median_null_bits\":{:.6},\"max_null_bits\":{:.6},\"used_synthetic_fallback\":{}}}",
        rep.threshold_bits,
        rep.base_floor_bits,
        rep.n_cal_bio,
        rep.median_null_bits,
        rep.max_null_bits,
        rep.used_synthetic_fallback,
    );
    Ok(())
}

fn run_crc_fallthrough(args: &[String]) -> Result<()> {
    anyhow::ensure!(
        args.len() == 2,
        "usage: crc-fallthrough <n_cal> <base_floor>"
    );
    let n_cal = parse_usize(&args[0], "n_cal")?;
    let base_floor = parse_f64(&args[1], "base_floor")?;
    #[allow(clippy::cast_precision_loss)]
    let scores: Vec<f64> = (1..=n_cal).map(|i| i as f64).collect();
    let mut rows = Vec::new();
    for alpha in [0.10, 0.07, 0.05, 0.02, 0.01, 0.005] {
        let lam = calibrate_crc_threshold(&scores, alpha, base_floor);
        #[allow(clippy::cast_precision_loss)]
        let min_feasible = 1.0 / (n_cal + 1) as f64;
        let hit_fallthrough = alpha < min_feasible;
        rows.push(format!(
            "{{\"alpha\":{alpha:.4},\"min_feasible_alpha\":{min_feasible:.6},\"hit_fallthrough\":{hit_fallthrough},\"lambda_hat\":{lam:.6}}}"
        ));
    }
    println!("[{}]", rows.join(","));
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    let Some((cmd, rest)) = args.split_first() else {
        bail!(
            "usage: phase3_bocpd_crc_eval <bocpd-steady|bocpd-switch|crc-report|crc-fallthrough> ..."
        );
    };
    match cmd.as_str() {
        "bocpd-steady" => run_bocpd_steady(rest),
        "bocpd-switch" => run_bocpd_switch(rest),
        "crc-report" => run_crc_report(rest),
        "crc-fallthrough" => run_crc_fallthrough(rest),
        other => bail!("unknown subcommand `{other}`"),
    }
}
