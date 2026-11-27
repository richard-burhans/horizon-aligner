//! Streaming MDL + BOCPD + SDUST filter for PAF alignments (e.g. `FastGA` / `wfmash` outputs).

use anyhow::{Context, Result, bail};
use horizon_aligner::mdl::{DEFAULT_CRC_ALPHA, DEFAULT_CRC_N_CAL};
use horizon_aligner::paf_filter::{NetScoreMode, PafFilterConfig, run_paf_filter};
use std::env;
use std::time::Instant;

struct CliArgs {
    filter: PafFilterConfig,
    threads: Option<usize>,
}

fn next_arg_val<'a>(args: &'a [String], i: usize, flag: &str) -> Result<&'a str> {
    args.get(i + 1)
        .map(String::as_str)
        .with_context(|| format!("missing value for {flag}"))
}

#[allow(clippy::too_many_lines)]
fn parse_args() -> Result<Option<CliArgs>> {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut paf_path: Option<String> = None;
    let mut query_fa: Option<String> = None;
    let mut target_fa: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut threads: Option<usize> = None;
    let mut min_bits_override: Option<f64> = None;
    let mut crc_alpha = DEFAULT_CRC_ALPHA;
    let mut crc_n_cal = DEFAULT_CRC_N_CAL;
    let mut min_nonrepeat_bp = 0usize;
    let mut min_aln_len = 0usize;
    let mut keep_all = false;
    let mut max_occ: Option<usize> = None;
    let mut net_score = NetScoreMode::Mdl;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!(
                    "usage: horizon_filter --paf <in.paf[.gz]|-> [--query <q.fa>] [--target <t.fa>] [--seqs <all.fa>] [--out <out.paf>] [-t N | --threads N] [--min-bits <f64>] [--crc-alpha <f64>] [--crc-n-cal <usize>] [--min-nonrepeat-bp <usize>] [--min-len <usize>] [--net] [--max-occ <usize>] [--net-score mdl|log-length-ani] [--keep-all]"
                );
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("horizon_filter {}", env!("CARGO_PKG_VERSION"));
                return Ok(None);
            }
            "-t" | "--threads" => {
                let raw = next_arg_val(&args, i, "--threads")?;
                threads = Some(horizon_aligner::parse_threads_arg(raw)?);
                i += 2;
            }
            "--paf" | "-i" => {
                paf_path = Some(next_arg_val(&args, i, "--paf")?.to_owned());
                i += 2;
            }
            "--query" | "-q" => {
                query_fa = Some(next_arg_val(&args, i, "--query")?.to_owned());
                i += 2;
            }
            "--target" => {
                target_fa = Some(next_arg_val(&args, i, "--target")?.to_owned());
                i += 2;
            }
            "--seqs" | "-s" => {
                let s = next_arg_val(&args, i, "--seqs")?.to_owned();
                query_fa = Some(s.clone());
                target_fa = Some(s);
                i += 2;
            }
            "--out" | "-o" => {
                out_path = Some(next_arg_val(&args, i, "--out")?.to_owned());
                i += 2;
            }
            "--min-bits" => {
                let val: f64 = next_arg_val(&args, i, "--min-bits")?
                    .parse()
                    .context("parsing --min-bits")?;
                anyhow::ensure!(val.is_finite(), "--min-bits must be finite, got {val}");
                min_bits_override = Some(val);
                i += 2;
            }
            "--crc-alpha" => {
                let val: f64 = next_arg_val(&args, i, "--crc-alpha")?
                    .parse()
                    .context("parsing --crc-alpha")?;
                anyhow::ensure!(
                    val.is_finite() && val > 0.0 && val <= 1.0,
                    "--crc-alpha must be in (0, 1], got {val}"
                );
                crc_alpha = val;
                i += 2;
            }
            "--crc-n-cal" => {
                crc_n_cal = next_arg_val(&args, i, "--crc-n-cal")?
                    .parse()
                    .context("parsing --crc-n-cal")?;
                i += 2;
            }
            "--min-nonrepeat-bp" => {
                min_nonrepeat_bp = next_arg_val(&args, i, "--min-nonrepeat-bp")?
                    .parse()
                    .context("parsing --min-nonrepeat-bp")?;
                i += 2;
            }
            "--min-len" => {
                min_aln_len = next_arg_val(&args, i, "--min-len")?
                    .parse()
                    .context("parsing --min-len")?;
                i += 2;
            }
            "--keep-all" => {
                keep_all = true;
                i += 1;
            }
            "--net" => {
                max_occ = Some(1);
                i += 1;
            }
            "--max-occ" => {
                let k: usize = next_arg_val(&args, i, "--max-occ")?
                    .parse()
                    .context("parsing --max-occ")?;
                anyhow::ensure!(k >= 1, "--max-occ must be >= 1");
                max_occ = Some(k);
                i += 2;
            }
            "--net-score" => {
                let raw = next_arg_val(&args, i, "--net-score")?;
                net_score = match raw {
                    "mdl" => NetScoreMode::Mdl,
                    "log-length-ani" | "log_length_ani" => NetScoreMode::LogLengthAni,
                    other => {
                        bail!("invalid --net-score `{other}` (expected `mdl` or `log-length-ani`)")
                    }
                };
                i += 2;
            }
            other => bail!("unknown argument '{other}'"),
        }
    }

    anyhow::ensure!(
        !(keep_all && max_occ.is_some()),
        "--keep-all cannot be combined with --net or --max-occ"
    );
    let paf_path = paf_path.context("missing required --paf <in.paf[.gz]|->")?;
    Ok(Some(CliArgs {
        filter: PafFilterConfig {
            paf_path,
            query_fa,
            target_fa,
            out_path,
            min_bits_override,
            crc_alpha,
            crc_n_cal,
            min_nonrepeat_bp,
            min_aln_len,
            keep_all,
            max_occ,
            net_score,
        },
        threads,
    }))
}

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    horizon_aligner::init_rayon_threads(cli.threads)?;
    let t0 = Instant::now();
    let stats = run_paf_filter(&cli.filter)?;
    eprintln!(
        "[horizon_filter] total={} kept={} filtered_mdl={} filtered_len={} elapsed={:.2}s",
        stats.total,
        stats.kept,
        stats.filtered_mdl,
        stats.filtered_len,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use horizon_aligner::mdl::DEFAULT_CRC_N_CAL;
    use horizon_aligner::paf::PafRecord;
    use horizon_aligner::paf_filter::{
        IndexedSeq, NetScoreMode, ScoreResult, calibrate_paf_operator_crc, net_top_k_records,
        net_top_k_records_with_mode, score_from_counts_only, score_paf_record_with_floor,
    };
    use horizon_aligner::splitmix64;
    use std::collections::HashMap;

    #[test]
    fn mdl_filter_rejects_simple_repeat_and_accepts_complex_match() {
        let complex =
            b"ACGTAGTGCTAGCATCGATCGTAGCTAGCTGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCA";
        let repeat =
            b"ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let mut q_map = HashMap::new();
        q_map.insert("q_c".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        q_map.insert("q_r".to_owned(), IndexedSeq::from_dna(repeat.to_vec()));

        let mut t_map = HashMap::new();
        t_map.insert("t_c".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        t_map.insert("t_r".to_owned(), IndexedSeq::from_dna(repeat.to_vec()));

        let paf_c =
            PafRecord::parse_line("q_c\t80\t0\t80\t+\tt_c\t80\t0\t80\t80\t80\t60\tcg:Z:80=")
                .expect("paf");
        let paf_r =
            PafRecord::parse_line("q_r\t80\t0\t80\t+\tt_r\t80\t0\t80\t80\t80\t60\tcg:Z:80=")
                .expect("paf");

        let res_c = score_paf_record_with_floor(&paf_c, Some(&q_map), Some(&t_map), 10.0);
        let res_r = score_paf_record_with_floor(&paf_r, Some(&q_map), Some(&t_map), 10.0);
        assert!(res_c.mdl_bits > res_c.threshold_bits);
        assert_eq!(res_c.nonrepeat_bp, 80);
        assert!(res_r.mdl_bits <= 0.0, "repeat should be SDUST-masked");
        assert_eq!(res_r.nonrepeat_bp, 0);

        let mut state = 0xCAFE_BABE_1234_5678_u64;
        let bases = *b"ACGT";
        let bio_x: Vec<u8> = (0..2_048)
            .map(|_| {
                state = splitmix64(state);
                bases[(state & 3) as usize]
            })
            .collect();
        let bio_y: Vec<u8> = (0..2_048)
            .map(|_| {
                state = splitmix64(state);
                bases[(state & 3) as usize]
            })
            .collect();
        let ix_x = IndexedSeq::from_dna(bio_x);
        let ix_y = IndexedSeq::from_dna(bio_y);
        let crc_floor = calibrate_paf_operator_crc(&ix_x, &ix_y, DEFAULT_CRC_N_CAL, 0.07, 10.0);
        assert!(
            crc_floor >= 10.0 && crc_floor.is_finite(),
            "expected finite operator-matched CRC floor >= 10.0, got {crc_floor}"
        );
    }

    #[test]
    fn counts_only_does_not_treat_match_op_as_100pct_exact() {
        let paf_m =
            PafRecord::parse_line("q\t100\t0\t100\t+\tt\t100\t0\t100\t50\t100\t60\tcg:Z:100M")
                .expect("paf");
        let paf_eq =
            PafRecord::parse_line("q\t100\t0\t100\t+\tt\t100\t0\t100\t100\t100\t60\tcg:Z:100=")
                .expect("paf");
        let res_m = score_from_counts_only(&paf_m, 10.0);
        let res_eq = score_from_counts_only(&paf_eq, 10.0);
        assert!(res_m.mdl_bits < res_eq.mdl_bits);
    }

    #[test]
    fn corrupt_or_oob_cigar_with_fasta_is_rejected_not_counted() {
        let complex =
            b"ACGTAGTGCTAGCATCGATCGTAGCTAGCTGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCA";
        let mut q_map = HashMap::new();
        q_map.insert("q".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        let mut t_map = HashMap::new();
        t_map.insert("t".to_owned(), IndexedSeq::from_dna(complex.to_vec()));

        let paf_oob = PafRecord::parse_line("q\t80\t0\t80\t+\tt\t80\t0\t80\t80\t80\t60\tcg:Z:200=")
            .expect("paf");
        let res = score_paf_record_with_floor(&paf_oob, Some(&q_map), Some(&t_map), 10.0);
        assert_eq!(res.mdl_bits, 0.0);
        assert_eq!(res.nonrepeat_bp, 0);
    }

    #[test]
    fn counts_only_m_cigar_does_not_double_penalize_indels_and_resets_on_corrupt() {
        let paf_m_del =
            PafRecord::parse_line("q\t100\t0\t50\t+\tt\t100\t0\t60\t50\t60\t60\tcg:Z:50M10D")
                .expect("paf");
        let paf_eq_del =
            PafRecord::parse_line("q\t100\t0\t50\t+\tt\t100\t0\t60\t50\t60\t60\tcg:Z:50=10D")
                .expect("paf");
        let r_m = score_from_counts_only(&paf_m_del, 10.0);
        let r_eq = score_from_counts_only(&paf_eq_del, 10.0);
        assert_eq!(r_m.nonrepeat_bp, 50);
        assert!((r_m.mdl_bits - r_eq.mdl_bits).abs() < 1e-9);

        let paf_corrupt =
            PafRecord::parse_line("q\t100\t0\t50\t+\tt\t100\t0\t60\t50\t60\t60\tcg:Z:10D???")
                .expect("paf");
        let paf_no_cg =
            PafRecord::parse_line("q\t100\t0\t50\t+\tt\t100\t0\t60\t50\t60\t60\t").expect("paf");
        let r_bad = score_from_counts_only(&paf_corrupt, 10.0);
        let r_none = score_from_counts_only(&paf_no_cg, 10.0);
        assert!((r_bad.mdl_bits - r_none.mdl_bits).abs() < 1e-9);
    }

    #[test]
    fn net_top_k_records_enforces_max_occ_and_50pct_overlap() {
        let mk = |qs: usize, qe: usize, ts: usize, te: usize, score: f64| {
            let line = format!("q1\t1000\t{qs}\t{qe}\t+\tt1\t1000\t{ts}\t{te}\t100\t100\t60");
            let rec = PafRecord::parse_line(&line).expect("paf");
            let res = ScoreResult {
                mdl_bits: score,
                nonrepeat_bp: qe - qs,
                threshold_bits: 10.0,
            };
            (rec, res)
        };
        let mk_recs = || {
            vec![
                mk(0, 100, 0, 100, 200.0),
                mk(40, 140, 200, 300, 150.0),
                mk(60, 160, 400, 500, 120.0),
                mk(200, 300, 0, 100, 100.0),
            ]
        };
        assert!(net_top_k_records(mk_recs(), 0).is_empty());
        let kept1 = net_top_k_records(mk_recs(), 1);
        assert_eq!(kept1.len(), 2);
        assert_eq!(kept1[0].0.qstart, 0);
        assert_eq!(kept1[1].0.qstart, 60);

        let kept2 = net_top_k_records(mk_recs(), 2);
        assert_eq!(kept2.len(), 4);
    }

    #[test]
    fn net_log_length_ani_prefers_high_identity_ortholog_over_longer_paralog() {
        // True ortholog: 10 kb at 99.5% identity (higher per-bp MDL density, shorter total bp)
        let ortho = PafRecord::parse_line(
            "q1\t100000\t10000\t20000\t+\tt1\t100000\t10000\t20000\t9950\t10000\t60",
        )
        .expect("paf");
        let ortho_res = ScoreResult {
            mdl_bits: 19_500.0,
            nonrepeat_bp: 10_000,
            threshold_bits: 20.0,
        };
        // Palindromic paralog: 20 kb overlapping the same query locus at 91% identity (higher raw
        // cumulative MDL bits purely due to 2x length, lower per-bp identity).
        let paralog = PafRecord::parse_line(
            "q1\t100000\t5000\t25000\t-\tt1\t100000\t50000\t70000\t18200\t20000\t60",
        )
        .expect("paf");
        let paralog_res = ScoreResult {
            mdl_bits: 25_000.0,
            nonrepeat_bp: 20_000,
            threshold_bits: 20.0,
        };

        let kept_mdl = net_top_k_records_with_mode(
            vec![(ortho.clone(), ortho_res), (paralog.clone(), paralog_res)],
            1,
            NetScoreMode::Mdl,
        );
        assert_eq!(kept_mdl.len(), 1);
        assert_eq!(kept_mdl[0].0.strand, '-');

        let kept_ani = net_top_k_records_with_mode(
            vec![(ortho, ortho_res), (paralog, paralog_res)],
            1,
            NetScoreMode::LogLengthAni,
        );
        assert_eq!(kept_ani.len(), 1);
        assert_eq!(
            kept_ani[0].0.strand, '+',
            "log-length-ani netting must prefer the 99.5% ortholog over the 91% palindromic paralog"
        );
    }

    #[test]
    fn net_log_length_ani_reciprocal_best_blocks_deleted_region_paralog_fallthrough() {
        let mk = |qs: usize, qe: usize, ts: usize, te: usize, m: usize, len: usize, bits: f64| {
            let line = format!("q1\t100000\t{qs}\t{qe}\t+\tt1\t100000\t{ts}\t{te}\t{m}\t{len}\t60");
            let rec = PafRecord::parse_line(&line).expect("paf");
            let res = ScoreResult {
                mdl_bits: bits,
                nonrepeat_bp: len,
                threshold_bits: 10.0,
            };
            (rec, res)
        };
        // A: true ortholog q1:[0,10k) -> t1:[0,10k) at 99% ANI.
        // E: paralog of deleted region q1:[20k,30k) -> t1:[0,10k) at 93% ANI (wins query [20k,30k), loses target [0,10k)).
        // F: 2nd-best spurious paralog q1:[20k,30k) -> t1:[50k,60k) at 90% ANI (free target!).
        // 1-pass MDL netting rejects E on target and lets F slip through; 2-axis reciprocal-best
        // LogLengthAni lets E block F on query while A blocks E on target, retaining only A.
        let recs = vec![
            mk(0, 10_000, 0, 10_000, 9_900, 10_000, 19_000.0),
            mk(20_000, 30_000, 0, 10_000, 9_300, 10_000, 14_000.0),
            mk(20_000, 30_000, 50_000, 60_000, 9_000, 10_000, 12_000.0),
        ];
        let kept_mdl = net_top_k_records_with_mode(recs.clone(), 1, NetScoreMode::Mdl);
        assert_eq!(kept_mdl.len(), 2, "1-pass MDL falls through to F");
        let kept_ani = net_top_k_records_with_mode(recs, 1, NetScoreMode::LogLengthAni);
        assert_eq!(
            kept_ani.len(),
            1,
            "2-axis reciprocal-best log-length-ani must retain only true ortholog A"
        );
        assert_eq!(kept_ani[0].0.qstart, 0);
        assert_eq!(kept_ani[0].0.tstart, 0);
    }

    #[test]
    fn cigar_exceeding_paf_qend_tend_bounds_is_rejected() {
        let complex =
            b"ACGTAGTGCTAGCATCGATCGTAGCTAGCTGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCA";
        let mut q_map = HashMap::new();
        q_map.insert("q".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        let mut t_map = HashMap::new();
        t_map.insert("t".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        let paf = PafRecord::parse_line("q\t80\t0\t40\t+\tt\t80\t0\t40\t40\t40\t60\tcg:Z:60=")
            .expect("paf");
        let res = score_paf_record_with_floor(&paf, Some(&q_map), Some(&t_map), 10.0);
        assert_eq!(res.mdl_bits, 0.0);
        assert_eq!(res.nonrepeat_bp, 0);
    }

    #[test]
    fn subrange_contig_name_applies_start_offset() {
        let mut seq = b"ATATATATATATATATATATATATATATATATATATATAT".to_vec();
        seq.extend_from_slice(b"ACGTAGTGCTAGCATCGATCGTAGCTAGCTGATCGATCGA");
        let mut q_map = HashMap::new();
        q_map.insert("chr1".to_owned(), IndexedSeq::from_dna(seq.clone()));
        let mut t_map = HashMap::new();
        t_map.insert("chr1".to_owned(), IndexedSeq::from_dna(seq));

        let paf = PafRecord::parse_line(
            "chr1:40-80\t40\t0\t40\t+\tchr1:40-80\t40\t0\t40\t40\t40\t60\tcg:Z:40=",
        )
        .expect("paf");
        let res = score_paf_record_with_floor(&paf, Some(&q_map), Some(&t_map), 10.0);
        assert!(res.nonrepeat_bp >= 35);
        assert!(res.mdl_bits > 40.0);
    }
}
