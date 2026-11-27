//! Measure live peak RSS and wall time of the 100-way progressive streaming Fitch HAL writer
//! (`HalAlignment::write_multi_contig_progressive_hdf5_hal`), including crossing the 4 MiB
//! streaming Fitch window boundary, and optionally emit the .hal for C++ `halStats` validation.

use horizon_aligner::hal::HalAlignment;
use std::path::PathBuf;
use std::time::Instant;

fn read_proc_kb(key: &str) -> Option<u64> {
    let txt = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in txt.lines() {
        if let Some(rest) = line.strip_prefix(key) {
            return rest.split_whitespace().next()?.parse::<u64>().ok();
        }
    }
    None
}

fn read_proc_mb(key: &str) -> f64 {
    let kb = read_proc_kb(key).unwrap_or(0);
    f64::from(u32::try_from(kb).unwrap_or(u32::MAX)) / 1024.0
}

fn synth_dna(len: usize, seed: u64) -> Vec<u8> {
    const BASES: [u8; 4] = *b"ACGT";
    let mut x = seed ^ 0x9E37_79B9_7F4A_7C15;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        out.push(BASES[(x & 3) as usize]);
    }
    out
}

fn mutate_dna(base: &[u8], seed: u64, rate_pct: u64) -> Vec<u8> {
    const BASES: [u8; 4] = *b"ACGT";
    let mut x = seed ^ 0xA076_1D64_78BD_642F;
    let mut out = base.to_vec();
    for b in &mut out {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        if (x % 100) < rate_pct {
            *b = BASES[((x >> 8) & 3) as usize];
        }
    }
    out
}

fn main() -> std::io::Result<()> {
    let mut args = std::env::args().skip(1);
    let n_species: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(100);
    let seq_len: usize = args
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5_000_000);
    let out_path: PathBuf = args.next().map_or_else(
        || PathBuf::from("/tmp/horizon_fitch_100way.hal"),
        PathBuf::from,
    );

    let n_queries = n_species.saturating_sub(1);
    let ref_seq = synth_dna(seq_len, 42);
    let ref_contigs = vec![("chr1".to_owned(), ref_seq.clone())];
    let mut queries = Vec::with_capacity(n_queries);
    let mut blocks = Vec::with_capacity(n_queries);
    for q in 0..n_queries {
        let qseq = mutate_dna(&ref_seq, 1000 + q as u64, 5);
        queries.push((format!("Sp{:03}", q + 1), vec![("chr1".to_owned(), qseq)]));
        blocks.push((q, 0usize, 0usize, 0usize, 0usize, seq_len, false));
    }

    let rss_before_mb = read_proc_mb("VmRSS:");
    let hwm_before_mb = read_proc_mb("VmHWM:");

    let t0 = Instant::now();
    let (_newick, root_bots) = HalAlignment::write_multi_contig_progressive_hdf5_hal(
        &out_path,
        "Sp000",
        &ref_contigs,
        &queries,
        &blocks,
    )?;
    let elapsed_s = t0.elapsed().as_secs_f64();

    let rss_after_mb = read_proc_mb("VmRSS:");
    let hwm_after_mb = read_proc_mb("VmHWM:");
    let hal_bytes = std::fs::metadata(&out_path)?.len();

    println!(
        "{{\"n_species\":{},\"n_hal_genomes\":{},\"seq_len_bp\":{},\"root_bottom_segments\":{},\"hal_bytes\":{},\"elapsed_s\":{:.4},\"input_rss_before_mb\":{:.2},\"hwm_before_mb\":{:.2},\"rss_after_mb\":{:.2},\"peak_hwm_after_mb\":{:.2},\"writer_delta_hwm_mb\":{:.2},\"out_path\":\"{}\"}}",
        n_species,
        1 + 2 * n_queries,
        seq_len,
        root_bots,
        hal_bytes,
        elapsed_s,
        rss_before_mb,
        hwm_before_mb,
        rss_after_mb,
        hwm_after_mb,
        (hwm_after_mb - hwm_before_mb).max(0.0),
        out_path.display()
    );
    Ok(())
}
