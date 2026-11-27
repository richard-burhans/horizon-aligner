//! GFA v1.1, MAF v1, FASTA I/O, and e-graph paralog resolution helpers.

pub mod hnex;
pub use hnex::{HnexArchive, HnexBlock, HnexBubble, HnexBubbleSpan, HnexRow, HnexSeqMeta};
pub mod ucsc;
mod ucsc_net;
pub use ucsc::{
    UcscChain, UcscChainBlock, UcscNetFill, UcscNetGap, UcscTargetNet, build_ucsc_chains,
    build_ucsc_nets, format_ucsc_chains, format_ucsc_nets, parse_ucsc_chains,
    write_ucsc_chain_and_net_files, write_ucsc_chains, write_ucsc_nets,
};

mod dup_resolve;
mod maf;

use crate::mdl::RunAlignment;
use crate::netting::partition_colinear_run_groups;
pub(crate) use dup_resolve::collect_windowed_duplications;
pub use dup_resolve::{
    ResolvedDuplications, prune_resolved_dup_clusters, resolve_alignment_duplications_mdl,
};
pub(crate) use maf::{MAX_MAF_LINE_BYTES, for_each_maf_stanza};
pub use maf::{MafRow, for_each_maf_stanza_reader, parse_maf_s_line, parse_two_row_maf_blocks};
use std::fmt::Write as _;

/// Read peak resident memory (`VmHWM` in KiB) from `/proc/self/status` on Linux.
#[must_use]
pub fn read_vm_hwm_kb() -> Option<usize> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            return rest.split_whitespace().next()?.parse::<usize>().ok();
        }
    }
    None
}

/// Check whether `name` equals `prefix` or starts with `prefix` followed by `.`, `#`, or `:`.
#[must_use]
pub fn matches_sample_prefix(name: &str, prefix: &str) -> bool {
    name == prefix
        || name.strip_prefix(prefix).is_some_and(|rest| {
            rest.starts_with('.') || rest.starts_with('#') || rest.starts_with(':')
        })
}

/// Strip a leading `prefix.` / `prefix#` / `prefix:` sample qualifier from `name` when present.
#[must_use]
pub fn strip_sample_prefix<'a>(name: &'a str, prefix: &str) -> &'a str {
    if let Some(rest) = name.strip_prefix(prefix)
        && (rest.starts_with('.') || rest.starts_with('#') || rest.starts_with(':'))
    {
        &rest[1..]
    } else {
        name
    }
}

fn emit_gfa_blocks(
    name_x: &str,
    name_y: &str,
    runs: &[(usize, usize, usize, usize, char)],
) -> String {
    let mut out = String::from("H\tVN:Z:1.1\n");
    for (idx, &(bx, by, blen, _x_end, bst)) in runs.iter().enumerate() {
        let seg_id = format!("{name_x}__{name_y}__blk_{idx}");
        let _ = writeln!(
            out,
            "S\t{seg_id}\t*\tLN:i:{blen}\tXS:i:{bx}\tYS:i:{by}\tST:Z:{bst}"
        );
        if idx > 0 && bx >= runs[idx - 1].3 {
            let prev_seg = format!("{name_x}__{name_y}__blk_{}", idx - 1);
            let prev_bst = runs[idx - 1].4;
            let gap_x = bx - runs[idx - 1].3;
            if gap_x >= 6 {
                let bub_id = format!("{name_x}__{name_y}__bub_{idx}");
                let _ = writeln!(out, "S\t{bub_id}\t*\tLN:i:{gap_x}\tHB:i:1");
                let _ = writeln!(out, "L\t{prev_seg}\t{prev_bst}\t{bub_id}\t+\t0M");
                let _ = writeln!(out, "L\t{bub_id}\t+\t{seg_id}\t{bst}\t0M");
            } else {
                let _ = writeln!(out, "L\t{prev_seg}\t{prev_bst}\t{seg_id}\t{bst}\t0M");
            }
        }
    }
    out
}

/// Format e-graph-resolved paralog annotations as GFA v1.1 comment/tag lines.
#[must_use]
pub fn format_gfa_paralog_tags(name_x: &str, name_y: &str, dups: &[String]) -> String {
    if dups.is_empty() {
        return String::new();
    }
    format!("# EG:Z:{name_x}__{name_y}\t{}\n", dups.join(","))
}

/// Emit a compact GFA v1.1 summary segment for a [`RunAlignment`].
#[must_use]
pub fn format_gfa_summary_runs(name_x: &str, name_y: &str, aln: &RunAlignment) -> String {
    let mut out = String::from("H\tVN:Z:1.1\n");
    let total_bubbles = aln.macro_bubbles + aln.internal_bubbles;
    let _ = writeln!(
        out,
        "S\t{name_x}__{name_y}\t*\tLN:i:{}\tLF:i:{}\tLR:i:{}\tHB:i:{total_bubbles}",
        aln.total_pairs, aln.lcbs_fwd, aln.lcbs_rc
    );
    out
}

/// Emit detailed GFA v1.1 block segments, Horizon-Bubble nodes, and links for a [`RunAlignment`].
#[must_use]
pub fn format_gfa_blocks_runs(name_x: &str, name_y: &str, aln: &RunAlignment) -> String {
    if aln.aligned_runs.is_empty() {
        return String::from("H\tVN:Z:1.1\n");
    }
    let blocks = partition_colinear_run_groups(&aln.aligned_runs);
    let runs: Vec<(usize, usize, usize, usize, char)> = blocks
        .into_iter()
        .map(|(grp, is_rev)| {
            let strand = if is_rev { '-' } else { '+' };
            let first = grp[0];
            let last = grp[grp.len() - 1];
            let total_bp: usize = grp.iter().map(|r| r.len).sum();
            let x_end = last.x_start.saturating_add(last.len);
            let y_lo = if is_rev {
                last.y_start.saturating_add(1).saturating_sub(last.len)
            } else {
                first.y_start
            };
            (first.x_start, y_lo, total_bp, x_end, strand)
        })
        .collect();
    emit_gfa_blocks(name_x, name_y, &runs)
}

mod realign;
pub(crate) use realign::materialize_run_group;

/// Group run-encoded alignments into colinear MAF blocks and emit standard MAF v1 text.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn format_maf_blocks_runs(
    name_x: &str,
    len_x: usize,
    name_y: &str,
    len_y: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &RunAlignment,
) -> String {
    format_maf_blocks_runs_avoiding_gaps(
        name_x, len_x, name_y, len_y, seq_x, seq_y, aln, None, None,
    )
}

/// Emit MAF v1 blocks for `aln`, refusing to bridge any gap overlapping `forbidden_x` or
/// `forbidden_y` (enforcing strict 0-bp target/query non-overlap under `--one-to-one`).
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn format_maf_blocks_runs_avoiding_gaps(
    name_x: &str,
    len_x: usize,
    name_y: &str,
    len_y: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &RunAlignment,
    forbidden_x: Option<&crate::coord::IntervalSet>,
    forbidden_y: Option<&crate::coord::IntervalSet>,
) -> String {
    let mut out = String::from("##maf version=1 scoring=horizon_mdl\n");
    if aln.aligned_runs.is_empty() {
        return out;
    }
    for (grp, is_rev) in crate::netting::partition_colinear_run_groups_avoiding_gaps(
        &aln.aligned_runs,
        forbidden_x,
        forbidden_y,
    ) {
        let Some(m) = materialize_run_group(len_x, len_y, seq_x, seq_y, &grp, is_rev) else {
            continue;
        };
        let ax = String::from_utf8_lossy(&m.ax);
        let ay = String::from_utf8_lossy(&m.ay);
        let _ = writeln!(out, "a score={}.0", m.score_bp);
        let _ = writeln!(out, "s {name_x} {} {} + {len_x} {ax}", m.x0, m.x_end - m.x0);
        let _ = writeln!(
            out,
            "s {name_y} {} {} {} {len_y} {ay}\n",
            m.y_start, m.y_span, m.y_strand
        );
    }
    out
}

/// Read up to `max_contigs` FASTA records `(contig_name, sequence_bytes)` from `path`
/// (`max_contigs == 0` reads all records).
///
/// # Errors
/// Returns an `std::io::Error` if `path` cannot be opened or read.
pub fn read_fasta_contigs(
    path: &std::path::Path,
    max_contigs: usize,
) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    let mut reader = crate::io_util::open_text_reader(path).map_err(std::io::Error::other)?;
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq = Vec::new();
    let mut line_buf = Vec::with_capacity(4096);
    loop {
        line_buf.clear();
        let n = std::io::BufRead::read_until(
            &mut std::io::Read::take(&mut reader, MAX_MAF_LINE_BYTES as u64 + 1),
            b'\n',
            &mut line_buf,
        )?;
        if n == 0 {
            break;
        }
        if line_buf.len() > MAX_MAF_LINE_BYTES {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "FASTA line exceeds 256 MiB safety limit",
            ));
        }
        let mut s = line_buf.as_slice();
        while s.last().is_some_and(u8::is_ascii_whitespace) {
            s = &s[..s.len() - 1];
        }
        while s.first().is_some_and(u8::is_ascii_whitespace) {
            s = &s[1..];
        }
        if let Some(rest) = s.strip_prefix(b">") {
            if !name.is_empty() {
                out.push((std::mem::take(&mut name), std::mem::take(&mut seq)));
                if max_contigs > 0 && out.len() >= max_contigs {
                    return Ok(out);
                }
            }
            let header = String::from_utf8_lossy(rest);
            let parsed_name = header.split_whitespace().next().unwrap_or("contig");
            if parsed_name.len() > 4096 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "FASTA contig header name exceeds 4096-byte safety limit ({} bytes)",
                        parsed_name.len()
                    ),
                ));
            }
            if parsed_name.bytes().any(|b| b < 0x20 || b == 0x7f) {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "FASTA contig header contains ASCII control character",
                ));
            }
            parsed_name.clone_into(&mut name);
        } else if !name.is_empty() && !s.is_empty() {
            if s.iter()
                .any(|&b| b.is_ascii_whitespace() || b < 0x20 || b == 0x7f)
            {
                for &b in s {
                    if b.is_ascii_whitespace() {
                        continue;
                    }
                    if b < 0x20 || b == 0x7f {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "FASTA sequence contains ASCII control character",
                        ));
                    }
                    seq.push(b);
                }
            } else {
                seq.extend_from_slice(s);
            }
        }
    }
    if !name.is_empty() && (max_contigs == 0 || out.len() < max_contigs) {
        out.push((name, seq));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::read_fasta_contigs;

    #[test]
    fn reads_fasta_contigs_with_limits_and_blank_lines() {
        let fa_path =
            std::env::temp_dir().join(format!("horizon_formats_fa_test_{}.fa", std::process::id()));
        std::fs::write(&fa_path, ">chr1 desc\nACGT\n\nTGCA\n>   \nGGGG\n").expect("write fa");
        let all = read_fasta_contigs(&fa_path, 0).expect("read all");
        assert_eq!(
            all,
            vec![
                ("chr1".to_owned(), b"ACGTTGCA".to_vec()),
                ("contig".to_owned(), b"GGGG".to_vec()),
            ]
        );
        let one = read_fasta_contigs(&fa_path, 1).expect("read one");
        assert_eq!(one, vec![("chr1".to_owned(), b"ACGTTGCA".to_vec())]);
        let _ = std::fs::remove_file(&fa_path);
    }

    #[test]
    fn maf_row_fwd_pos_and_parse_s_line_handle_both_strands_and_errors() {
        let (fwd, len_f) = super::parse_maf_s_line("s hg38.chr1 10 4 + 100 ACGT").expect("fwd");
        assert_eq!(len_f, 4);
        assert_eq!(fwd.fwd_pos(0), 10);
        assert_eq!(fwd.fwd_pos(3), 13);

        let (rev, len_r) = super::parse_maf_s_line("s mm10.chr2 10 4 - 100 ACGT").expect("rev");
        assert_eq!(len_r, 4);
        // UCSC '-' start=10, size=4, src_size=100 -> forward [86, 90), so idx 0 -> 89, idx 3 -> 86.
        assert_eq!(rev.fwd_pos(0), 89);
        assert_eq!(rev.fwd_pos(3), 86);

        assert!(super::parse_maf_s_line("a score=1").is_none());
        assert!(super::parse_maf_s_line("s chr1 0 4 ? 100 ACGT").is_none());
        assert!(super::parse_maf_s_line("s chr1 bad 4 + 100 ACGT").is_none());
    }

    #[test]
    fn sample_prefix_helpers_and_two_row_maf_blocks_work() {
        assert!(super::matches_sample_prefix("hg38.chr1", "hg38"));
        assert!(super::matches_sample_prefix("hg38#0#chr1", "hg38"));
        assert!(!super::matches_sample_prefix("hg38_alt.chr1", "hg38"));
        assert_eq!(super::strip_sample_prefix("hg38#0#chr1", "hg38"), "0#chr1");
        assert_eq!(super::strip_sample_prefix("hg38#chr1", "hg38"), "chr1");
        assert_eq!(super::strip_sample_prefix("hg38.chr1", "hg38"), "chr1");

        let maf = "a score=1\ns hg38.chr1 0 4 + 10 ACGT\ns mm10.chr1 2 4 - 20 TGCA\n";
        let blocks = super::parse_two_row_maf_blocks(maf);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].0.src, "hg38.chr1");
        assert_eq!(blocks[0].1.src, "mm10.chr1");
    }

    #[test]
    fn format_gfa_blocks_runs_emits_forward_start_for_minus_strand() {
        // R11 regression: format_gfa_blocks_runs emits 0-based forward start for '-' YS:i.
        let rc_aln = crate::mdl::RunAlignment {
            aligned_runs: vec![crate::mdl::AlignedRun {
                x_start: 10,
                y_start: 49,
                len: 20,
                is_rc: true,
            }],
            total_pairs: 20,
            lcbs_fwd: 0,
            lcbs_rc: 1,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let gfa = super::format_gfa_blocks_runs("chr1", "q1", &rc_aln);
        assert!(
            gfa.contains("YS:i:30\tST:Z:-"),
            "minus-strand GFA YS:i must be 0-based forward start (30), got:\n{gfa}"
        );
    }

    #[test]
    fn prune_resolved_dup_clusters_keeps_highest_identity_ortholog_and_prunes_paralog() {
        let mut seq_x = vec![b'A'; 2_500];
        let mut seq_y = vec![b'T'; 6_500];
        // Unique flanking region in B0 [100..200)
        for i in 100..200 {
            let b = b"ACGT"[i % 4];
            seq_x[i] = b;
            seq_y[i] = b;
        }
        // Duplicated locus [1000..1100): exact ortholog at y=1000, ~10% diverged paralog at y=5000.
        for i in 0..100 {
            let b = b"ACGTTGCATG"[i % 10];
            seq_x[1000 + i] = b;
            seq_y[1000 + i] = b; // 100% identity ortholog
            seq_y[5000 + i] = if i % 10 == 0 { b'G' } else { b }; // ~10% diverged paralog
        }

        let aln = crate::mdl::RunAlignment {
            aligned_runs: vec![
                crate::mdl::AlignedRun {
                    x_start: 100,
                    y_start: 100,
                    len: 100,
                    is_rc: false,
                },
                crate::mdl::AlignedRun {
                    x_start: 1000,
                    y_start: 5000,
                    len: 100,
                    is_rc: false,
                },
                crate::mdl::AlignedRun {
                    x_start: 1000,
                    y_start: 1000,
                    len: 100,
                    is_rc: false,
                },
            ],
            total_pairs: 300,
            lcbs_fwd: 3,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };

        let resolved = super::resolve_alignment_duplications_mdl(&aln, &seq_x, &seq_y);
        assert_eq!(resolved.dup_cluster_ids.len(), 1);

        let pruned = super::prune_resolved_dup_clusters(&aln, &resolved, &seq_x, &seq_y);
        assert_eq!(pruned.aligned_runs.len(), 2);
        assert_eq!(pruned.total_pairs, 200);
        assert_eq!(pruned.lcbs_fwd, 2);
        assert_eq!(pruned.aligned_runs[0].x_start, 100);
        assert_eq!(pruned.aligned_runs[1].x_start, 1000);
        assert_eq!(
            pruned.aligned_runs[1].y_start, 1000,
            "must keep 100%-identical ortholog (y=1000) and prune diverged paralog (y=5000)"
        );

        let maf = super::format_maf_blocks_runs(
            "ref",
            seq_x.len(),
            "qry",
            seq_y.len(),
            &seq_x,
            &seq_y,
            &pruned,
        );
        let n_blocks = maf.lines().filter(|l| l.starts_with("a score=")).count();
        assert_eq!(n_blocks, 2);
        assert!(!maf.contains("s qry 5000 "));
    }

    #[test]
    fn resolve_alignment_duplications_mdl_avoids_modulo_false_collision() {
        // Disjoint blocks [1050..1250) and [1750..1950) previously collided into fixed bin "B1".
        let seq_x: Vec<u8> = (0..2_500).map(|i| b"ACGT"[(i * 7 + 3) % 4]).collect();
        let seq_y = seq_x.clone();
        let aln = crate::mdl::RunAlignment {
            aligned_runs: vec![
                crate::mdl::AlignedRun {
                    x_start: 1050,
                    y_start: 1050,
                    len: 200,
                    is_rc: false,
                },
                crate::mdl::AlignedRun {
                    x_start: 1750,
                    y_start: 1750,
                    len: 200,
                    is_rc: false,
                },
            ],
            total_pairs: 400,
            lcbs_fwd: 2,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let resolved = super::resolve_alignment_duplications_mdl(&aln, &seq_x, &seq_y);
        assert!(
            resolved.dup_cluster_ids.is_empty(),
            "disjoint blocks in the same 1kb modulo bin must NOT collide as duplications: {:?}",
            resolved.dup_labels
        );
    }

    #[test]
    fn resolve_alignment_duplications_mdl_detects_boundary_straddling_duplication() {
        // Two copies straddling the 1000 bp boundary: [990..1890) and [1010..1910).
        let mut seq_x: Vec<u8> = (0..3_000).map(|i| b"ACGT"[(i * 5 + 1) % 4]).collect();
        let mut seq_y = vec![b'A'; 6_000];
        // Ortholog at y=990 (100% identical to x[990..1890))
        seq_y[990..1890].copy_from_slice(&seq_x[990..1890]);
        // Paralog at y=4010 (5% diverged from x[1010..1910))
        for i in 0..900 {
            let b = seq_x[1010 + i];
            seq_y[4010 + i] = if i % 20 == 0 {
                if b == b'A' { b'C' } else { b'A' }
            } else {
                b
            };
        }
        // Keep seq_x mutable warning free
        seq_x[0] = b'A';

        let aln = crate::mdl::RunAlignment {
            aligned_runs: vec![
                crate::mdl::AlignedRun {
                    x_start: 990,
                    y_start: 990,
                    len: 900,
                    is_rc: false,
                },
                crate::mdl::AlignedRun {
                    x_start: 1010,
                    y_start: 4010,
                    len: 900,
                    is_rc: false,
                },
            ],
            total_pairs: 1_800,
            lcbs_fwd: 2,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };

        let resolved = super::resolve_alignment_duplications_mdl(&aln, &seq_x, &seq_y);
        assert_eq!(
            resolved.dup_cluster_ids.len(),
            1,
            "overlapping copies straddling 1000 bp must cluster and collapse under MDL: {:?}",
            resolved.dup_labels
        );
        let pruned = super::prune_resolved_dup_clusters(&aln, &resolved, &seq_x, &seq_y);
        assert_eq!(pruned.aligned_runs.len(), 1);
        assert_eq!(
            pruned.aligned_runs[0].y_start, 990,
            "must keep 100%-identical ortholog and drop diverged paralog"
        );
    }

    #[test]
    fn resolve_alignment_duplications_mdl_rejects_high_divergence_spurious_overlap() {
        // Copy 1 is 100% identical; Copy 2 is 75% diverged (JC-saturated noise -> must NOT collapse).
        let seq_x: Vec<u8> = (0..2_000).map(|i| b"ACGT"[i % 4]).collect();
        let mut seq_y = vec![b'A'; 4_000];
        seq_y[500..1000].copy_from_slice(&seq_x[500..1000]);
        for i in 0..500 {
            let b = seq_x[500 + i];
            // 75% mismatches (3 of every 4 bases mutated)
            seq_y[2500 + i] = if i % 4 == 0 {
                b
            } else {
                match b {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                }
            };
        }
        let aln = crate::mdl::RunAlignment {
            aligned_runs: vec![
                crate::mdl::AlignedRun {
                    x_start: 500,
                    y_start: 500,
                    len: 500,
                    is_rc: false,
                },
                crate::mdl::AlignedRun {
                    x_start: 500,
                    y_start: 2500,
                    len: 500,
                    is_rc: false,
                },
            ],
            total_pairs: 1_000,
            lcbs_fwd: 2,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let resolved = super::resolve_alignment_duplications_mdl(&aln, &seq_x, &seq_y);
        assert!(
            resolved.dup_cluster_ids.is_empty(),
            "75% diverged overlap must NOT be collapsed by MDL: {:?}",
            resolved.dup_labels
        );
    }
}
