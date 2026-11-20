//! BED and GFF3 feature parsing for `horizon_lift`.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub(crate) struct InputFeature {
    pub chrom: String,
    pub start: usize, // 0-based half-open
    pub end: usize,
    pub name: String,
    pub score: String,
    pub strand: char,
}

pub(crate) fn parse_features(path: &Path, is_gff: bool) -> Result<Vec<InputFeature>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut out = Vec::new();
    for (idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed == "track"
            || trimmed.starts_with("track ")
            || trimmed.starts_with("track\t")
        {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if is_gff {
            if cols.len() < 9 {
                continue;
            }
            let chrom = cols[0].to_string();
            let start_1: usize = cols[3]
                .parse()
                .with_context(|| format!("invalid GFF3 start at line {}", idx + 1))?;
            let end_1: usize = cols[4]
                .parse()
                .with_context(|| format!("invalid GFF3 end at line {}", idx + 1))?;
            let start = start_1.saturating_sub(1);
            let end = end_1;
            anyhow::ensure!(
                start <= end,
                "invalid GFF3 interval {start_1}..{end_1} at line {}",
                idx + 1
            );
            let score = cols[5].to_string();
            let strand = match cols[6].chars().next() {
                Some('-') => '-',
                Some('+') => '+',
                _ => '.',
            };
            let name = cols[8]
                .split(';')
                .find_map(|kv| kv.strip_prefix("ID=").or_else(|| kv.strip_prefix("Name=")))
                .unwrap_or(cols[2])
                .to_string();
            out.push(InputFeature {
                chrom,
                start,
                end,
                name,
                score,
                strand,
            });
        } else {
            if cols.len() < 3 {
                continue;
            }
            let chrom = cols[0].to_string();
            let start: usize = cols[1]
                .parse()
                .with_context(|| format!("invalid BED start at line {}", idx + 1))?;
            let end: usize = cols[2]
                .parse()
                .with_context(|| format!("invalid BED end at line {}", idx + 1))?;
            anyhow::ensure!(
                start <= end,
                "invalid BED interval {start}..{end} at line {}",
                idx + 1
            );
            let name = match cols.get(3).copied() {
                Some(n) if !n.is_empty() => n.to_string(),
                _ => ".".to_string(),
            };
            let score = cols.get(4).copied().unwrap_or("0").to_string();
            let strand = match cols.get(5).and_then(|s| s.chars().next()) {
                Some('-') => '-',
                Some('+') => '+',
                _ => '.',
            };
            out.push(InputFeature {
                chrom,
                start,
                end,
                name,
                score,
                strand,
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn write_temp(name: &str, body: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("{name}_{}.txt", std::process::id()));
        std::fs::write(&p, body).unwrap();
        p
    }

    #[test]
    fn parse_bed_fallbacks_track_header_and_unstranded() {
        let p = write_temp(
            "horizon_lift_bed_test",
            "# comment\ntrack name=demo\nchr1\t10\t20\nchr1\t30\t40\t\t100\t.\nchr1\t50\t60\tf1\t0\t-\n",
        );
        let feats = parse_features(&p, false).unwrap();
        let _ = std::fs::remove_file(&p);
        assert_eq!(feats.len(), 3);
        assert_eq!(feats[0].name, ".");
        assert_eq!(feats[0].score, "0");
        assert_eq!(feats[0].strand, '.');
        assert_eq!(feats[1].name, ".");
        assert_eq!(feats[1].score, "100");
        assert_eq!(feats[1].strand, '.');
        assert_eq!(feats[2].name, "f1");
        assert_eq!(feats[2].strand, '-');
    }

    #[test]
    fn parse_gff3_name_fallbacks_and_unstranded() {
        let p = write_temp(
            "horizon_lift_gff_test",
            "##gff-version 3\ntrack\nchr1\t.\texon\t1\t10\t.\t.\t.\tName=ex1\nchr1\t.\tgene\t11\t20\t.\t+\t.\tnote=none\n",
        );
        let feats = parse_features(&p, true).unwrap();
        let _ = std::fs::remove_file(&p);
        assert_eq!(feats.len(), 2);
        assert_eq!(feats[0].start, 0);
        assert_eq!(feats[0].end, 10);
        assert_eq!(feats[0].name, "ex1");
        assert_eq!(feats[0].strand, '.');
        assert_eq!(feats[1].name, "gene");
        assert_eq!(feats[1].strand, '+');
    }

    #[test]
    fn parse_features_rejects_inverted_interval() {
        let p = write_temp("horizon_lift_bad_interval", "chr1\t25\t10\n");
        let err = parse_features(&p, false).unwrap_err();
        let _ = std::fs::remove_file(&p);
        assert!(err.to_string().contains("invalid BED interval"));
    }
}
