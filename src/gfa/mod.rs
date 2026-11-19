//! Buffered line-streaming parser for Graphical Fragment Assembly (GFA v1.0 and v1.1).
//!
//! Provides a streaming reader capable of parsing multi-gigabyte (e.g. 12 GB HPRC)
//! pangenome graphs record-by-record without loading the full graph topology into RAM.

use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Represents a GFA Segment (S line)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaSegment {
    pub name: String,
    pub sequence: Vec<u8>,
}

/// Represents a GFA Link (L line)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaLink {
    pub from: String,
    pub from_orient: u8, // '+' or '-'
    pub to: String,
    pub to_orient: u8,
    pub overlap: String,
}

/// Represents a GFA v1.0 Path (P line)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaPath {
    pub name: String,
    pub segment_names: String, // Comma separated, e.g. "14+,15-"
    pub overlaps: String,
}

/// Represents a GFA v1.1 Walk (W line), used by Minigraph-Cactus
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaWalk {
    pub sample: String,
    pub haplotype: String,
    pub seq_id: String,
    pub start: usize,
    pub end: usize,
    pub walk: String, // e.g. ">1>2<3>4"
}

/// Parse a GFA v1.0 path `segment_names` field such as `14+,15-` into `(segment_name, is_forward)` steps.
#[must_use]
pub fn parse_path_steps(seg_names: &str) -> Vec<(&str, bool)> {
    seg_names
        .split(',')
        .filter_map(|tok| {
            let tok = tok.trim();
            if tok.len() < 2 {
                return None;
            }
            let (name, orient) = tok.split_at(tok.len() - 1);
            match orient {
                "+" => Some((name, true)),
                "-" => Some((name, false)),
                _ => None,
            }
        })
        .collect()
}

/// Parse a GFA v1.1 walk string such as `>12<34>56` into `(segment_name, is_forward)` steps.
#[must_use]
pub fn parse_walk_steps(walk: &str) -> Vec<(&str, bool)> {
    let bytes = walk.as_bytes();
    let mut steps = Vec::new();
    let mut i = 0usize;
    while i < bytes.len() {
        let dir = bytes[i];
        if dir != b'>' && dir != b'<' {
            i += 1;
            continue;
        }
        let is_fwd = dir == b'>';
        let start = i + 1;
        let mut end = start;
        while end < bytes.len() && bytes[end] != b'>' && bytes[end] != b'<' {
            end += 1;
        }
        if end > start {
            steps.push((&walk[start..end], is_fwd));
        }
        i = end;
    }
    steps
}

/// Stream a GFA from a file path (`.gfa` or `.gfa.gz`) or `"-"` for `stdin`.
pub fn stream_gfa<P, FS, FL, FP, FW>(
    path: P,
    on_segment: FS,
    on_link: FL,
    on_path: FP,
    on_walk: FW,
) -> Result<()>
where
    P: AsRef<Path>,
    FS: FnMut(GfaSegment) -> Result<()>,
    FL: FnMut(GfaLink) -> Result<()>,
    FP: FnMut(GfaPath) -> Result<()>,
    FW: FnMut(GfaWalk) -> Result<()>,
{
    let path_ref = path.as_ref();
    if path_ref == Path::new("-") {
        let stdin = io::stdin();
        let reader = BufReader::new(stdin.lock());
        return stream_gfa_reader(reader, on_segment, on_link, on_path, on_walk);
    }

    let file = File::open(path_ref)?;
    let is_gz = path_ref.extension().is_some_and(|ext| ext == "gz");
    if is_gz {
        let reader = BufReader::new(MultiGzDecoder::new(file));
        stream_gfa_reader(reader, on_segment, on_link, on_path, on_walk)
    } else {
        let reader = BufReader::new(file);
        stream_gfa_reader(reader, on_segment, on_link, on_path, on_walk)
    }
}

/// Stream GFA records from any buffered reader.
pub fn stream_gfa_reader<R, FS, FL, FP, FW>(
    mut reader: R,
    mut on_segment: FS,
    mut on_link: FL,
    mut on_path: FP,
    mut on_walk: FW,
) -> Result<()>
where
    R: BufRead,
    FS: FnMut(GfaSegment) -> Result<()>,
    FL: FnMut(GfaLink) -> Result<()>,
    FP: FnMut(GfaPath) -> Result<()>,
    FW: FnMut(GfaWalk) -> Result<()>,
{
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with('H') {
            line.clear();
            continue;
        }

        let mut parts = trimmed.split('\t');
        let record_type = parts.next().unwrap_or("");

        match record_type {
            "S" => {
                let name = parts.next().unwrap_or("").to_string();
                let sequence = parts.next().unwrap_or("").as_bytes().to_vec();
                on_segment(GfaSegment { name, sequence })?;
            }
            "L" => {
                let from = parts.next().unwrap_or("").to_string();
                let from_orient = parts
                    .next()
                    .and_then(|s| s.as_bytes().first().copied())
                    .unwrap_or(b'+');
                let to = parts.next().unwrap_or("").to_string();
                let to_orient = parts
                    .next()
                    .and_then(|s| s.as_bytes().first().copied())
                    .unwrap_or(b'+');
                let overlap = parts.next().unwrap_or("").to_string();
                on_link(GfaLink {
                    from,
                    from_orient,
                    to,
                    to_orient,
                    overlap,
                })?;
            }
            "P" => {
                on_path(GfaPath {
                    name: parts.next().unwrap_or("").to_string(),
                    segment_names: parts.next().unwrap_or("").to_string(),
                    overlaps: parts.next().unwrap_or("").to_string(),
                })?;
            }
            "W" => {
                // W sample hap seq_id start end walk
                on_walk(GfaWalk {
                    sample: parts.next().unwrap_or("").to_string(),
                    haplotype: parts.next().unwrap_or("").to_string(),
                    seq_id: parts.next().unwrap_or("").to_string(),
                    start: parts.next().unwrap_or("0").parse().unwrap_or(0),
                    end: parts.next().unwrap_or("0").parse().unwrap_or(0),
                    walk: parts.next().unwrap_or("").to_string(),
                })?;
            }
            _ => {} // Ignore custom tags/lines like C, J, U, O
        }
        line.clear();
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::{Cursor, Write};

    #[test]
    fn streams_all_record_types_and_handles_empty_link_orients() {
        let gfa = "\
H\tVN:Z:1.1
# comment line
S\t1\tACGT
S\t2\tTTAA
L\t1\t+\t2\t-\t0M
L\t1\t\t2\t\t0M
P\tp1\t1+,2-\t0M
W\tGRCh38\t0\tchr22\t0\t8\t>1<2
";
        let mut segs = Vec::new();
        let mut links = Vec::new();
        let mut paths = Vec::new();
        let mut walks = Vec::new();

        stream_gfa_reader(
            Cursor::new(gfa.as_bytes()),
            |s| {
                segs.push(s);
                Ok(())
            },
            |l| {
                links.push(l);
                Ok(())
            },
            |p| {
                paths.push(p);
                Ok(())
            },
            |w| {
                walks.push(w);
                Ok(())
            },
        )
        .expect("stream_gfa_reader should succeed");

        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0].name, "1");
        assert_eq!(segs[0].sequence, b"ACGT");

        assert_eq!(links.len(), 2);
        assert_eq!(links[0].from_orient, b'+');
        assert_eq!(links[0].to_orient, b'-');
        // Empty orientation fields must not panic and safely default to '+'
        assert_eq!(links[1].from_orient, b'+');
        assert_eq!(links[1].to_orient, b'+');

        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].name, "p1");
        assert_eq!(paths[0].segment_names, "1+,2-");

        assert_eq!(walks.len(), 1);
        assert_eq!(walks[0].sample, "GRCh38");
        assert_eq!(walks[0].seq_id, "chr22");
        assert_eq!(walks[0].end, 8);
        assert_eq!(walks[0].walk, ">1<2");
    }

    #[test]
    fn streams_gzipped_gfa_file() {
        let dir = std::env::temp_dir();
        let gz_path = dir.join("horizon_test_stream.gfa.gz");
        {
            // Write two concatenated gzip members (simulating BGZF / multi-member gzip)
            let mut f = File::create(&gz_path).unwrap();
            {
                let mut enc1 = GzEncoder::new(&mut f, Compression::fast());
                enc1.write_all(b"H\tVN:Z:1.0\nS\tn1\tGATTACA\n").unwrap();
                enc1.finish().unwrap();
            }
            {
                let mut enc2 = GzEncoder::new(&mut f, Compression::fast());
                enc2.write_all(b"W\tCHM13\t0\tchr22\t0\t7\t>n1\n").unwrap();
                enc2.finish().unwrap();
            }
        }

        let mut seg_count = 0usize;
        let mut walk_count = 0usize;
        stream_gfa(
            &gz_path,
            |s| {
                assert_eq!(s.sequence, b"GATTACA");
                seg_count += 1;
                Ok(())
            },
            |_| Ok(()),
            |_| Ok(()),
            |w| {
                assert_eq!(w.sample, "CHM13");
                walk_count += 1;
                Ok(())
            },
        )
        .unwrap();
        let _ = std::fs::remove_file(&gz_path);
        assert_eq!(seg_count, 1);
        assert_eq!(walk_count, 1);
    }
}
