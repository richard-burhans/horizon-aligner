//! Buffered line-streaming parser for Graphical Fragment Assembly (GFA v1.0 and v1.1).
//!
//! Provides a streaming reader capable of parsing multi-gigabyte (e.g. 12 GB HPRC)
//! pangenome graphs record-by-record without loading the full graph topology into RAM.

use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static STDIN_SPOOL_SEQ: AtomicU64 = AtomicU64::new(0);

/// RAII handle for a GFA input path that may be `"-"` (`stdin`).
///
/// When `path == "-"`, `stdin` is spooled once into an exclusively created temporary file
/// that is automatically deleted on drop, allowing multi-pass GFA tools to re-read `stdin`
/// safely without leaking temp files on early error returns.
pub struct SpoolableGfaInput {
    path: String,
    temp_file: Option<PathBuf>,
}

impl SpoolableGfaInput {
    /// Open or spool a GFA path. If `path == "-"`, copies `stdin` into a unique temp file.
    ///
    /// # Errors
    /// Returns an error if creating the temp file or reading `stdin` fails.
    pub fn new(path: &str) -> Result<Self> {
        if path != "-" {
            return Ok(Self {
                path: path.to_string(),
                temp_file: None,
            });
        }
        let pid = std::process::id();
        let temp_dir = std::env::temp_dir();
        let (tmp_path, mut file) = loop {
            let seq = STDIN_SPOOL_SEQ.fetch_add(1, Ordering::Relaxed);
            let candidate = temp_dir.join(format!("horizon_stdin_{pid}_{seq}.gfa"));
            let mut opts = OpenOptions::new();
            opts.write(true).create_new(true);
            #[cfg(unix)]
            {
                use std::os::unix::fs::OpenOptionsExt as _;
                opts.mode(0o600);
            }
            match opts.open(&candidate) {
                Ok(f) => break (candidate, f),
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
                Err(e) => return Err(e.into()),
            }
        };
        let spooled = Self {
            path: tmp_path.to_string_lossy().into_owned(),
            temp_file: Some(tmp_path),
        };
        let stdin = io::stdin();
        let mut locked = stdin.lock();
        io::copy(&mut locked, &mut file)?;
        file.flush()?;
        Ok(spooled)
    }

    /// Path to pass to single- or multi-pass GFA readers.
    #[must_use]
    pub fn path(&self) -> &str {
        &self.path
    }
}

impl Drop for SpoolableGfaInput {
    fn drop(&mut self) {
        if let Some(ref tmp) = self.temp_file {
            let _ = std::fs::remove_file(tmp);
        }
    }
}

/// Open a buffered writer to `out_path`, or to standard output when `out_path` is `None`.
///
/// # Errors
/// Returns an error if `out_path` cannot be created.
pub fn open_output_writer(out_path: Option<&str>) -> io::Result<io::BufWriter<Box<dyn io::Write>>> {
    let sink: Box<dyn io::Write> = match out_path {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(io::stdout().lock()),
    };
    Ok(io::BufWriter::new(sink))
}

/// Represents a GFA Segment (S line)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaSegment {
    pub name: String,
    pub sequence: Vec<u8>,
    pub tags: Option<String>,
}

impl GfaSegment {
    /// Return segment length in bp, reading `LN:i:` from tab-separated tags when sequence is `*`.
    #[must_use]
    pub fn seq_len(&self) -> usize {
        if self.sequence != b"*" {
            return self.sequence.len();
        }
        self.tags
            .as_deref()
            .into_iter()
            .flat_map(|t| t.split('\t'))
            .find_map(|tag| tag.strip_prefix("LN:i:")?.parse::<usize>().ok())
            .unwrap_or(0)
    }
}

fn open_reader_path(path: &Path) -> Result<Box<dyn BufRead>> {
    if path == Path::new("-") {
        return Ok(Box::new(BufReader::new(io::stdin())));
    }
    let file = File::open(path)?;
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz") || e.eq_ignore_ascii_case("bgz"));
    if is_gz {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Open a GFA/text path (plain, `.gz`/`.bgz`, or `"-"` for stdin) as a buffered line reader.
///
/// # Errors
/// Returns an error if the file cannot be opened.
pub fn open_line_reader(path: &str) -> Result<Box<dyn BufRead>> {
    open_reader_path(Path::new(path))
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
    pub tags: Option<String>,
}

/// Parse a GFA v1.0 path `segment_names` field such as `14+,15-` into `(segment_name, is_forward)` steps.
#[must_use]
pub fn parse_path_steps(seg_names: &str) -> Vec<(&str, bool)> {
    seg_names
        .split(',')
        .filter_map(|tok| {
            let tok = tok.trim();
            if let Some(name) = tok.strip_suffix('+') {
                (!name.is_empty()).then_some((name, true))
            } else if let Some(name) = tok.strip_suffix('-') {
                (!name.is_empty()).then_some((name, false))
            } else {
                None
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
///
/// # Errors
/// Returns an error if the file cannot be opened/decompressed or if any callback fails.
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
    let reader = open_reader_path(path.as_ref())?;
    stream_gfa_reader(reader, on_segment, on_link, on_path, on_walk)
}

/// A single parsed GFA v1.1 record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GfaRecord {
    Segment(GfaSegment),
    Link(GfaLink),
    Path(GfaPath),
    Walk(GfaWalk),
}

const MAX_GFA_LINE_BYTES: usize = 256 * 1024 * 1024;

/// Stream GFA records from a path with a single callback that may stop early via `ControlFlow::Break(())`.
///
/// # Errors
/// Returns an error if the file cannot be opened, reading fails, a line exceeds 256 MiB, or the callback errors.
pub fn stream_gfa_records<P, F>(path: P, on_record: F) -> Result<()>
where
    P: AsRef<Path>,
    F: FnMut(GfaRecord) -> Result<std::ops::ControlFlow<()>>,
{
    let reader = open_reader_path(path.as_ref())?;
    stream_gfa_records_reader(reader, on_record)
}

fn parse_walk_record(mut parts: std::str::Split<'_, char>) -> std::io::Result<GfaWalk> {
    let parse_coord = |tok: Option<&str>, field: &str| -> std::io::Result<usize> {
        match tok {
            Some("*") | None => Ok(0),
            Some(s) => s.parse::<usize>().map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("invalid GFA walk {field} coordinate: `{s}`"),
                )
            }),
        }
    };
    let sample = parts.next().unwrap_or("").to_string();
    let haplotype = parts.next().unwrap_or("").to_string();
    let seq_id = parts.next().unwrap_or("").to_string();
    let start = parse_coord(parts.next(), "start")?;
    let end = parse_coord(parts.next(), "end")?;
    let walk = parts.next().unwrap_or("").to_string();
    let rest: Vec<&str> = parts.collect();
    let tags = (!rest.is_empty()).then(|| rest.join("\t"));
    Ok(GfaWalk {
        sample,
        haplotype,
        seq_id,
        start,
        end,
        walk,
        tags,
    })
}

/// Stream GFA records from a buffered reader with a single callback that may stop early.
///
/// # Errors
/// Returns an error if reading fails, a line exceeds 256 MiB, or the callback errors.
pub fn stream_gfa_records_reader<R, F>(mut reader: R, mut on_record: F) -> Result<()>
where
    R: BufRead,
    F: FnMut(GfaRecord) -> Result<std::ops::ControlFlow<()>>,
{
    let mut line = String::new();
    loop {
        line.clear();
        let n =
            std::io::Read::take(&mut reader, MAX_GFA_LINE_BYTES as u64 + 1).read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.len() > MAX_GFA_LINE_BYTES && !line.ends_with('\n') {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "GFA line exceeds 256 MiB safety limit",
            )
            .into());
        }
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed == "H"
            || trimmed.starts_with("H\t")
        {
            continue;
        }

        let mut parts = trimmed.split('\t');
        let record_type = parts.next().unwrap_or("");

        let flow = match record_type {
            "S" => {
                let mut s_parts = trimmed.splitn(4, '\t');
                let _ = s_parts.next();
                let name = s_parts.next().unwrap_or("").to_string();
                let sequence = s_parts.next().unwrap_or("").as_bytes().to_vec();
                let tags = s_parts
                    .next()
                    .and_then(|t| (!t.is_empty()).then(|| t.to_string()));
                on_record(GfaRecord::Segment(GfaSegment {
                    name,
                    sequence,
                    tags,
                }))?
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
                on_record(GfaRecord::Link(GfaLink {
                    from,
                    from_orient,
                    to,
                    to_orient,
                    overlap,
                }))?
            }
            "P" => on_record(GfaRecord::Path(GfaPath {
                name: parts.next().unwrap_or("").to_string(),
                segment_names: parts.next().unwrap_or("").to_string(),
                overlaps: parts.next().unwrap_or("").to_string(),
            }))?,
            "W" => on_record(GfaRecord::Walk(parse_walk_record(parts)?))?,
            _ => std::ops::ControlFlow::Continue(()),
        };
        if flow.is_break() {
            break;
        }
    }

    Ok(())
}

/// Stream GFA records from any buffered reader.
///
/// # Errors
/// Returns an error if reading fails or if any callback returns an error.
pub fn stream_gfa_reader<R, FS, FL, FP, FW>(
    reader: R,
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
    stream_gfa_records_reader(reader, |rec| {
        match rec {
            GfaRecord::Segment(s) => on_segment(s)?,
            GfaRecord::Link(l) => on_link(l)?,
            GfaRecord::Path(p) => on_path(p)?,
            GfaRecord::Walk(w) => on_walk(w)?,
        }
        Ok(std::ops::ControlFlow::Continue(()))
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
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
        let gz_path = dir.join(format!(
            "horizon_test_stream_{}_{:?}.gfa.gz",
            std::process::id(),
            std::thread::current().id()
        ));
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

    #[test]
    fn seq_len_path_steps_break_and_invalid_walk_coord() {
        let seg = GfaSegment {
            name: "s1".to_owned(),
            sequence: b"*".to_vec(),
            tags: Some("LN:i:42".to_owned()),
        };
        assert_eq!(seg.seq_len(), 42);
        assert_eq!(
            parse_path_steps("a+,b-, ,c+"),
            vec![("a", true), ("b", false), ("c", true)]
        );

        let data = "H\tVN:Z:1.1\nS\ts1\tACGT\nS\ts2\tTGCA\n";
        let mut seen = 0usize;
        stream_gfa_records_reader(std::io::Cursor::new(data), |rec| {
            if matches!(rec, GfaRecord::Segment(_)) {
                seen += 1;
                return Ok(std::ops::ControlFlow::Break(()));
            }
            Ok(std::ops::ControlFlow::Continue(()))
        })
        .expect("break cleanly");
        assert_eq!(seen, 1);

        let bad = "W\tHG002\t0\tchr1\tbad\t10\t>s1\n";
        assert!(
            stream_gfa_records_reader(std::io::Cursor::new(bad), |_| {
                Ok(std::ops::ControlFlow::Continue(()))
            })
            .is_err()
        );
    }
}
