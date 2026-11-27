//! Buffered line-streaming parser for Graphical Fragment Assembly (GFA v1.0 and v1.1).
//!
//! Provides a streaming reader capable of parsing multi-gigabyte (e.g. 12 GB HPRC)
//! pangenome graphs record-by-record without loading the full graph topology into RAM.

use crate::io_util::{open_text_reader, open_text_writer};
use anyhow::Result;
use std::fs::OpenOptions;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static STDIN_SPOOL_SEQ: AtomicU64 = AtomicU64::new(0);

/// RAII handle for a GFA input path that may be `"-"` (`stdin`).
///
/// When `path == "-"`, `stdin` is spooled once into an exclusively created temporary file
/// that is automatically deleted on drop, allowing multi-pass GFA tools to re-read `stdin`
/// safely without leaking temp files on early error returns.
pub struct SpoolableGfaInput {
    path: PathBuf,
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
                path: PathBuf::from(path),
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
            path: tmp_path.clone(),
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
    pub fn path(&self) -> &Path {
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

/// Open a buffered writer to `out_path` (plain, `.gz`, or `.zst`), or to standard output when `out_path` is `None`.
///
/// # Errors
/// Returns an error if `out_path` cannot be created.
pub fn open_output_writer(out_path: Option<&str>) -> io::Result<crate::io_util::TextWriter> {
    open_text_writer(out_path).map_err(io::Error::other)
}

/// Represents a GFA Segment (S line)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaSegment {
    pub name: String,
    pub sequence: Vec<u8>,
    pub tags: Option<String>,
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

/// Parse a GFA v1.0/v1.1 link overlap CIGAR string (e.g. `"0M"`, `"12M"`, `"*"`) into overlap bp.
#[must_use]
pub fn parse_cigar_overlap_bp(cigar: &str) -> usize {
    if cigar.is_empty() || cigar == "*" || cigar == "0M" {
        return 0;
    }
    let mut ov = 0usize;
    let mut num = 0usize;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            num = num.saturating_mul(10).saturating_add(usize::from(b - b'0'));
        } else {
            if matches!(b, b'M' | b'=' | b'X') {
                ov = ov.saturating_add(num);
            }
            num = 0;
        }
    }
    ov
}

/// Split a `PanSN`-style `sample#hap#contig` GFA `P`-line path name into `(sample, hap, contig)`,
/// falling back to `(name, "0", name)` when not in 3-field `PanSN` format.
#[must_use]
pub fn parse_path_sample_fields(name: &str) -> (&str, &str, &str) {
    let mut parts = name.splitn(3, '#');
    match (parts.next(), parts.next(), parts.next()) {
        (Some(sample), Some(hap), Some(contig))
            if !sample.is_empty() && !hap.is_empty() && !contig.is_empty() =>
        {
            (sample, hap, contig)
        }
        _ => (name, "0", name),
    }
}

/// A single parsed GFA v1.1 record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GfaRecord {
    Segment(GfaSegment),
    Link(GfaLink),
    Path(GfaPath),
    Walk(GfaWalk),
}

/// Zero-copy borrowed view of a GFA Segment (`S` line).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GfaSegmentRef<'a> {
    pub name: &'a str,
    pub sequence: &'a [u8],
    pub tags: Option<&'a str>,
}

impl GfaSegmentRef<'_> {
    /// Return segment length in bp, reading `LN:i:` from tab-separated tags when sequence is `*`.
    #[must_use]
    pub fn seq_len(&self) -> usize {
        if self.sequence != b"*" {
            return self.sequence.len();
        }
        self.tags
            .into_iter()
            .flat_map(|t| t.split('\t'))
            .find_map(|tag| tag.strip_prefix("LN:i:")?.parse::<usize>().ok())
            .unwrap_or(0)
    }
}

/// Zero-copy borrowed view of a GFA Link (`L` line).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GfaLinkRef<'a> {
    pub from: &'a str,
    pub from_orient: u8,
    pub to: &'a str,
    pub to_orient: u8,
    pub overlap: &'a str,
}

/// Zero-copy borrowed view of a GFA v1.0 Path (`P` line).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GfaPathRef<'a> {
    pub name: &'a str,
    pub segment_names: &'a str,
    pub overlaps: &'a str,
}

/// Zero-copy borrowed view of a GFA v1.1 Walk (`W` line).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GfaWalkRef<'a> {
    pub sample: &'a str,
    pub haplotype: &'a str,
    pub seq_id: &'a str,
    pub start: usize,
    pub end: usize,
    pub walk: &'a str,
    pub tags: Option<&'a str>,
}

/// Zero-copy borrowed view of a single GFA v1.0/v1.1 record, valid for the duration of the
/// streaming callback without heap allocations per line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GfaRecordRef<'a> {
    Segment(GfaSegmentRef<'a>),
    Link(GfaLinkRef<'a>),
    Path(GfaPathRef<'a>),
    Walk(GfaWalkRef<'a>),
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
    let reader = open_text_reader(path.as_ref())?;
    stream_gfa_records_reader(reader, on_record)
}

fn req_gfa_field<'a>(
    tok: Option<&'a str>,
    rec_type: &str,
    field: &str,
) -> std::io::Result<&'a str> {
    tok.ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("truncated GFA {rec_type} record: missing {field}"),
        )
    })
}

fn parse_walk_record_ref<'a>(
    mut parts: impl Iterator<Item = &'a str>,
) -> std::io::Result<GfaWalkRef<'a>> {
    let parse_coord = |tok: Option<&str>, field: &str| -> std::io::Result<usize> {
        match tok {
            Some("*") => Ok(0),
            None => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("truncated GFA W record: missing {field} coordinate"),
            )),
            Some(s) => s.parse::<usize>().map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("invalid GFA walk {field} coordinate: `{s}`"),
                )
            }),
        }
    };
    let sample = req_gfa_field(parts.next(), "W", "sample")?;
    let haplotype = req_gfa_field(parts.next(), "W", "haplotype")?;
    let seq_id = req_gfa_field(parts.next(), "W", "seq_id")?;
    let start = parse_coord(parts.next(), "start")?;
    let end = parse_coord(parts.next(), "end")?;
    let walk = req_gfa_field(parts.next(), "W", "walk")?;
    let tags = parts.next().filter(|rest| !rest.is_empty());
    Ok(GfaWalkRef {
        sample,
        haplotype,
        seq_id,
        start,
        end,
        walk,
        tags,
    })
}

fn parse_link_orient(tok: Option<&str>) -> std::io::Result<u8> {
    match tok {
        None | Some("" | "+") => Ok(b'+'),
        Some("-") => Ok(b'-'),
        Some(other) => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("invalid GFA link orientation: `{other}`"),
        )),
    }
}

/// Zero-copy streaming GFA reader from a file path (`.gfa`, `.gfa.gz`, `.gfa.zst`, or `"-"`).
///
/// Borrows all record fields directly from a single reused line buffer, avoiding per-record heap
/// allocations on multi-gigabyte pangenome graphs.
///
/// # Errors
/// Returns an error if the file cannot be opened, reading fails, a line exceeds 256 MiB, or the callback errors.
pub fn stream_gfa_record_refs<P, F>(path: P, on_record: F) -> Result<()>
where
    P: AsRef<Path>,
    F: FnMut(GfaRecordRef<'_>) -> Result<std::ops::ControlFlow<()>>,
{
    let reader = open_text_reader(path.as_ref())?;
    stream_gfa_record_refs_reader(reader, on_record)
}

/// Zero-copy streaming GFA reader from any [`BufRead`].
///
/// # Errors
/// Returns an error if reading fails, a line exceeds 256 MiB, or the callback errors.
pub fn stream_gfa_record_refs_reader<R, F>(mut reader: R, mut on_record: F) -> Result<()>
where
    R: BufRead,
    F: FnMut(GfaRecordRef<'_>) -> Result<std::ops::ControlFlow<()>>,
{
    let mut line = String::new();
    loop {
        line.clear();
        let n =
            std::io::Read::take(&mut reader, MAX_GFA_LINE_BYTES as u64 + 1).read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.len() > MAX_GFA_LINE_BYTES {
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
                let name = req_gfa_field(s_parts.next(), "S", "name")?;
                let raw_seq = req_gfa_field(s_parts.next(), "S", "sequence")?;
                let raw_tags = s_parts.next().filter(|t| !t.is_empty());
                on_record(GfaRecordRef::Segment(GfaSegmentRef {
                    name,
                    sequence: raw_seq.as_bytes(),
                    tags: raw_tags,
                }))?
            }
            "L" => {
                let from = req_gfa_field(parts.next(), "L", "from")?;
                let from_orient =
                    parse_link_orient(Some(req_gfa_field(parts.next(), "L", "from_orient")?))?;
                let to = req_gfa_field(parts.next(), "L", "to")?;
                let to_orient =
                    parse_link_orient(Some(req_gfa_field(parts.next(), "L", "to_orient")?))?;
                let overlap = req_gfa_field(parts.next(), "L", "overlap")?;
                on_record(GfaRecordRef::Link(GfaLinkRef {
                    from,
                    from_orient,
                    to,
                    to_orient,
                    overlap,
                }))?
            }
            "P" => on_record(GfaRecordRef::Path(GfaPathRef {
                name: req_gfa_field(parts.next(), "P", "name")?,
                segment_names: req_gfa_field(parts.next(), "P", "segment_names")?,
                overlaps: req_gfa_field(parts.next(), "P", "overlaps")?,
            }))?,
            "W" => {
                let mut w_parts = trimmed.splitn(8, '\t');
                let _ = w_parts.next();
                on_record(GfaRecordRef::Walk(parse_walk_record_ref(w_parts)?))?
            }
            "J" | "C" => std::ops::ControlFlow::Continue(()),
            other => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("unrecognized GFA record type `{other}`"),
                )
                .into());
            }
        };
        if flow.is_break() {
            break;
        }
    }

    Ok(())
}

/// Stream GFA records from a buffered reader with a single callback that may stop early.
///
/// # Errors
/// Returns an error if reading fails, a line exceeds 256 MiB, or the callback errors.
pub fn stream_gfa_records_reader<R, F>(reader: R, mut on_record: F) -> Result<()>
where
    R: BufRead,
    F: FnMut(GfaRecord) -> Result<std::ops::ControlFlow<()>>,
{
    stream_gfa_record_refs_reader(reader, |rec| {
        let owned = match rec {
            GfaRecordRef::Segment(s) => GfaRecord::Segment(GfaSegment {
                name: s.name.to_string(),
                sequence: s.sequence.to_vec(),
                tags: s.tags.map(str::to_string),
            }),
            GfaRecordRef::Link(l) => GfaRecord::Link(GfaLink {
                from: l.from.to_string(),
                from_orient: l.from_orient,
                to: l.to.to_string(),
                to_orient: l.to_orient,
                overlap: l.overlap.to_string(),
            }),
            GfaRecordRef::Path(p) => GfaRecord::Path(GfaPath {
                name: p.name.to_string(),
                segment_names: p.segment_names.to_string(),
                overlaps: p.overlaps.to_string(),
            }),
            GfaRecordRef::Walk(w) => GfaRecord::Walk(GfaWalk {
                sample: w.sample.to_string(),
                haplotype: w.haplotype.to_string(),
                seq_id: w.seq_id.to_string(),
                start: w.start,
                end: w.end,
                walk: w.walk.to_string(),
                tags: w.tags.map(str::to_string),
            }),
        };
        on_record(owned)
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::fs::File;
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

        stream_gfa_records_reader(Cursor::new(gfa.as_bytes()), |rec| {
            match rec {
                GfaRecord::Segment(s) => segs.push(s),
                GfaRecord::Link(l) => links.push(l),
                GfaRecord::Path(p) => paths.push(p),
                GfaRecord::Walk(w) => walks.push(w),
            }
            Ok(std::ops::ControlFlow::Continue(()))
        })
        .expect("stream_gfa_records_reader should succeed");

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
        stream_gfa_record_refs(&gz_path, |rec| {
            match rec {
                GfaRecordRef::Segment(s) => {
                    assert_eq!(s.sequence, b"GATTACA");
                    seg_count += 1;
                }
                GfaRecordRef::Walk(w) => {
                    assert_eq!(w.sample, "CHM13");
                    walk_count += 1;
                }
                _ => {}
            }
            Ok(std::ops::ControlFlow::Continue(()))
        })
        .unwrap();
        let _ = std::fs::remove_file(&gz_path);
        assert_eq!(seg_count, 1);
        assert_eq!(walk_count, 1);
    }

    #[test]
    fn seq_len_path_steps_break_and_invalid_walk_coord() {
        let seg = GfaSegmentRef {
            name: "s1",
            sequence: b"*",
            tags: Some("LN:i:42"),
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
