//! Streaming compressed I/O helpers supporting plain, `.gz`/`.bgz`, `.zst`/`.zstd`, and `stdin`/`stdout`.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::{Compress, Compression, Crc, FlushCompress, Status};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

fn has_ext(path: &Path, exts: &[&str]) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| exts.iter().any(|want| e.eq_ignore_ascii_case(want)))
}

struct ZstdReadChild {
    child: Child,
    stdout: Option<ChildStdout>,
    exited: bool,
}

impl Read for ZstdReadChild {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let Some(ref mut stdout) = self.stdout else {
            return Ok(0);
        };
        let n = stdout.read(buf)?;
        if n == 0 && !self.exited {
            self.exited = true;
            self.stdout = None;
            let status = self.child.wait()?;
            if !status.success() {
                return Err(io::Error::other(format!(
                    "`zstd -dc` failed with exit status {status}"
                )));
            }
        }
        Ok(n)
    }
}

impl Drop for ZstdReadChild {
    fn drop(&mut self) {
        drop(self.stdout.take());
        if !self.exited {
            let _ = self.child.kill();
            let _ = self.child.wait();
        }
    }
}

struct ZstdWriteChild {
    child: Child,
    stdin: Option<ChildStdin>,
    exited: bool,
}

impl Write for ZstdWriteChild {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.stdin
            .as_mut()
            .ok_or_else(|| io::Error::new(io::ErrorKind::BrokenPipe, "zstd stdin closed"))?
            .write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Some(ref mut stdin) = self.stdin {
            stdin.flush()?;
        }
        Ok(())
    }
}

impl Drop for ZstdWriteChild {
    fn drop(&mut self) {
        if let Some(mut stdin) = self.stdin.take() {
            let _ = stdin.flush();
        }
        if !self.exited {
            let _ = self.child.wait();
        }
    }
}

const BGZF_MAX_BLOCK_LEN: usize = 65_280;
const BGZF_MAX_CDATA_LEN: usize = 65_536 - 26;
const BGZF_EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Pure-Rust RFC 1952 + `SAMv1` BGZF (Blocked GNU Zip Format) writer compatible with `tabix` and `htslib`.
struct BgzfWriter<W: Write> {
    inner: W,
    buf: Vec<u8>,
    cdata: Vec<u8>,
}

impl<W: Write> BgzfWriter<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            buf: Vec::with_capacity(BGZF_MAX_BLOCK_LEN),
            cdata: Vec::with_capacity(BGZF_MAX_CDATA_LEN + 64),
        }
    }

    fn write_chunk(&mut self, chunk: &[u8]) -> io::Result<()> {
        if chunk.is_empty() {
            return Ok(());
        }
        let mut comp = Compress::new(Compression::default(), false);
        self.cdata.clear();
        self.cdata
            .resize(chunk.len().saturating_add(128).max(256), 0);
        let mut total_out = 0usize;
        let mut total_in = 0usize;
        loop {
            if total_out == self.cdata.len() {
                self.cdata.resize(self.cdata.len().saturating_mul(2), 0);
            }
            let before_in = comp.total_in();
            let before_out = comp.total_out();
            let status = comp
                .compress(
                    &chunk[total_in..],
                    &mut self.cdata[total_out..],
                    FlushCompress::Finish,
                )
                .map_err(|e| io::Error::other(format!("BGZF deflate error: {e}")))?;
            total_in += usize::try_from(comp.total_in() - before_in).unwrap_or(0);
            total_out += usize::try_from(comp.total_out() - before_out).unwrap_or(0);
            match status {
                Status::StreamEnd => break,
                Status::Ok | Status::BufError => {}
            }
        }
        // Incompressible chunks exceeding the 64 KiB BGZF block limit split recursively in half.
        if total_out > BGZF_MAX_CDATA_LEN && chunk.len() > 1 {
            let mid = chunk.len() / 2;
            self.write_chunk(&chunk[..mid])?;
            return self.write_chunk(&chunk[mid..]);
        }
        let block_len = 18usize.saturating_add(total_out).saturating_add(8);
        let bsize = u16::try_from(block_len.saturating_sub(1))
            .map_err(|_| io::Error::other("BGZF block exceeds 65536 bytes"))?;
        let mut header = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C',
            0x02, 0x00, 0x00, 0x00,
        ];
        header[16..18].copy_from_slice(&bsize.to_le_bytes());

        let mut crc = Crc::new();
        crc.update(chunk);
        let mut trailer = [0u8; 8];
        trailer[0..4].copy_from_slice(&crc.sum().to_le_bytes());
        let isize_u32 = u32::try_from(chunk.len())
            .map_err(|_| io::Error::other("BGZF uncompressed chunk exceeds u32::MAX"))?;
        trailer[4..8].copy_from_slice(&isize_u32.to_le_bytes());

        self.inner.write_all(&header)?;
        self.inner.write_all(&self.cdata[..total_out])?;
        self.inner.write_all(&trailer)?;
        Ok(())
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let chunk = std::mem::take(&mut self.buf);
        let res = self.write_chunk(&chunk);
        self.buf = chunk;
        self.buf.clear();
        res
    }

    fn finish(mut self) -> io::Result<W> {
        self.flush_block()?;
        self.inner.write_all(&BGZF_EOF_BLOCK)?;
        self.inner.flush()?;
        Ok(self.inner)
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, mut src: &[u8]) -> io::Result<usize> {
        let orig_len = src.len();
        while !src.is_empty() {
            let space = BGZF_MAX_BLOCK_LEN.saturating_sub(self.buf.len());
            let take = src.len().min(space);
            self.buf.extend_from_slice(&src[..take]);
            src = &src[take..];
            if self.buf.len() >= BGZF_MAX_BLOCK_LEN {
                self.flush_block()?;
            }
        }
        Ok(orig_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_block()?;
        self.inner.flush()
    }
}

fn is_bgzf_output_path(path: &Path) -> bool {
    if has_ext(path, &["bgz"]) {
        return true;
    }
    path.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|name| {
            let lower = name.to_ascii_lowercase();
            lower.ends_with(".vcf.gz") || lower.ends_with(".bed.gz")
        })
}

enum TextSink {
    Stdout(io::StdoutLock<'static>),
    File(File),
    Gz(Box<GzEncoder<File>>),
    Bgzf(Box<BgzfWriter<File>>),
    Zstd(Box<ZstdWriteChild>),
}

impl Write for TextSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Stdout(w) => w.write(buf),
            Self::File(w) => w.write(buf),
            Self::Gz(w) => w.write(buf),
            Self::Bgzf(w) => w.write(buf),
            Self::Zstd(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Stdout(w) => w.flush(),
            Self::File(w) => w.flush(),
            Self::Gz(w) => w.flush(),
            Self::Bgzf(w) => w.flush(),
            Self::Zstd(w) => w.flush(),
        }
    }
}

/// Buffered text writer that propagates compression/subprocess finalization errors via [`Self::finish`].
pub struct TextWriter {
    inner: Option<BufWriter<TextSink>>,
}

impl TextWriter {
    fn new(sink: TextSink) -> Self {
        Self {
            inner: Some(BufWriter::with_capacity(256 * 1024, sink)),
        }
    }

    /// Flush buffered bytes and finalize any `.gz`/`.zst` encoder or child process.
    ///
    /// # Errors
    /// Returns an error if flushing, `GzEncoder::finish()`, or `zstd -c` exit status fails.
    pub fn finish(mut self) -> Result<()> {
        let Some(bw) = self.inner.take() else {
            return Ok(());
        };
        let sink = bw
            .into_inner()
            .map_err(std::io::IntoInnerError::into_error)?;
        match sink {
            TextSink::Stdout(mut w) => w.flush()?,
            TextSink::File(mut w) => w.flush()?,
            TextSink::Gz(gz) => {
                let mut f = gz.finish()?;
                f.flush()?;
            }
            TextSink::Bgzf(bgzf) => {
                let mut f = bgzf.finish()?;
                f.flush()?;
            }
            TextSink::Zstd(mut zstd) => {
                if let Some(mut stdin) = zstd.stdin.take() {
                    stdin.flush()?;
                    drop(stdin);
                }
                zstd.exited = true;
                let status = zstd.child.wait()?;
                anyhow::ensure!(
                    status.success(),
                    "`zstd -c` failed with exit status {status}"
                );
            }
        }
        Ok(())
    }
}

impl Write for TextWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.inner
            .as_mut()
            .ok_or_else(|| io::Error::new(io::ErrorKind::BrokenPipe, "writer already finished"))?
            .write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Some(ref mut bw) = self.inner {
            bw.flush()?;
        }
        Ok(())
    }
}

/// Open a buffered text reader for plain, `.gz`/`.bgz`, `.zst`/`.zstd`, or `"-"` (`stdin`).
///
/// # Errors
/// Returns an error if the file cannot be opened or `zstd` cannot be spawned for `.zst` inputs.
pub fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    if path == Path::new("-") {
        return Ok(Box::new(BufReader::with_capacity(256 * 1024, io::stdin())));
    }
    let file = File::open(path).with_context(|| format!("opening '{}'", path.display()))?;
    if has_ext(path, &["zst", "zstd"]) {
        let mut child = Command::new("zstd")
            .arg("-dc")
            .stdin(Stdio::from(file))
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .with_context(|| format!("spawning `zstd -dc` for '{}'", path.display()))?;
        let stdout = child
            .stdout
            .take()
            .context("capturing `zstd -dc` stdout pipe")?;
        return Ok(Box::new(BufReader::with_capacity(
            256 * 1024,
            ZstdReadChild {
                child,
                stdout: Some(stdout),
                exited: false,
            },
        )));
    }
    if has_ext(path, &["gz", "bgz"]) {
        Ok(Box::new(BufReader::with_capacity(
            256 * 1024,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(256 * 1024, file)))
    }
}

/// Read one newline-delimited UTF-8 line from `reader` into `buf` (clearing `buf` first),
/// refusing any line longer than `max_bytes` to guard against zip-bomb / unbounded-line OOM.
/// Returns `0` on EOF.
///
/// # Errors
/// Returns an error on I/O failure, invalid UTF-8, or if a single line exceeds `max_bytes`.
pub fn read_bounded_line<R: BufRead>(
    reader: &mut R,
    buf: &mut String,
    max_bytes: usize,
    label: &str,
) -> Result<usize> {
    buf.clear();
    let cap = u64::try_from(max_bytes)
        .unwrap_or(u64::MAX - 1)
        .saturating_add(1);
    let n = BufRead::read_line(&mut Read::take(reader, cap), buf)?;
    if buf.len() > max_bytes {
        anyhow::bail!("{label} line exceeds {max_bytes}-byte safety limit");
    }
    Ok(n)
}

/// Open a buffered text writer for `out_path` (plain, `.gz`/`.bgz`, `.zst`/`.zstd`), or `stdout`.
///
/// # Errors
/// Returns an error if `out_path` cannot be created or `zstd` cannot be spawned for `.zst` outputs.
pub fn open_text_writer(out_path: Option<&str>) -> Result<TextWriter> {
    let Some(raw) = out_path.filter(|&p| p != "-") else {
        return Ok(TextWriter::new(TextSink::Stdout(io::stdout().lock())));
    };
    let path = Path::new(raw);
    if has_ext(path, &["zst", "zstd"]) {
        // Probe zstd availability before truncating/creating `path` so a missing binary never
        // destroys an existing target file.
        let status = Command::new("zstd")
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .with_context(|| format!("spawning `zstd` for '{}'", path.display()))?;
        anyhow::ensure!(
            status.success(),
            "`zstd --version` failed with exit status {status}"
        );
        let file = File::create(path).with_context(|| format!("creating '{}'", path.display()))?;
        let mut child = Command::new("zstd")
            .arg("-c")
            .arg("-T0")
            .stdin(Stdio::piped())
            .stdout(Stdio::from(file))
            .stderr(Stdio::inherit())
            .spawn()
            .with_context(|| format!("spawning `zstd -c` for '{}'", path.display()))?;
        let stdin = child
            .stdin
            .take()
            .context("capturing `zstd -c` stdin pipe")?;
        return Ok(TextWriter::new(TextSink::Zstd(Box::new(ZstdWriteChild {
            child,
            stdin: Some(stdin),
            exited: false,
        }))));
    }
    let file = File::create(path).with_context(|| format!("creating '{}'", path.display()))?;
    if is_bgzf_output_path(path) {
        return Ok(TextWriter::new(TextSink::Bgzf(Box::new(BgzfWriter::new(
            file,
        )))));
    }
    if has_ext(path, &["gz", "bgz"]) {
        return Ok(TextWriter::new(TextSink::Gz(Box::new(GzEncoder::new(
            file,
            Compression::default(),
        )))));
    }
    Ok(TextWriter::new(TextSink::File(file)))
}

/// Timestamped `stderr` progress logger shared across multi-stage CLI binaries.
#[derive(Clone, Copy)]
pub struct ProgressLogger {
    tool: &'static str,
    t0: std::time::Instant,
    quiet: bool,
}

impl ProgressLogger {
    /// Create a new logger anchored at `t0`; when `quiet` is `true`, `log` is a no-op.
    #[must_use]
    pub const fn new(tool: &'static str, t0: std::time::Instant, quiet: bool) -> Self {
        Self { tool, t0, quiet }
    }

    /// Emit a timestamped status line to `stderr` unless `quiet` was set.
    pub fn log(&self, msg: impl std::fmt::Display) {
        if !self.quiet {
            eprintln!(
                "[{} {:>7.1}s] {msg}",
                self.tool,
                self.t0.elapsed().as_secs_f64()
            );
        }
    }
}

/// Exclusive, owner-only (`0o600` on Unix) temporary file guard that auto-deletes on drop
/// unless atomically persisted to its final destination via [`Self::persist`].
pub struct AtomicTempFile {
    path: std::path::PathBuf,
    persisted: bool,
}

static ATOMIC_TMP_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

impl AtomicTempFile {
    /// Create a unique sibling temporary file next to `target`, returning the RAII guard and
    /// a `256 KiB` buffered writer.
    ///
    /// # Errors
    /// Returns an error if the sibling temporary file cannot be created.
    pub fn create_exclusive(target: &Path) -> Result<(Self, BufWriter<File>)> {
        let pid = std::process::id();
        let file_name = target
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("horizon_out");
        let parent = target.parent().unwrap_or_else(|| Path::new("."));
        loop {
            let seq = ATOMIC_TMP_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let candidate = parent.join(format!(".{file_name}.tmp.{pid}.{seq}"));
            let mut opts = std::fs::OpenOptions::new();
            opts.write(true).create_new(true);
            #[cfg(unix)]
            {
                use std::os::unix::fs::OpenOptionsExt as _;
                opts.mode(0o600);
            }
            match opts.open(&candidate) {
                Ok(f) => {
                    return Ok((
                        Self {
                            path: candidate,
                            persisted: false,
                        },
                        BufWriter::with_capacity(256 * 1024, f),
                    ));
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
                Err(e) => {
                    return Err(e).with_context(|| {
                        format!("creating exclusive temp file `{}`", candidate.display())
                    });
                }
            }
        }
    }

    /// Path of the underlying temporary file on disk.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Atomically rename the temporary file to `target`, disarming drop-time cleanup.
    ///
    /// # Errors
    /// Returns an error if the filesystem rename fails.
    pub fn persist(mut self, target: &Path) -> Result<()> {
        std::fs::rename(&self.path, target).with_context(|| {
            format!(
                "renaming `{}` -> `{}`",
                self.path.display(),
                target.display()
            )
        })?;
        self.persisted = true;
        Ok(())
    }
}

impl Drop for AtomicTempFile {
    fn drop(&mut self) {
        if !self.persisted {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atomic_temp_file_cleans_up_on_drop_and_persists_on_rename() {
        let dir = std::env::temp_dir().join(format!("horizon_io_atomic_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("mkdir");
        let target = dir.join("sp1__sp2.maf");

        let dropped_path = {
            let (guard, mut w) = AtomicTempFile::create_exclusive(&target).expect("create");
            w.write_all(b"partial").expect("write");
            let p = guard.path().to_path_buf();
            assert!(p.exists());
            drop(w);
            drop(guard);
            p
        };
        assert!(!dropped_path.exists(), "temp file must be deleted on drop");
        assert!(!target.exists());

        let (guard, mut w) = AtomicTempFile::create_exclusive(&target).expect("create2");
        w.write_all(b"##maf version=1\n").expect("write2");
        w.flush().expect("flush");
        drop(w);
        guard.persist(&target).expect("persist");
        assert_eq!(
            std::fs::read_to_string(&target).expect("read target"),
            "##maf version=1\n"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn roundtrip_plain_and_gz_text_io() {
        let dir = std::env::temp_dir().join(format!("horizon_io_util_{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("mkdir");
        for ext in ["txt", "txt.gz"] {
            let path = dir.join(format!("sample.{ext}"));
            let path_str = path.to_string_lossy().into_owned();
            {
                let mut w = open_text_writer(Some(&path_str)).expect("open writer");
                writeln!(w, "line1\nline2").expect("write");
                w.finish().expect("finish");
            }
            let mut r = open_text_reader(&path).expect("open reader");
            let mut got = String::new();
            r.read_to_string(&mut got).expect("read");
            assert_eq!(got, "line1\nline2\n");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn zst_roundtrip_and_early_drop_do_not_deadlock() {
        let has_zstd = Command::new("zstd")
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|s| s.success());
        if !has_zstd {
            return;
        }
        let dir = std::env::temp_dir().join(format!("horizon_io_zst_{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("mkdir");
        let path = dir.join("big.txt.zst");
        let path_str = path.to_string_lossy().into_owned();
        {
            let mut w = open_text_writer(Some(&path_str)).expect("open zst writer");
            for i in 0..20_000 {
                writeln!(w, "record_{i:06}\tpayload_0123456789abcdef").expect("write");
            }
            w.finish().expect("finish");
        }
        // Read just 1 line and drop early to verify ZstdReadChild::drop does not deadlock.
        {
            let mut r = open_text_reader(&path).expect("open zst reader");
            let mut first = String::new();
            r.read_line(&mut first).expect("read_line");
            assert!(first.starts_with("record_000000\t"));
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_bounded_line_rejects_oversized_line() {
        let data = b"abcdefg\n";
        let mut cursor = io::Cursor::new(data.as_slice());
        let mut buf = String::new();
        assert!(read_bounded_line(&mut cursor, &mut buf, 4, "test").is_err());

        let mut cursor_ok = io::Cursor::new(b"abc\n".as_slice());
        let n = read_bounded_line(&mut cursor_ok, &mut buf, 4, "test").expect("fits");
        assert_eq!(n, 4);
        assert_eq!(buf, "abc\n");
    }
}
