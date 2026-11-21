//! Streaming compressed I/O helpers supporting plain, `.gz`/`.bgz`, `.zst`/`.zstd`, and `stdin`/`stdout`.

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
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
        let _ = self.child.wait();
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
    if has_ext(path, &["zst", "zstd"]) {
        let mut child = Command::new("zstd")
            .arg("-dc")
            .arg("--")
            .arg(path)
            .stdin(Stdio::null())
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
    let file = File::open(path).with_context(|| format!("opening '{}'", path.display()))?;
    if has_ext(path, &["gz", "bgz"]) {
        Ok(Box::new(BufReader::with_capacity(
            256 * 1024,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(256 * 1024, file)))
    }
}

/// Open a buffered text writer for `out_path` (plain, `.gz`/`.bgz`, `.zst`/`.zstd`), or `stdout`.
///
/// # Errors
/// Returns an error if `out_path` cannot be created or `zstd` cannot be spawned for `.zst` outputs.
pub fn open_text_writer(out_path: Option<&str>) -> Result<BufWriter<Box<dyn Write>>> {
    let Some(raw) = out_path.filter(|&p| p != "-") else {
        return Ok(BufWriter::new(Box::new(io::stdout().lock())));
    };
    let path = Path::new(raw);
    let file = File::create(path).with_context(|| format!("creating '{}'", path.display()))?;
    if has_ext(path, &["zst", "zstd"]) {
        let mut child = match Command::new("zstd")
            .arg("-c")
            .arg("-T0")
            .stdin(Stdio::piped())
            .stdout(Stdio::from(file))
            .stderr(Stdio::inherit())
            .spawn()
        {
            Ok(c) => c,
            Err(e) => {
                let _ = std::fs::remove_file(path);
                return Err(e).with_context(|| format!("spawning `zstd -c` for '{}'", path.display()));
            }
        };
        let stdin = child.stdin.take().context("capturing `zstd -c` stdin pipe")?;
        return Ok(BufWriter::new(Box::new(ZstdWriteChild {
            child,
            stdin: Some(stdin),
        })));
    }
    if has_ext(path, &["gz", "bgz"]) {
        return Ok(BufWriter::new(Box::new(GzEncoder::new(
            file,
            Compression::default(),
        ))));
    }
    Ok(BufWriter::new(Box::new(file)))
}

#[cfg(test)]
mod tests {
    use super::*;

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
                w.flush().expect("flush");
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
            w.flush().expect("flush");
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
}
