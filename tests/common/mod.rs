//! Shared deterministic DNA generators and RAII temporary-directory helper for E2E CLI tests.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_DIR_SEQ: AtomicU64 = AtomicU64::new(0);

/// RAII temporary directory that automatically cleans up on drop.
pub struct TestTempDir {
    path: PathBuf,
}

impl TestTempDir {
    /// Create a fresh unique temporary directory prefixed with `prefix`.
    pub fn new(prefix: &str) -> Self {
        let seq = TEMP_DIR_SEQ.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!("{prefix}_{}_{seq}", std::process::id()));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).expect("create test temp dir");
        Self { path }
    }

    /// Borrow the directory path.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Join a relative file or subdirectory name onto this temp directory.
    pub fn join(&self, rel: impl AsRef<Path>) -> PathBuf {
        self.path.join(rel)
    }
}

impl Drop for TestTempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

/// Deterministic pseudo-random ACGT sequence generator for integration tests.
pub fn synthetic_dna(seed: u64, len: usize) -> Vec<u8> {
    const BASES: [u8; 4] = *b"ACGT";
    let mut state = seed;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        out.push(BASES[((state >> 33) & 3) as usize]);
    }
    out
}

/// Mutate `src` deterministically at roughly `rate_per_100` percent of positions.
pub fn mutate_dna(src: &[u8], seed: u64, rate_per_100: u64) -> Vec<u8> {
    const BASES: [u8; 4] = *b"ACGT";
    let mut state = seed;
    let mut out = src.to_vec();
    for b in &mut out {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        if (state >> 33) % 100 < rate_per_100 {
            *b = BASES[((state >> 45) & 3) as usize];
        }
    }
    out
}
