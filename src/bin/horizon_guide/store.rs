//! Adaptive in-memory vs. streaming species index store for `horizon_guide`.
//!
//! Small-to-medium multi-species panels (`<= 512 MiB` estimated input FASTA) keep all
//! `ContigIndex` instances resident in memory for zero re-reading across sister and ancestral
//! rescue. Multi-gigabase panels (e.g., 26 nematode genomes = `3.22 Gbp`, 27 insect genomes, or
//! 100 vertebrate genomes) automatically switch to streaming mode, keeping only the hinge
//! `ContigIndex` permanently resident and loading/dropping query species on demand so peak RSS
//! scales as `O(max_species)` rather than `O(sum_species)`.

use crate::{GuideCliArgs, IndexedSpecies, SpeciesContigs, select_central_hinge};
use anyhow::{Context, Result};
use horizon_aligner::{
    ProgressLogger, effective_contig_count, merge_contig_sketches, prepare_contig_index_owned,
    read_fasta_contigs,
};
use rayon::prelude::*;
use std::fs;
use std::path::{Path, PathBuf};

/// Estimated uncompressed FASTA byte budget above which `horizon_guide` streams query species
/// on demand instead of holding every species' `ContigIndex` (~14 B/bp) in RAM simultaneously.
const STREAMING_THRESHOLD_BYTES: u64 = 512 * 1024 * 1024;
pub(crate) const HINGE_SKETCH_K: usize = 512;

pub(crate) enum SpeciesStore {
    InMemory(Vec<IndexedSpecies>),
    Streaming {
        specs: Vec<(String, PathBuf)>,
        max_contigs: usize,
        respect_softmask: bool,
        hinge_idx: usize,
        hinge: IndexedSpecies,
    },
}

fn estimate_uncompressed_bytes(specs: &[(String, PathBuf)]) -> u64 {
    specs
        .iter()
        .map(|(_, path)| {
            let len = fs::metadata(path).map_or(0, |m| m.len());
            let is_gz = path
                .extension()
                .and_then(|e| e.to_str())
                .is_some_and(|e| e.eq_ignore_ascii_case("gz") || e.eq_ignore_ascii_case("bgz"));
            if is_gz { len.saturating_mul(4) } else { len }
        })
        .sum()
}

pub(crate) fn load_indexed_species(
    sp: &str,
    path: &Path,
    max_contigs: usize,
    respect_softmask: bool,
) -> Result<IndexedSpecies> {
    let contigs = read_fasta_contigs(path, max_contigs)?;
    let idx = contigs
        .into_par_iter()
        .map(|(n, s)| Ok((n, prepare_contig_index_owned(s, respect_softmask)?)))
        .collect::<Result<Vec<_>>>()?;
    Ok((sp.to_owned(), idx))
}

impl SpeciesStore {
    /// Open inputs in either `InMemory` or `Streaming` mode, returning `(store, hinge_idx, sketches)`.
    ///
    /// In `Streaming` mode with an explicit `--hinge`, only the hinge is indexed upfront; query
    /// species sketches are populated during the direct-hinge pass at zero extra I/O cost.
    #[allow(clippy::too_many_lines)]
    pub(crate) fn open(
        cli: &GuideCliArgs,
        logger: ProgressLogger,
    ) -> Result<(Self, usize, Vec<Vec<u64>>)> {
        let est_bytes = estimate_uncompressed_bytes(&cli.specs);
        let use_progressive = cli.progressive && cli.specs.len() >= 3;
        if est_bytes <= STREAMING_THRESHOLD_BYTES {
            let indexed: Vec<IndexedSpecies> = cli
                .specs
                .par_iter()
                .map(|(sp, path)| {
                    load_indexed_species(sp, path, cli.max_contigs, cli.respect_softmask)
                })
                .collect::<Result<Vec<_>>>()?;
            let total_contigs: usize = indexed.iter().map(|(_, idx)| idx.len()).sum();
            let total_bp: usize = indexed
                .iter()
                .flat_map(|(_, idx)| idx.iter())
                .map(|(_, ix)| ix.seq().len())
                .sum();
            #[allow(clippy::cast_precision_loss)]
            let total_mb = (total_bp as f64) / 1_000_000.0;
            logger.log(format_args!(
                "indexed {} species ({total_contigs} contigs, {total_mb:.2} Mb)",
                indexed.len()
            ));
            let sketches = if cli.hinge.is_none() || use_progressive {
                indexed
                    .par_iter()
                    .map(|(_, idx_list)| merge_contig_sketches(idx_list, HINGE_SKETCH_K))
                    .collect()
            } else {
                Vec::new()
            };
            let hinge_idx = if let Some(ref hinge_name) = cli.hinge {
                indexed
                    .iter()
                    .position(|(sp, _)| sp == hinge_name)
                    .with_context(|| {
                        format!("--hinge species `{hinge_name}` not found in input specs")
                    })?
            } else {
                let eff_contigs: Vec<f64> = indexed
                    .iter()
                    .map(|(_, idx)| {
                        effective_contig_count(idx.iter().map(|(_, ix)| ix.seq().len()))
                    })
                    .collect();
                select_central_hinge(&sketches, &eff_contigs)
            };
            return Ok((Self::InMemory(indexed), hinge_idx, sketches));
        }

        #[allow(clippy::cast_precision_loss)]
        let est_mb = (est_bytes as f64) / (1024.0 * 1024.0);
        logger.log(format_args!(
            "streaming mode enabled ({est_mb:.1} MiB > 512 MiB threshold): indexing query species on demand"
        ));

        let (hinge_idx, sketches, hinge) = if let Some(ref hinge_name) = cli.hinge {
            let h_idx = cli
                .specs
                .iter()
                .position(|(sp, _)| sp == hinge_name)
                .with_context(|| {
                    format!("--hinge species `{hinge_name}` not found in input specs")
                })?;
            let (ref h_sp, ref h_path) = cli.specs[h_idx];
            let hinge = load_indexed_species(h_sp, h_path, cli.max_contigs, cli.respect_softmask)?;
            let mut sketches = vec![Vec::new(); cli.specs.len()];
            if use_progressive {
                sketches[h_idx] = merge_contig_sketches(&hinge.1, HINGE_SKETCH_K);
            }
            (h_idx, sketches, hinge)
        } else {
            let mut sketches = Vec::with_capacity(cli.specs.len());
            let mut eff_contigs = Vec::with_capacity(cli.specs.len());
            for (sp, path) in &cli.specs {
                let loaded = load_indexed_species(sp, path, cli.max_contigs, cli.respect_softmask)?;
                sketches.push(merge_contig_sketches(&loaded.1, HINGE_SKETCH_K));
                eff_contigs.push(effective_contig_count(
                    loaded.1.iter().map(|(_, ix)| ix.seq().len()),
                ));
            }
            let h_idx = select_central_hinge(&sketches, &eff_contigs);
            let (ref h_sp, ref h_path) = cli.specs[h_idx];
            let hinge = load_indexed_species(h_sp, h_path, cli.max_contigs, cli.respect_softmask)?;
            (h_idx, sketches, hinge)
        };

        let h_bp: usize = hinge.1.iter().map(|(_, ix)| ix.seq().len()).sum();
        #[allow(clippy::cast_precision_loss)]
        let h_mb = (h_bp as f64) / 1_000_000.0;
        logger.log(format_args!(
            "indexed hinge `{}` ({} contigs, {h_mb:.2} Mb); {} query species will stream",
            hinge.0,
            hinge.1.len(),
            cli.specs.len().saturating_sub(1),
        ));

        Ok((
            Self::Streaming {
                specs: cli.specs.clone(),
                max_contigs: cli.max_contigs,
                respect_softmask: cli.respect_softmask,
                hinge_idx,
                hinge,
            },
            hinge_idx,
            sketches,
        ))
    }

    pub(crate) const fn len(&self) -> usize {
        match self {
            Self::InMemory(v) => v.len(),
            Self::Streaming { specs, .. } => specs.len(),
        }
    }

    pub(crate) fn in_memory(&self) -> Option<&[IndexedSpecies]> {
        match self {
            Self::InMemory(v) => Some(v),
            Self::Streaming { .. } => None,
        }
    }

    pub(crate) fn species_name(&self, idx: usize) -> &str {
        match self {
            Self::InMemory(v) => &v[idx].0,
            Self::Streaming { specs, .. } => &specs[idx].0,
        }
    }

    pub(crate) fn hinge(&self, hinge_idx: usize) -> &IndexedSpecies {
        match self {
            Self::InMemory(v) => &v[hinge_idx],
            Self::Streaming { hinge, .. } => hinge,
        }
    }

    /// Load and index species `idx` on demand in `Streaming` mode (uses all Rayon threads within
    /// the species, then lets the caller drop the returned `IndexedSpecies` immediately).
    pub(crate) fn load_streaming_species(&self, idx: usize) -> Result<IndexedSpecies> {
        match self {
            Self::Streaming {
                specs,
                max_contigs,
                respect_softmask,
                ..
            } => {
                let (ref sp, ref path) = specs[idx];
                load_indexed_species(sp, path, *max_contigs, *respect_softmask)
            }
            Self::InMemory(_) => unreachable!("load_streaming_species called on InMemory store"),
        }
    }

    /// Drop heavy `ContigIndex` tables and produce raw DNA contigs for downstream MAF/GFA/HAL
    /// export. In `Streaming` mode, raw query FASTAs are only reloaded when `need_sequences` is
    /// true (i.e., `--hal` export), keeping peak RAM bounded to raw DNA bytes with zero indices.
    pub(crate) fn into_species_contigs(self, need_sequences: bool) -> Result<Vec<SpeciesContigs>> {
        match self {
            Self::InMemory(v) => Ok(v
                .into_iter()
                .map(|(sp, idxs)| {
                    let contigs = idxs
                        .into_iter()
                        .map(|(name, ix)| (name, ix.into_seq()))
                        .collect();
                    (sp, contigs)
                })
                .collect()),
            Self::Streaming {
                specs,
                max_contigs,
                hinge_idx,
                hinge,
                ..
            } => {
                let hinge_contigs: Vec<(String, Vec<u8>)> = hinge
                    .1
                    .into_iter()
                    .map(|(name, ix)| (name, ix.into_seq()))
                    .collect();
                let mut out: Vec<SpeciesContigs> = Vec::with_capacity(specs.len());
                for (i, (sp, path)) in specs.into_iter().enumerate() {
                    if i == hinge_idx {
                        out.push((sp, Vec::new()));
                    } else if need_sequences {
                        let contigs = read_fasta_contigs(&path, max_contigs)?;
                        out.push((sp, contigs));
                    } else {
                        out.push((sp, Vec::new()));
                    }
                }
                out[hinge_idx].1 = hinge_contigs;
                Ok(out)
            }
        }
    }
}
