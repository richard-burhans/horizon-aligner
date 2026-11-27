use anyhow::{Context, Result, bail};
use horizon_aligner::{
    PairSchedule, parse_species_fasta_spec, parse_threads_arg, validate_species_name,
};
use std::collections::HashSet;
use std::path::PathBuf;

const TBA_USAGE: &str = "\
usage:
  horizon_tba [OPTIONS] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>...
  horizon_tba --from-mafs <maf_dir> [OPTIONS] <out_prefix> <sp1> <sp2>...
  horizon_tba --extend-maf <backbone.maf> [OPTIONS] <out_prefix> <max_contigs> <hinge=fa>... <tip=fa>...

options:
  --tree <newick|file>     Explicit Newick guide tree (all species, or hinges in two-tier mode)
  --hinges <sp1,sp2,...>   Two-tier mode: all-vs-all on hinges + O(N-H) tip-to-hinge projection
  --auto-hinges <H>        Automatically select H clade hinges + tip map via Eisenstein sketches
  --tip-map <tsv>          Two-column <tip> <hinge> TSV for explicit two-tier / extend assignments
  --extend-maf <maf>       Post-hoc: project new tip species onto an existing multi-species TBA MAF
  --progressive-ancestors  Figure 1D Progressive Cactus-style O(N) bottom-up ancestral alignment
  --ancestral-model <m>    Ancestral model for --progressive-ancestors: fitch (default) | ml (HKY85 + 2-state CTMC indel)
  --max-outgroups <K>      Outgroup leaves per internal subproblem in --progressive-ancestors (default: 3)
  --emit-ancestors <dir>   Write reconstructed Anc*.fa genomes & ancestors_summary.tsv into <dir>
  --hal                    Write native hierarchical <out_prefix>.hal (with --progressive-ancestors)
  --variation-gfa          Emit sequence-bearing base-level variation GFA v1.1 (for resolve --vcf)
  --out-gz                 Gzip-compress final <out_prefix>.maf.gz and <out_prefix>.gfa.gz directly
  --from-mafs <dir>        Merge-only mode reading precomputed pairwise MAFs from <dir>
  --pairwise-dir <dir>     Save/resume pairwise MAFs (<spA>__<spB>.maf[.gz]) in <dir> + lazy merge
  --pairwise-gz            Gzip-compress pairwise MAFs written into --pairwise-dir (.maf.gz)
  --respect-softmask       Skip softmasked lowercase bases during seeding
  --tantan                 Softmask tandem repeats with Frith's tantan HMM before seeding
  --rbest | --top-k K | --all-pairs   Contig-pair scheduling policy (default: rbest when max_contigs!=1)
  -t, --threads <N>        Limit Rayon worker threads (default: all cores)
  -q, --quiet              Suppress timestamped stderr progress logging";

#[allow(clippy::struct_excessive_bools)]
pub(crate) struct TbaCliArgs {
    pub(crate) tree_arg: Option<String>,
    pub(crate) hinges_arg: Option<String>,
    pub(crate) auto_hinges: Option<usize>,
    pub(crate) tip_map_path: Option<PathBuf>,
    pub(crate) extend_maf_path: Option<PathBuf>,
    pub(crate) progressive_ancestors: bool,
    pub(crate) ancestral_ml: bool,
    pub(crate) max_outgroups: usize,
    pub(crate) emit_ancestors_dir: Option<PathBuf>,
    pub(crate) write_hal: bool,
    pub(crate) variation_gfa: bool,
    pub(crate) out_gz: bool,
    pub(crate) from_mafs_dir: Option<PathBuf>,
    pub(crate) pairwise_dir: Option<PathBuf>,
    pub(crate) pairwise_gz: bool,
    pub(crate) respect_softmask: bool,
    pub(crate) use_tantan: bool,
    pub(crate) quiet: bool,
    pub(crate) threads: Option<usize>,
    pub(crate) top_k: usize,
    pub(crate) out_prefix: String,
    pub(crate) max_contigs: usize,
    pub(crate) species_names: Vec<String>,
    pub(crate) fasta_specs: Vec<(String, PathBuf)>,
}

#[allow(clippy::too_many_lines)]
pub(crate) fn parse_args() -> Result<Option<TbaCliArgs>> {
    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    let mut tree_arg = None;
    let mut hinges_arg = None;
    let mut auto_hinges = None;
    let mut tip_map_path = None;
    let mut extend_maf_path = None;
    let mut progressive_ancestors = false;
    let mut ancestral_ml = false;
    let mut max_outgroups = 3usize;
    let mut emit_ancestors_dir = None;
    let mut write_hal = false;
    let mut variation_gfa = false;
    let mut out_gz = false;
    let mut from_mafs_dir = None;
    let mut pairwise_dir = None;
    let mut pairwise_gz = false;
    let mut respect_softmask = false;
    let mut use_tantan = false;
    let mut quiet = false;
    let mut threads = None;
    let mut schedule_opt: Option<PairSchedule> = None;
    let mut pos = Vec::new();
    let mut i = 0usize;
    while i < raw_args.len() {
        match raw_args[i].as_str() {
            "-h" | "--help" => {
                println!("{TBA_USAGE}");
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("horizon_tba {}", env!("CARGO_PKG_VERSION"));
                return Ok(None);
            }
            "-t" | "--threads" => {
                let raw = raw_args
                    .get(i + 1)
                    .context("--threads requires an integer argument")?;
                threads = Some(parse_threads_arg(raw)?);
                i += 2;
            }
            "--tree" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--tree requires a Newick string or file path")?;
                tree_arg = Some(val.clone());
                i += 2;
            }
            "--hinges" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--hinges requires a comma-separated list of species")?;
                hinges_arg = Some(val.clone());
                i += 2;
            }
            "--auto-hinges" => {
                let raw = raw_args
                    .get(i + 1)
                    .context("--auto-hinges requires a positive integer H >= 2")?;
                let parsed = raw
                    .parse::<usize>()
                    .with_context(|| format!("invalid --auto-hinges value `{raw}`"))?;
                anyhow::ensure!(parsed >= 2, "--auto-hinges requires H >= 2 (got {parsed})");
                auto_hinges = Some(parsed);
                i += 2;
            }
            "--tip-map" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--tip-map requires a file path")?;
                tip_map_path = Some(PathBuf::from(val));
                i += 2;
            }
            "--extend-maf" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--extend-maf requires a MAF file path")?;
                extend_maf_path = Some(PathBuf::from(val));
                i += 2;
            }
            "--progressive-ancestors" => {
                progressive_ancestors = true;
                i += 1;
            }
            "--ancestral-model" => {
                let raw = raw_args
                    .get(i + 1)
                    .context("--ancestral-model requires `fitch` or `ml`")?;
                match raw.as_str() {
                    "fitch" => ancestral_ml = false,
                    "ml" | "transducer" => ancestral_ml = true,
                    other => {
                        bail!("invalid --ancestral-model `{other}` (expected `fitch` or `ml`)")
                    }
                }
                i += 2;
            }
            "--max-outgroups" => {
                let raw = raw_args
                    .get(i + 1)
                    .context("--max-outgroups requires a non-negative integer")?;
                max_outgroups = raw
                    .parse::<usize>()
                    .with_context(|| format!("invalid --max-outgroups value `{raw}`"))?;
                i += 2;
            }
            "--emit-ancestors" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--emit-ancestors requires a directory path")?;
                emit_ancestors_dir = Some(PathBuf::from(val));
                i += 2;
            }
            "--hal" => {
                write_hal = true;
                i += 1;
            }
            "--variation-gfa" => {
                variation_gfa = true;
                i += 1;
            }
            "--out-gz" => {
                out_gz = true;
                i += 1;
            }
            "--from-mafs" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--from-mafs requires a directory path")?;
                from_mafs_dir = Some(PathBuf::from(val));
                i += 2;
            }
            "--pairwise-dir" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--pairwise-dir requires a directory path")?;
                pairwise_dir = Some(PathBuf::from(val));
                i += 2;
            }
            "--pairwise-gz" => {
                pairwise_gz = true;
                i += 1;
            }
            "--respect-softmask" => {
                respect_softmask = true;
                i += 1;
            }
            "--tantan" => {
                use_tantan = true;
                i += 1;
            }
            "-q" | "--quiet" => {
                quiet = true;
                i += 1;
            }
            "--rbest" => {
                schedule_opt = Some(PairSchedule::TopK(1));
                i += 1;
            }
            "--all-pairs" => {
                schedule_opt = Some(PairSchedule::AllPairs);
                i += 1;
            }
            "--top-k" => {
                let val = raw_args
                    .get(i + 1)
                    .context("--top-k requires an integer argument")?;
                let parsed: usize = val
                    .parse()
                    .with_context(|| format!("invalid --top-k value `{val}`"))?;
                schedule_opt = Some(if parsed == 0 {
                    PairSchedule::AllPairs
                } else {
                    PairSchedule::TopK(parsed)
                });
                i += 2;
            }
            other if other.starts_with('-') => bail!("unknown flag `{other}`"),
            other => {
                pos.push(other.to_owned());
                i += 1;
            }
        }
    }

    let (out_prefix, max_contigs, top_k, species_names, fasta_specs) = if from_mafs_dir.is_some() {
        if pos.len() < 3 {
            bail!("{TBA_USAGE}");
        }
        let mut seen = HashSet::new();
        let mut names = Vec::new();
        for sp in &pos[1..] {
            validate_species_name(sp)?;
            anyhow::ensure!(seen.insert(sp.clone()), "duplicate species name: {sp}");
            names.push(sp.clone());
        }
        (
            pos[0].clone(),
            0,
            PairSchedule::resolve_top_k(schedule_opt, 0),
            names,
            Vec::new(),
        )
    } else {
        if pos.len() < 4 {
            bail!("{TBA_USAGE}");
        }
        let max_c: usize = pos[1].parse().context("parsing max_contigs")?;
        let tk = PairSchedule::resolve_top_k(schedule_opt, max_c);
        let mut seen = HashSet::new();
        let mut names = Vec::new();
        let mut specs = Vec::new();
        for spec in &pos[2..] {
            let (sp, path) = parse_species_fasta_spec(spec)?;
            anyhow::ensure!(seen.insert(sp.clone()), "duplicate species name: {sp}");
            names.push(sp.clone());
            specs.push((sp, path));
        }
        (pos[0].clone(), max_c, tk, names, specs)
    };
    let is_two_tier = hinges_arg.is_some()
        || auto_hinges.is_some()
        || tip_map_path.is_some()
        || extend_maf_path.is_some();
    if let Some(num_h) = auto_hinges {
        anyhow::ensure!(
            hinges_arg.is_none() && tip_map_path.is_none() && extend_maf_path.is_none(),
            "--auto-hinges cannot be combined with --hinges, --tip-map, or --extend-maf"
        );
        anyhow::ensure!(
            from_mafs_dir.is_none(),
            "--auto-hinges requires FASTA inputs (cannot be combined with --from-mafs)"
        );
        anyhow::ensure!(
            num_h <= species_names.len(),
            "--auto-hinges {num_h} exceeds input species count ({})",
            species_names.len()
        );
    }
    if from_mafs_dir.is_some() {
        anyhow::ensure!(
            !use_tantan,
            "--tantan requires FASTA inputs (cannot be combined with --from-mafs)"
        );
    }
    if progressive_ancestors {
        anyhow::ensure!(
            !is_two_tier,
            "--progressive-ancestors cannot be combined with --hinges, --auto-hinges, --tip-map, or --extend-maf"
        );
        anyhow::ensure!(
            pairwise_dir.is_none(),
            "--pairwise-dir is not supported with --progressive-ancestors"
        );
    } else {
        anyhow::ensure!(
            !ancestral_ml && emit_ancestors_dir.is_none() && !write_hal,
            "--ancestral-model ml, --emit-ancestors, and --hal require --progressive-ancestors"
        );
    }
    if pairwise_gz {
        anyhow::ensure!(
            pairwise_dir.is_some(),
            "--pairwise-gz requires --pairwise-dir <dir>"
        );
    }
    if from_mafs_dir.is_some() {
        anyhow::ensure!(
            pairwise_dir.is_none(),
            "--pairwise-dir cannot be combined with --from-mafs"
        );
    }
    if write_hal {
        for sp in &species_names {
            anyhow::ensure!(
                !matches!(sp.as_str(), "Meta" | "Phylogeny" | "Genomes" | "Verison"),
                "species name `{sp}` collides with reserved HAL/HDF5 root group name"
            );
        }
    }
    if pairwise_dir.is_some() || from_mafs_dir.is_some() {
        for sp in &species_names {
            anyhow::ensure!(
                !sp.contains("__"),
                "species name `{sp}` cannot contain `__` when using --pairwise-dir or --from-mafs"
            );
        }
    }
    Ok(Some(TbaCliArgs {
        tree_arg,
        hinges_arg,
        auto_hinges,
        tip_map_path,
        extend_maf_path,
        progressive_ancestors,
        ancestral_ml,
        max_outgroups,
        emit_ancestors_dir,
        write_hal,
        variation_gfa,
        out_gz,
        from_mafs_dir,
        pairwise_dir,
        pairwise_gz,
        respect_softmask,
        use_tantan,
        quiet,
        threads,
        top_k,
        out_prefix,
        max_contigs,
        species_names,
        fasta_specs,
    }))
}
