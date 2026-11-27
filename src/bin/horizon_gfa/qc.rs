//! Streaming GFA v1.0/v1.1 quality-control audit (`horizon_gfa qc`):
//! dialect classification (`W_only` / `P_only` / `MIXED_PW`), per-`(sample, hap)` retained bp
//! and placeholder-step accounting, and budgeted Union-Find connected-component statistics.

use anyhow::{Context, Result, bail};
use horizon_aligner::gfa::{
    GfaRecordRef, parse_path_steps, parse_walk_steps, stream_gfa_record_refs,
};
use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

fn dsu_find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        let gp = parent[parent[x as usize] as usize];
        parent[x as usize] = gp;
        x = gp;
    }
    x
}

fn dsu_union(parent: &mut [u32], size: &mut [u32], a: u32, b: u32) {
    let ra = dsu_find(parent, a);
    let rb = dsu_find(parent, b);
    if ra == rb {
        return;
    }
    if size[ra as usize] < size[rb as usize] {
        parent[ra as usize] = rb;
        size[rb as usize] += size[ra as usize];
    } else {
        parent[rb as usize] = ra;
        size[ra as usize] += size[rb as usize];
    }
}

struct QcGraphIndex {
    max_segments: usize,
    over_budget: bool,
    seg_to_id: HashMap<String, u32>,
    seg_lens: Vec<u64>,
    seg_is_star: Vec<bool>,
    dsu_parent: Vec<u32>,
    dsu_size: Vec<u32>,
}

impl QcGraphIndex {
    fn new(max_segments: usize) -> Self {
        Self {
            max_segments,
            over_budget: false,
            seg_to_id: HashMap::new(),
            seg_lens: Vec::new(),
            seg_is_star: Vec::new(),
            dsu_parent: Vec::new(),
            dsu_size: Vec::new(),
        }
    }

    fn intern(&mut self, name: &str, len: u64, is_star: bool) -> Option<u32> {
        if let Some(&id) = self.seg_to_id.get(name) {
            if len > 0 && self.seg_lens[id as usize] == 0 {
                self.seg_lens[id as usize] = len;
            }
            if is_star {
                self.seg_is_star[id as usize] = true;
            }
            return Some(id);
        }
        if self.over_budget || self.seg_lens.len() >= self.max_segments {
            if !self.over_budget {
                self.over_budget = true;
                self.seg_to_id.clear();
                self.seg_to_id.shrink_to_fit();
                self.seg_lens.clear();
                self.seg_lens.shrink_to_fit();
                self.seg_is_star.clear();
                self.seg_is_star.shrink_to_fit();
                self.dsu_parent.clear();
                self.dsu_parent.shrink_to_fit();
                self.dsu_size.clear();
                self.dsu_size.shrink_to_fit();
            }
            return None;
        }
        #[allow(clippy::cast_possible_truncation)]
        let id = self.seg_lens.len() as u32;
        self.seg_to_id.insert(name.to_owned(), id);
        self.seg_lens.push(len);
        self.seg_is_star.push(is_star);
        self.dsu_parent.push(id);
        self.dsu_size.push(1);
        Some(id)
    }

    fn union(&mut self, a: u32, b: u32) {
        dsu_union(&mut self.dsu_parent, &mut self.dsu_size, a, b);
    }
}

#[derive(Default)]
struct QcHapStats {
    records: u64,
    retained_bp: u64,
    placeholder_steps: u64,
    bp_unknown: bool,
    stars_unknown: bool,
}

#[allow(clippy::too_many_lines)]
pub(crate) fn cmd_qc(args: &[String]) -> Result<()> {
    let t0 = Instant::now();
    let mut gfa_path: Option<String> = None;
    let mut max_segments: usize = 5_000_000;
    let mut i = 2usize;
    while i < args.len() {
        match args[i].as_str() {
            "--max-segments" => {
                max_segments = args
                    .get(i + 1)
                    .context("missing value for --max-segments")?
                    .parse()
                    .context("parsing --max-segments")?;
                anyhow::ensure!(
                    u32::try_from(max_segments).is_ok(),
                    "--max-segments must be <= {} (got {max_segments})",
                    u32::MAX
                );
                i += 2;
            }
            other if other.starts_with('-') && other != "-" => {
                bail!("unknown flag in qc: {other}");
            }
            other => {
                anyhow::ensure!(
                    gfa_path.is_none(),
                    "qc accepts only one <file.gfa[.gz]|-> argument, got extra `{other}`"
                );
                gfa_path = Some(other.to_owned());
                i += 1;
            }
        }
    }
    let gfa_path = gfa_path.context("qc requires <file.gfa[.gz]|->")?;

    let mut seg_count: u64 = 0;
    let mut seg_bp: u64 = 0;
    let mut placeholder_segs: u64 = 0;
    let mut link_count: u64 = 0;
    let mut path_count: u64 = 0;
    let mut walk_count: u64 = 0;
    let mut unresolved_steps: u64 = 0;
    let mut graph_idx = QcGraphIndex::new(max_segments);
    let mut hap_stats: BTreeMap<(String, String), QcHapStats> = BTreeMap::new();

    stream_gfa_record_refs(&gfa_path, |rec| {
        match rec {
            GfaRecordRef::Segment(s) => {
                let len = s.seq_len() as u64;
                let is_star = s.sequence == b"*";
                seg_count += 1;
                seg_bp = seg_bp.saturating_add(len);
                if is_star {
                    placeholder_segs += 1;
                }
                let _ = graph_idx.intern(s.name, len, is_star);
            }
            GfaRecordRef::Link(l) => {
                link_count += 1;
                if !graph_idx.over_budget
                    && let (Some(u), Some(v)) = (
                        graph_idx.intern(l.from, 0, false),
                        graph_idx.intern(l.to, 0, false),
                    )
                {
                    graph_idx.union(u, v);
                }
            }
            GfaRecordRef::Path(p) => {
                path_count += 1;
                let (sample, hap, _) = horizon_aligner::gfa::parse_path_sample_fields(p.name);
                let entry = hap_stats
                    .entry((sample.to_owned(), hap.to_owned()))
                    .or_default();
                entry.records += 1;
                if graph_idx.over_budget {
                    entry.bp_unknown = true;
                    entry.stars_unknown = true;
                } else {
                    for (seg, _) in parse_path_steps(p.segment_names) {
                        if let Some(&id) = graph_idx.seg_to_id.get(seg) {
                            entry.retained_bp = entry
                                .retained_bp
                                .saturating_add(graph_idx.seg_lens[id as usize]);
                            if graph_idx.seg_is_star[id as usize] {
                                entry.placeholder_steps += 1;
                            }
                        } else {
                            unresolved_steps += 1;
                        }
                    }
                }
            }
            GfaRecordRef::Walk(w) => {
                walk_count += 1;
                let span_bp = w.end.saturating_sub(w.start) as u64;
                let entry = hap_stats
                    .entry((w.sample.to_owned(), w.haplotype.to_owned()))
                    .or_default();
                entry.records += 1;
                if graph_idx.over_budget {
                    if span_bp > 0 {
                        entry.retained_bp = entry.retained_bp.saturating_add(span_bp);
                    } else {
                        entry.bp_unknown = true;
                    }
                    entry.stars_unknown = true;
                } else {
                    let mut sum_bp = 0u64;
                    for (seg, _) in parse_walk_steps(w.walk) {
                        if let Some(&id) = graph_idx.seg_to_id.get(seg) {
                            sum_bp = sum_bp.saturating_add(graph_idx.seg_lens[id as usize]);
                            if graph_idx.seg_is_star[id as usize] {
                                entry.placeholder_steps += 1;
                            }
                        } else {
                            unresolved_steps += 1;
                        }
                    }
                    entry.retained_bp = entry.retained_bp.saturating_add(if span_bp > 0 {
                        span_bp
                    } else {
                        sum_bp
                    });
                }
            }
        }
        Ok(std::ops::ControlFlow::Continue(()))
    })?;

    let dialect = match (path_count > 0, walk_count > 0) {
        (true, true) => "MIXED_PW",
        (true, false) => "P_only",
        (false, true) => "W_only",
        (false, false) => "none",
    };
    let (comp_str, largest_str) = if graph_idx.over_budget {
        (
            "skipped_over_budget".to_owned(),
            "skipped_over_budget".to_owned(),
        )
    } else {
        let mut n_comp = 0usize;
        let mut largest = 0u32;
        for id in 0..graph_idx.dsu_parent.len() {
            #[allow(clippy::cast_possible_truncation)]
            if graph_idx.dsu_parent[id] == id as u32 {
                n_comp += 1;
                largest = largest.max(graph_idx.dsu_size[id]);
            }
        }
        (n_comp.to_string(), largest.to_string())
    };
    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    println!(
        "gfa={gfa_path} dialect={dialect} segments={seg_count} segment_bp={seg_bp} placeholder_segments={placeholder_segs} links={link_count} paths={path_count} walks={walk_count} components={comp_str} largest_component_segs={largest_str} haplotypes={} unresolved_steps={unresolved_steps} elapsed={:.2}s peak_rss_mb={hwm_mb:.1}",
        hap_stats.len(),
        t0.elapsed().as_secs_f64()
    );
    println!("sample\thaplotype\trecords\tretained_bp\tplaceholder_steps");
    for ((sample, hap), st) in hap_stats {
        let bp_cell = if st.bp_unknown {
            "skipped_over_budget".to_owned()
        } else {
            st.retained_bp.to_string()
        };
        let star_cell = if st.stars_unknown {
            "skipped_over_budget".to_owned()
        } else {
            st.placeholder_steps.to_string()
        };
        println!("{sample}\t{hap}\t{}\t{bp_cell}\t{star_cell}", st.records);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qc_index_handles_l_before_s_star_and_over_budget_cleanup() {
        let mut idx = QcGraphIndex::new(3);
        let u = idx.intern("s1", 0, false).expect("s1");
        let v = idx.intern("s2", 0, false).expect("s2");
        idx.union(u, v);
        // L preceded unsized `*` segment `s2`: updating s2 must still set seg_is_star = true.
        let s2 = idx.intern("s2", 0, true).expect("s2 update");
        assert_eq!(s2, v);
        assert!(idx.seg_is_star[v as usize]);
        assert!(idx.intern("s3", 10, false).is_some());
        // 4th segment trips max_segments=3 and frees all vectors including seg_lens/seg_is_star.
        assert!(idx.intern("s4", 10, false).is_none());
        assert!(idx.over_budget);
        assert_eq!(idx.seg_lens.len(), 0);
        assert_eq!(idx.seg_is_star.len(), 0);
        assert_eq!(idx.dsu_parent.len(), 0);
    }
}
