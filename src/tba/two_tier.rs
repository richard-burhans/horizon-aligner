//! Two-tier hinge + tip projection planning (`O(H^2 + (N - H))`) and post-hoc tip extension.

use anyhow::{Result, bail, ensure};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

/// Validated two-tier partition of input species into Tier-1 hinges and Tier-2 tips.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TwoTierPlan {
    /// Tier-1 backbone hinge species (aligned all-vs-all in one-shot mode, or loaded from
    /// `--extend-maf` in post-hoc mode).
    pub hinges: Vec<String>,
    /// Explicitly mapped `(tip, hinge)` pairs from `--tip-map`.
    pub explicit_tips: Vec<(String, String)>,
    /// Remaining non-hinge tip species to be auto-assigned to their closest hinge by similarity.
    pub unassigned_tips: Vec<String>,
}

/// Parse a comma-separated `--hinges` list and validate species names.
///
/// # Errors
/// Returns an error if the list is empty, contains invalid species names, or has duplicates.
pub fn parse_hinge_list(raw: &str) -> Result<Vec<String>> {
    let mut hinges = Vec::new();
    let mut seen = HashSet::new();
    for tok in raw.split(',') {
        let name = tok.trim();
        if name.is_empty() {
            continue;
        }
        crate::validate_species_name(name)?;
        ensure!(
            seen.insert(name.to_owned()),
            "duplicate hinge species `{name}` in --hinges"
        );
        hinges.push(name.to_owned());
    }
    ensure!(!hinges.is_empty(), "--hinges requires at least one species");
    Ok(hinges)
}

/// Parse a two-column whitespace/tab-separated `--tip-map` (`tip<TAB>hinge`), ignoring blank lines,
/// `#` comments, and an optional `tip hinge` header row.
///
/// # Errors
/// Returns an error if reading fails, a row does not have two columns, a species name is invalid,
/// a tip maps to itself, or a tip appears more than once.
pub fn parse_tip_map_reader<R: BufRead>(mut reader: R) -> Result<Vec<(String, String)>> {
    const MAX_TIP_MAP_LINE_BYTES: usize = 64 * 1024;
    let mut out = Vec::new();
    let mut seen_tips = HashSet::new();
    let mut line = String::new();
    let mut line_no = 0usize;
    loop {
        let n = crate::io_util::read_bounded_line(
            &mut reader,
            &mut line,
            MAX_TIP_MAP_LINE_BYTES,
            "--tip-map",
        )?;
        if n == 0 {
            break;
        }
        line_no += 1;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split_ascii_whitespace().collect();
        if cols.len() != 2 {
            bail!(
                "--tip-map line {line_no}: expected 2 columns `<tip> <hinge>`, got {} columns",
                cols.len()
            );
        }
        if out.is_empty()
            && cols[0].eq_ignore_ascii_case("tip")
            && cols[1].eq_ignore_ascii_case("hinge")
        {
            continue;
        }
        let tip = cols[0];
        let hinge = cols[1];
        crate::validate_species_name(tip)?;
        crate::validate_species_name(hinge)?;
        ensure!(
            tip != hinge,
            "--tip-map line {line_no}: tip `{tip}` cannot map to itself"
        );
        ensure!(
            seen_tips.insert(tip.to_owned()),
            "--tip-map line {line_no}: duplicate tip `{tip}`"
        );
        out.push((tip.to_owned(), hinge.to_owned()));
    }
    ensure!(!out.is_empty(), "--tip-map file contains no tip mappings");
    Ok(out)
}

impl TwoTierPlan {
    /// Resolve and validate hinges vs. tips from CLI flags (`--hinges`, `--tip-map`, and optional
    /// backbone species from `--extend-maf`).
    ///
    /// # Errors
    /// Returns an error if hinges/tips are inconsistent, overlap, reference missing species, or
    /// fail to define at least one hinge (or at least two hinges in one-shot mode).
    pub fn resolve(
        species_names: &[String],
        hinges_opt: Option<&[String]>,
        tip_map_opt: Option<&[(String, String)]>,
        backbone_species: Option<&HashSet<String>>,
    ) -> Result<Self> {
        let cli_set: HashSet<&str> = species_names.iter().map(String::as_str).collect();
        let explicit_tips: Vec<(String, String)> = tip_map_opt.unwrap_or(&[]).to_vec();
        let tip_set: HashSet<&str> = explicit_tips.iter().map(|(t, _)| t.as_str()).collect();

        // Determine hinge set: explicit --hinges > backbone species > all non-tip CLI species.
        let hinges: Vec<String> = if let Some(hs) = hinges_opt {
            hs.to_vec()
        } else if let Some(bb) = backbone_species {
            species_names
                .iter()
                .filter(|s| bb.contains(s.as_str()) && !tip_set.contains(s.as_str()))
                .cloned()
                .collect()
        } else {
            species_names
                .iter()
                .filter(|s| !tip_set.contains(s.as_str()))
                .cloned()
                .collect()
        };

        let min_hinges = if backbone_species.is_some() { 1 } else { 2 };
        ensure!(
            hinges.len() >= min_hinges,
            "two-tier mode requires at least {min_hinges} hinge species (found {})",
            hinges.len()
        );
        let hinge_set: HashSet<&str> = hinges.iter().map(String::as_str).collect();
        for h in &hinges {
            ensure!(
                cli_set.contains(h.as_str()),
                "hinge `{h}` is not listed among input species"
            );
            if let Some(bb) = backbone_species {
                ensure!(
                    bb.contains(h.as_str()),
                    "hinge `{h}` was not found in --extend-maf backbone"
                );
            }
        }

        for (tip, hinge) in &explicit_tips {
            ensure!(
                cli_set.contains(tip.as_str()),
                "--tip-map tip `{tip}` is not listed among input species"
            );
            ensure!(
                hinge_set.contains(hinge.as_str()),
                "--tip-map assigns tip `{tip}` to `{hinge}`, which is not a hinge"
            );
            ensure!(
                !hinge_set.contains(tip.as_str()),
                "species `{tip}` cannot be both a hinge and a tip"
            );
            if let Some(bb) = backbone_species {
                ensure!(
                    !bb.contains(tip.as_str()),
                    "tip `{tip}` already exists in --extend-maf backbone"
                );
            }
        }

        let unassigned_tips: Vec<String> = species_names
            .iter()
            .filter(|s| !hinge_set.contains(s.as_str()) && !tip_set.contains(s.as_str()))
            .cloned()
            .collect();
        if let Some(bb) = backbone_species {
            for tip in &unassigned_tips {
                ensure!(
                    !bb.contains(tip.as_str()),
                    "species `{tip}` already exists in --extend-maf backbone (include it in --hinges or omit it)"
                );
            }
            ensure!(
                !explicit_tips.is_empty() || !unassigned_tips.is_empty(),
                "--extend-maf requires at least one new tip species not already in the backbone"
            );
        }

        Ok(Self {
            hinges,
            explicit_tips,
            unassigned_tips,
        })
    }

    /// Assign any `unassigned_tips` to their highest-similarity hinge (breaking ties by hinge order)
    /// and return the combined `(tip, hinge, similarity)` assignments.
    ///
    /// # Errors
    /// Returns an error if any species is missing from `sp_to_idx`.
    pub fn finalize_assignments(
        &self,
        sim: &[Vec<f64>],
        sp_to_idx: &HashMap<&str, usize>,
    ) -> Result<Vec<(String, String, f64)>> {
        let mut out = Vec::with_capacity(self.explicit_tips.len() + self.unassigned_tips.len());
        for (tip, hinge) in &self.explicit_tips {
            let ti = *sp_to_idx
                .get(tip.as_str())
                .ok_or_else(|| anyhow::anyhow!("missing index for tip `{tip}`"))?;
            let hi = *sp_to_idx
                .get(hinge.as_str())
                .ok_or_else(|| anyhow::anyhow!("missing index for hinge `{hinge}`"))?;
            out.push((tip.clone(), hinge.clone(), sim[ti][hi]));
        }
        for tip in &self.unassigned_tips {
            let ti = *sp_to_idx
                .get(tip.as_str())
                .ok_or_else(|| anyhow::anyhow!("missing index for tip `{tip}`"))?;
            let mut best_h = &self.hinges[0];
            let mut best_sim = -1.0_f64;
            for h in &self.hinges {
                let hi = *sp_to_idx
                    .get(h.as_str())
                    .ok_or_else(|| anyhow::anyhow!("missing index for hinge `{h}`"))?;
                let s = sim[ti][hi];
                if s > best_sim {
                    best_sim = s;
                    best_h = h;
                }
            }
            out.push((tip.clone(), best_h.clone(), best_sim.max(0.0)));
        }
        Ok(out)
    }

    /// Build the deduplicated pairwise species index list `(min_idx, max_idx)` for two-tier mode:
    /// `H(H-1)/2` hinge pairs (when `include_hinge_pairs == true`) plus `(N - H)` tip-to-hinge pairs.
    #[must_use]
    pub fn schedule_pairs(
        &self,
        final_tips: &[(String, String)],
        sp_to_idx: &HashMap<&str, usize>,
        include_hinge_pairs: bool,
    ) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if include_hinge_pairs {
            for a in 0..self.hinges.len() {
                let ia = sp_to_idx[self.hinges[a].as_str()];
                for b in (a + 1)..self.hinges.len() {
                    let ib = sp_to_idx[self.hinges[b].as_str()];
                    pairs.push((ia.min(ib), ia.max(ib)));
                }
            }
        }
        for (tip, hinge) in final_tips {
            let it = sp_to_idx[tip.as_str()];
            let ih = sp_to_idx[hinge.as_str()];
            pairs.push((it.min(ih), it.max(ih)));
        }
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tip_map_parser_handles_comments_header_and_rejects_self_maps() {
        let raw = "# Plasmodium tips\ntip\thinge\npf_iso1\tpf3d7\npv_iso1  pv_p01\n";
        let map = parse_tip_map_reader(raw.as_bytes()).expect("valid tip map");
        assert_eq!(
            map,
            vec![
                ("pf_iso1".to_owned(), "pf3d7".to_owned()),
                ("pv_iso1".to_owned(), "pv_p01".to_owned()),
            ]
        );
        assert!(parse_tip_map_reader("pf3d7\tpf3d7\n".as_bytes()).is_err());
        assert!(parse_tip_map_reader("a\tb\na\tc\n".as_bytes()).is_err());
    }

    #[test]
    fn two_tier_plan_reduces_pair_count_to_hinge_quad_plus_linear_tips() {
        // 3 hinges + 4 tips = 7 species -> all-vs-all = 21 pairs; two-tier = 3 + 4 = 7 pairs;
        // extend-maf = 4 pairs only!
        let species: Vec<String> = ["h1", "h2", "h3", "t1", "t2", "t3", "t4"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let hinges = parse_hinge_list("h1,h2,h3").expect("hinges");
        let plan = TwoTierPlan::resolve(&species, Some(&hinges), None, None).expect("plan");
        assert_eq!(plan.hinges.len(), 3);
        assert_eq!(plan.unassigned_tips.len(), 4);

        let sp_to_idx: HashMap<&str, usize> = species
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        let mut sim = vec![vec![0.1; 7]; 7];
        sim[3][0] = 0.9; // t1 -> h1
        sim[4][0] = 0.8; // t2 -> h1
        sim[5][1] = 0.9; // t3 -> h2
        sim[6][2] = 0.9; // t4 -> h3
        let assigned = plan.finalize_assignments(&sim, &sp_to_idx).expect("assign");
        let pairs_only: Vec<(String, String)> = assigned
            .iter()
            .map(|(t, h, _)| (t.clone(), h.clone()))
            .collect();
        assert_eq!(
            pairs_only,
            vec![
                ("t1".to_owned(), "h1".to_owned()),
                ("t2".to_owned(), "h1".to_owned()),
                ("t3".to_owned(), "h2".to_owned()),
                ("t4".to_owned(), "h3".to_owned()),
            ]
        );
        assert_eq!(plan.schedule_pairs(&pairs_only, &sp_to_idx, true).len(), 7);
        assert_eq!(plan.schedule_pairs(&pairs_only, &sp_to_idx, false).len(), 4);
    }

    #[test]
    fn parse_tip_map_rejects_oversized_line() {
        let huge = format!("{} h1\n", "a".repeat(70 * 1024));
        assert!(parse_tip_map_reader(huge.as_bytes()).is_err());
    }
}
