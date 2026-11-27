//! Eisenstein bottom-k sampler sketches and O(N * H) farthest-first auto-hinge clade selection.

use super::EisensteinHasher;
use std::collections::HashSet;

/// Bottom-`s` Eisenstein sampler sketch for arbitrary-`k` strand-symmetric Jaccard estimation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EisensteinSamplerSketch {
    pub k: usize,
    pub hashes: Vec<u64>,
}

impl EisensteinSamplerSketch {
    /// Build a bottom-`max_size` strand-symmetric Eisenstein sketch across all `contigs`.
    #[must_use]
    pub fn from_contigs(
        hasher: &EisensteinHasher,
        contigs: &[&[u8]],
        k: usize,
        max_size: usize,
    ) -> Self {
        if max_size == 0 || k == 0 {
            return Self {
                k,
                hashes: Vec::new(),
            };
        }
        let mut heap = std::collections::BinaryHeap::<u64>::with_capacity(max_size + 1);
        let mut seen = HashSet::<u64>::with_capacity(max_size * 2);
        for &seq in contigs {
            hasher.for_each_kmer_hash(seq, k, |_, h| {
                if heap.len() == max_size
                    && let Some(&top) = heap.peek()
                    && h >= top
                {
                    return;
                }
                if !seen.insert(h) {
                    return;
                }
                heap.push(h);
                if heap.len() > max_size
                    && let Some(evicted) = heap.pop()
                {
                    seen.remove(&evicted);
                }
            });
        }
        let mut vals: Vec<u64> = heap.into_iter().collect();
        vals.sort_unstable();
        Self { k, hashes: vals }
    }

    /// Unbiased Broder `MIN_s(A ∪ B)` Jaccard similarity estimate.
    #[must_use]
    pub fn jaccard(&self, other: &Self) -> f64 {
        crate::seeding::sketch_jaccard(&self.hashes, &other.hashes)
    }
}

/// `(hinge_species, tip_to_hinge_assignments)` returned by [`select_auto_hinges_eisenstein`].
pub type AutoHingePlan = (Vec<String>, Vec<(String, String)>);

/// Automatically select `num_hinges` maximally diverse clade hinges and assign every remaining
/// tip species to its nearest hinge in `O(N * H)` time using Eisenstein sketches.
///
/// Uses deterministic farthest-first (k-center) selection starting from the first species
/// (or the most central medoid), iteratively adding the species with lowest max-similarity to
/// the selected hinge set.
///
/// # Errors
/// Returns an error if `num_hinges == 0` or `num_hinges > species_sketches.len()`.
pub fn select_auto_hinges_eisenstein(
    species_sketches: &[(String, EisensteinSamplerSketch)],
    num_hinges: usize,
) -> anyhow::Result<AutoHingePlan> {
    let n = species_sketches.len();
    anyhow::ensure!(num_hinges >= 1, "num_hinges must be at least 1");
    anyhow::ensure!(
        num_hinges <= n,
        "num_hinges ({num_hinges}) cannot exceed species count ({n})"
    );

    let mut chosen: Vec<usize> = Vec::with_capacity(num_hinges);
    let mut is_chosen = vec![false; n];
    // Deterministic seed hinge: index 0 (first assembly in input list).
    chosen.push(0);
    is_chosen[0] = true;

    // Track each species' maximum similarity to any chosen hinge (updated in O(N) per hinge!).
    let mut max_sim_to_hinges: Vec<f64> = (0..n)
        .map(|i| {
            if i == 0 {
                1.0
            } else {
                species_sketches[i].1.jaccard(&species_sketches[0].1)
            }
        })
        .collect();
    let mut best_hinge_idx: Vec<usize> = vec![0; n];

    while chosen.len() < num_hinges {
        let mut next_h = None;
        let mut lowest_sim = f64::INFINITY;
        for (i, &sim) in max_sim_to_hinges.iter().enumerate() {
            if !is_chosen[i] && sim < lowest_sim {
                lowest_sim = sim;
                next_h = Some(i);
            }
        }
        let Some(h_idx) = next_h else {
            break;
        };
        chosen.push(h_idx);
        is_chosen[h_idx] = true;
        max_sim_to_hinges[h_idx] = 1.0;
        best_hinge_idx[h_idx] = h_idx;

        for i in 0..n {
            if is_chosen[i] {
                continue;
            }
            let sim = species_sketches[i].1.jaccard(&species_sketches[h_idx].1);
            if sim > max_sim_to_hinges[i] {
                max_sim_to_hinges[i] = sim;
                best_hinge_idx[i] = h_idx;
            }
        }
    }

    let hinges: Vec<String> = chosen
        .iter()
        .map(|&idx| species_sketches[idx].0.clone())
        .collect();
    let mut tip_assignments = Vec::with_capacity(n.saturating_sub(hinges.len()));
    for i in 0..n {
        if !is_chosen[i] {
            tip_assignments.push((
                species_sketches[i].0.clone(),
                species_sketches[best_hinge_idx[i]].0.clone(),
            ));
        }
    }
    Ok((hinges, tip_assignments))
}
