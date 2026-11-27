//! Newick parser and UPGMA guide-tree construction for reference-free TBA alignment.

use anyhow::{Result, bail};

/// Binary species guide tree used to schedule bottom-up TBA blockset merges.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GuideTree {
    /// Leaf species name.
    Leaf(String),
    /// Internal merge node joining two clades.
    Internal {
        /// Left child clade.
        left: Box<Self>,
        /// Right child clade.
        right: Box<Self>,
    },
}

impl GuideTree {
    /// Collect all leaf species names in left-to-right order iteratively (stack-safe).
    #[must_use]
    pub fn leaves(&self) -> Vec<&str> {
        let mut out = Vec::new();
        let mut stack = vec![self];
        while let Some(node) = stack.pop() {
            match node {
                Self::Leaf(name) => out.push(name.as_str()),
                Self::Internal { left, right } => {
                    stack.push(right);
                    stack.push(left);
                }
            }
        }
        out
    }

    /// Format this guide tree as a Newick string terminated by `;` iteratively (stack-safe).
    #[must_use]
    pub fn to_newick(&self) -> String {
        enum Frame<'a> {
            Node(&'a GuideTree),
            Lit(char),
        }
        let mut out = String::new();
        let mut stack = vec![Frame::Node(self)];
        while let Some(frame) = stack.pop() {
            match frame {
                Frame::Lit(ch) => out.push(ch),
                Frame::Node(Self::Leaf(name)) => out.push_str(name),
                Frame::Node(Self::Internal { left, right }) => {
                    stack.push(Frame::Lit(')'));
                    stack.push(Frame::Node(right));
                    stack.push(Frame::Lit(','));
                    stack.push(Frame::Node(left));
                    stack.push(Frame::Lit('('));
                }
            }
        }
        out.push(';');
        out
    }
}

const MAX_NEWICK_DEPTH: usize = 512;

/// Parse a Newick guide tree string (ignoring optional branch lengths and internal node labels,
/// and left-binarizing any multifurcations).
///
/// # Errors
/// Returns an error if the Newick syntax is malformed, exceeds maximum nesting depth, or contains
/// empty/duplicate leaves.
pub fn parse_newick(input: &str) -> Result<GuideTree> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        bail!("empty Newick tree");
    }
    let bytes = trimmed.as_bytes();
    let mut pos = 0usize;
    let tree = parse_subtree(bytes, &mut pos, 0)?;
    skip_ws(bytes, &mut pos);
    if pos < bytes.len() && bytes[pos] == b';' {
        pos += 1;
    }
    skip_ws(bytes, &mut pos);
    if pos != bytes.len() {
        bail!(
            "unexpected trailing characters in Newick tree at byte {pos}: `{}`",
            &trimmed[pos..]
        );
    }
    let leaves = tree.leaves();
    let mut seen = std::collections::HashSet::with_capacity(leaves.len());
    for leaf in leaves {
        if !seen.insert(leaf) {
            bail!("duplicate species `{leaf}` in Newick tree");
        }
    }
    Ok(tree)
}

fn skip_ws(bytes: &[u8], pos: &mut usize) {
    while *pos < bytes.len() && bytes[*pos].is_ascii_whitespace() {
        *pos += 1;
    }
}

fn skip_branch_len(bytes: &[u8], pos: &mut usize) {
    skip_ws(bytes, pos);
    if *pos < bytes.len() && bytes[*pos] == b':' {
        *pos += 1;
        while *pos < bytes.len() && !matches!(bytes[*pos], b',' | b')' | b';' | b'(') {
            *pos += 1;
        }
    }
}

fn parse_label(bytes: &[u8], pos: &mut usize) -> Result<String> {
    skip_ws(bytes, pos);
    if *pos < bytes.len() && (bytes[*pos] == b'\'' || bytes[*pos] == b'"') {
        let quote = bytes[*pos];
        *pos += 1;
        let start = *pos;
        while *pos < bytes.len() && bytes[*pos] != quote {
            *pos += 1;
        }
        if *pos >= bytes.len() {
            bail!("unterminated quoted label in Newick tree");
        }
        let label = std::str::from_utf8(&bytes[start..*pos])
            .map_err(|e| anyhow::anyhow!("invalid UTF-8 in Newick label: {e}"))?
            .to_owned();
        *pos += 1;
        return Ok(label);
    }
    let start = *pos;
    while *pos < bytes.len()
        && !bytes[*pos].is_ascii_whitespace()
        && !matches!(bytes[*pos], b'(' | b')' | b',' | b':' | b';')
    {
        *pos += 1;
    }
    let label = std::str::from_utf8(&bytes[start..*pos])
        .map_err(|e| anyhow::anyhow!("invalid UTF-8 in Newick label: {e}"))?
        .to_owned();
    Ok(label)
}

fn parse_subtree(bytes: &[u8], pos: &mut usize, depth: usize) -> Result<GuideTree> {
    if depth > MAX_NEWICK_DEPTH {
        bail!("Newick tree exceeds maximum nesting depth ({MAX_NEWICK_DEPTH})");
    }
    skip_ws(bytes, pos);
    if *pos >= bytes.len() {
        bail!("unexpected end of Newick tree");
    }
    if bytes[*pos] == b'(' {
        *pos += 1;
        let mut children = Vec::new();
        loop {
            children.push(parse_subtree(bytes, pos, depth + 1)?);
            skip_ws(bytes, pos);
            if *pos >= bytes.len() {
                bail!("unclosed `(` in Newick tree");
            }
            if bytes[*pos] == b',' {
                *pos += 1;
            } else if bytes[*pos] == b')' {
                *pos += 1;
                break;
            } else {
                bail!(
                    "expected `,` or `)` in Newick tree at byte {}, got `{}`",
                    *pos,
                    bytes[*pos] as char
                );
            }
        }
        // Optional internal node label and branch length.
        let _ = parse_label(bytes, pos)?;
        skip_branch_len(bytes, pos);
        if depth.saturating_add(children.len().saturating_sub(1)) > MAX_NEWICK_DEPTH {
            bail!(
                "Newick multifurcation binarized depth exceeds maximum nesting depth ({MAX_NEWICK_DEPTH})"
            );
        }
        let mut iter = children.into_iter();
        let Some(mut acc) = iter.next() else {
            bail!("empty internal node in Newick tree");
        };
        for child in iter {
            acc = GuideTree::Internal {
                left: Box::new(acc),
                right: Box::new(child),
            };
        }
        Ok(acc)
    } else {
        let name = parse_label(bytes, pos)?;
        if name.is_empty() {
            bail!("empty leaf name in Newick tree at byte {}", *pos);
        }
        skip_branch_len(bytes, pos);
        Ok(GuideTree::Leaf(name))
    }
}

/// Build a deterministic UPGMA binary guide tree from species names and a symmetric similarity
/// matrix `sim[i][j]` (higher values indicate closer species).
///
/// # Errors
/// Returns an error if `names` is empty or `sim` dimensions do not match `names.len()`.
pub fn build_upgma_tree(names: &[String], sim: &[Vec<f64>]) -> Result<GuideTree> {
    let n = names.len();
    if n == 0 {
        bail!("cannot build UPGMA tree from 0 species");
    }
    if sim.len() != n || sim.iter().any(|row| row.len() != n) {
        bail!("UPGMA similarity matrix dimensions do not match species count {n}");
    }
    let mut seen_names = std::collections::HashSet::with_capacity(n);
    for name in names {
        if name.trim().is_empty() {
            bail!("UPGMA species name must not be empty");
        }
        if !seen_names.insert(name.as_str()) {
            bail!("duplicate species `{name}` in UPGMA input");
        }
    }
    for (i, row) in sim.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            if !v.is_finite() {
                bail!("non-finite similarity value `{v}` at ({i}, {j}) in UPGMA matrix");
            }
        }
    }
    let mut clades: Vec<(GuideTree, Vec<usize>)> = names
        .iter()
        .enumerate()
        .map(|(i, name)| (GuideTree::Leaf(name.clone()), vec![i]))
        .collect();

    while clades.len() >= 2 {
        let mut best: Option<(f64, usize, usize)> = None;
        for a in 0..clades.len() {
            for b in (a + 1)..clades.len() {
                let mut sum = 0.0_f64;
                let mut cnt = 0usize;
                for &i in &clades[a].1 {
                    for &j in &clades[b].1 {
                        sum += sim[i][j];
                        cnt += 1;
                    }
                }
                #[allow(clippy::cast_precision_loss)]
                let avg = if cnt == 0 { 0.0 } else { sum / (cnt as f64) };
                if best.is_none_or(|(best_sim, _, _)| {
                    avg.total_cmp(&best_sim) == std::cmp::Ordering::Greater
                }) {
                    best = Some((avg, a, b));
                }
            }
        }
        let Some((_, a, b)) = best else {
            break;
        };
        let (right_tree, right_members) = clades.remove(b);
        let (left_tree, mut left_members) = clades.remove(a);
        left_members.extend(right_members);
        let merged = GuideTree::Internal {
            left: Box::new(left_tree),
            right: Box::new(right_tree),
        };
        clades.push((merged, left_members));
    }

    clades
        .into_iter()
        .next()
        .map(|(tree, _)| tree)
        .ok_or_else(|| anyhow::anyhow!("UPGMA produced empty tree"))
}

/// Build an `N x N` symmetric similarity matrix in `[0.0, 1.0]` from pairwise TBA blocks by
/// normalizing total aligned bases against the shorter genome's observed contig length sum.
#[must_use]
pub fn aligned_bp_similarity_matrix<S: std::hash::BuildHasher>(
    species: &[String],
    pair_blocks: &std::collections::HashMap<(String, String), Vec<super::PairwiseTbaBlock>, S>,
) -> Vec<Vec<f64>> {
    let n = species.len();
    let sp_idx: std::collections::HashMap<&str, usize> = species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut contig_lens: Vec<std::collections::HashMap<&str, usize>> =
        vec![std::collections::HashMap::new(); n];
    for blocks in pair_blocks.values() {
        for blk in blocks {
            for row in [&blk.row_a, &blk.row_b] {
                if let Some(&si) = sp_idx.get(row.species.as_str()) {
                    let entry = contig_lens[si].entry(row.contig.as_str()).or_insert(0);
                    *entry = (*entry).max(row.src_size);
                }
            }
        }
    }
    let sp_total_len: Vec<usize> = contig_lens
        .iter()
        .map(|m| m.values().copied().sum())
        .collect();
    let mut sim = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        sim[i][i] = 1.0;
        for j in (i + 1)..n {
            let key = super::canonical_pair_key(&species[i], &species[j]);
            let bp: usize = pair_blocks.get(&key).map_or(0, |v| {
                v.iter().map(|b| b.row_a.size.min(b.row_b.size)).sum()
            });
            let denom = sp_total_len[i].min(sp_total_len[j]).max(1);
            #[allow(clippy::cast_precision_loss)]
            let val = ((bp as f64) / (denom as f64)).clamp(0.0, 1.0);
            sim[i][j] = val;
            sim[j][i] = val;
        }
    }
    sim
}

/// Extract the `M x M` similarity submatrix induced by `subset` (using `sp_to_idx` lookups into
/// `full_sim`).
#[must_use]
pub fn similarity_submatrix<S: std::hash::BuildHasher>(
    full_sim: &[Vec<f64>],
    subset: &[String],
    sp_to_idx: &std::collections::HashMap<&str, usize, S>,
) -> Vec<Vec<f64>> {
    let m = subset.len();
    let mut sub = vec![vec![0.0_f64; m]; m];
    for a in 0..m {
        let ia = sp_to_idx[subset[a].as_str()];
        sub[a][a] = 1.0;
        for b in (a + 1)..m {
            let ib = sp_to_idx[subset[b].as_str()];
            sub[a][b] = full_sim[ia][ib];
            sub[b][a] = full_sim[ib][ia];
        }
    }
    sub
}

/// Build a `BioNJ` (Gascuel 1997) variance-weighted neighbor-joining guide tree from a symmetric
/// pairwise distance matrix `dist[i][j]` (`>= 0.0`, finite, zero diagonal).
///
/// # Errors
/// Returns an error if `species` is empty, contains duplicates, or `dist` is malformed/non-finite.
#[allow(clippy::too_many_lines)]
pub fn build_bionj_tree(species: &[String], dist: &[Vec<f64>]) -> Result<GuideTree> {
    let n_sp = species.len();
    if n_sp == 0 {
        bail!("cannot build BioNJ guide tree from zero species");
    }
    if dist.len() != n_sp || dist.iter().any(|row| row.len() != n_sp) {
        bail!("BioNJ distance matrix must be {n_sp}x{n_sp}");
    }
    let mut seen = std::collections::HashSet::with_capacity(n_sp);
    for sp in species {
        if sp.trim().is_empty() {
            bail!("empty species name in BioNJ input");
        }
        if !seen.insert(sp.as_str()) {
            bail!("duplicate species `{sp}` in BioNJ input");
        }
    }
    for (i, row) in dist.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            if !val.is_finite() || val < 0.0 {
                bail!("invalid non-finite or negative BioNJ distance at ({i}, {j}): {val}");
            }
        }
    }
    if n_sp == 1 {
        return Ok(GuideTree::Leaf(species[0].clone()));
    }

    let mut nodes: Vec<GuideTree> = species.iter().cloned().map(GuideTree::Leaf).collect();
    let mut dmat: Vec<Vec<f64>> = (0..n_sp)
        .map(|i| {
            (0..n_sp)
                .map(|j| {
                    if i == j {
                        0.0
                    } else {
                        f64::midpoint(dist[i][j], dist[j][i])
                    }
                })
                .collect()
        })
        .collect();
    let mut vmat = dmat.clone();

    while nodes.len() > 2 {
        let active = nodes.len();
        #[allow(clippy::cast_precision_loss)]
        let r_m2 = (active - 2) as f64;
        let row_sums: Vec<f64> = dmat.iter().map(|row| row.iter().sum()).collect();

        let mut best_pair = (0usize, 1usize);
        let mut best_q = f64::INFINITY;
        for i in 0..active {
            for j in (i + 1)..active {
                let q_ij = r_m2 * dmat[i][j] - row_sums[i] - row_sums[j];
                if q_ij < best_q {
                    best_q = q_ij;
                    best_pair = (i, j);
                }
            }
        }
        let (idx_a, idx_b) = best_pair;
        let d_ab = dmat[idx_a][idx_b];
        let v_ab = vmat[idx_a][idx_b];
        // Use unclamped branch estimates in distance reduction (Gascuel 1997 eq. 4) to preserve
        // exact additivity and avoid asymmetric bias on noisy matrices.
        let branch_a = 0.5 * d_ab + (row_sums[idx_a] - row_sums[idx_b]) / (2.0 * r_m2);
        let branch_b = d_ab - branch_a;

        // Gascuel (1997) variance-minimizing reduction weight λ ∈ [0, 1].
        let lambda = if v_ab <= 1e-12 {
            0.5
        } else {
            let mut diff_sum = 0.0_f64;
            for k in 0..active {
                if k != idx_a && k != idx_b {
                    diff_sum += vmat[idx_b][k] - vmat[idx_a][k];
                }
            }
            (0.5 + diff_sum / (2.0 * r_m2 * v_ab)).clamp(0.0, 1.0)
        };

        let mut next_dists = Vec::with_capacity(active - 1);
        let mut next_vars = Vec::with_capacity(active - 1);
        for k in 0..active {
            if k == idx_a || k == idx_b {
                continue;
            }
            let d_uk = (lambda * (dmat[idx_a][k] - branch_a)
                + (1.0 - lambda) * (dmat[idx_b][k] - branch_b))
                .max(0.0);
            let v_uk = (lambda * vmat[idx_a][k] + (1.0 - lambda) * vmat[idx_b][k]
                - lambda * (1.0 - lambda) * v_ab)
                .max(0.0);
            next_dists.push(d_uk);
            next_vars.push(v_uk);
        }

        // Remove idx_b then idx_a (idx_b > idx_a).
        let right = nodes.remove(idx_b);
        let left = nodes.remove(idx_a);
        nodes.push(GuideTree::Internal {
            left: Box::new(left),
            right: Box::new(right),
        });

        dmat.remove(idx_b);
        dmat.remove(idx_a);
        vmat.remove(idx_b);
        vmat.remove(idx_a);
        for (idx, (drow, vrow)) in dmat.iter_mut().zip(vmat.iter_mut()).enumerate() {
            drow.remove(idx_b);
            drow.remove(idx_a);
            vrow.remove(idx_b);
            vrow.remove(idx_a);
            drow.push(next_dists[idx]);
            vrow.push(next_vars[idx]);
        }
        next_dists.push(0.0);
        next_vars.push(0.0);
        dmat.push(next_dists);
        vmat.push(next_vars);
    }

    let (Some(right), Some(left)) = (nodes.pop(), nodes.pop()) else {
        bail!("internal BioNJ error: expected two final clusters");
    };
    Ok(GuideTree::Internal {
        left: Box::new(left),
        right: Box::new(right),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligned_bp_similarity_matrix_normalizes_by_shorter_species_length() {
        let species = vec!["sp1".to_owned(), "sp2".to_owned()];
        let maf = "a score=10\ns chr1 0 50 + 100 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\ns chr1 0 50 + 200 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\n";
        let blks = super::super::parse_pairwise_tba_reader(
            "sp1",
            "sp2",
            std::io::Cursor::new(maf.as_bytes()),
        )
        .expect("parse");
        let mut map = std::collections::HashMap::new();
        map.insert(super::super::canonical_pair_key("sp1", "sp2"), blks);
        let sim = aligned_bp_similarity_matrix(&species, &map);
        assert!((sim[0][0] - 1.0).abs() < 1e-12);
        assert!((sim[0][1] - 0.5).abs() < 1e-12);
        assert!((sim[1][0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn parses_newick_with_branch_lengths_labels_and_multifurcations() {
        let tree =
            parse_newick("((sp1:0.1,sp2:0.2)Anc0:0.05,sp3:0.3,sp4:0.4)Root;").expect("parse");
        assert_eq!(tree.leaves(), vec!["sp1", "sp2", "sp3", "sp4"]);
        assert_eq!(tree.to_newick(), "(((sp1,sp2),sp3),sp4);");
    }

    #[test]
    fn upgma_clusters_closest_clades_first() {
        let names = vec![
            "lav1".to_owned(),
            "lav2".to_owned(),
            "vin1".to_owned(),
            "vin2".to_owned(),
        ];
        let sim = vec![
            vec![1.0, 0.9, 0.1, 0.1],
            vec![0.9, 1.0, 0.1, 0.1],
            vec![0.1, 0.1, 1.0, 0.85],
            vec![0.1, 0.1, 0.85, 1.0],
        ];
        let tree = build_upgma_tree(&names, &sim).expect("upgma");
        assert_eq!(tree.to_newick(), "((lav1,lav2),(vin1,vin2));");
    }

    #[test]
    fn build_upgma_tree_rejects_nan_and_duplicate_or_empty_names() {
        let names = vec!["a".to_owned(), "b".to_owned(), "c".to_owned()];
        let nan_sim = vec![
            vec![1.0, f64::NAN, 0.2],
            vec![f64::NAN, 1.0, 0.3],
            vec![0.2, 0.3, 1.0],
        ];
        assert!(build_upgma_tree(&names, &nan_sim).is_err());

        let dup_names = vec!["a".to_owned(), "a".to_owned()];
        let ok_sim = vec![vec![1.0, 0.5], vec![0.5, 1.0]];
        assert!(build_upgma_tree(&dup_names, &ok_sim).is_err());

        let empty_names = vec![String::new(), "b".to_owned()];
        assert!(build_upgma_tree(&empty_names, &ok_sim).is_err());
    }

    #[test]
    fn parse_newick_rejects_open_paren_inside_branch_length_and_deep_nesting() {
        assert!(parse_newick("(a:0.1(b);").is_err());
        let mut deep = String::new();
        for _ in 0..600 {
            deep.push('(');
        }
        deep.push('a');
        for _ in 0..600 {
            deep.push(')');
        }
        deep.push(';');
        assert!(parse_newick(&deep).is_err());
    }

    #[test]
    fn bionj_recovers_non_ultrametric_additive_tree_topology_where_upgma_fails() {
        // Classic non-ultrametric 4-taxon additive tree: ((A:0.05, B:0.60):0.25, (C:0.05, D:0.05):0.25)
        // True split is {A, B} | {C, D}, even though d(A, C) = 0.60 < d(A, B) = 0.65 (so UPGMA wrongly joins A with C/D!).
        // BioNJ (variance-weighted NJ) MUST recover ((A,B),(C,D))!
        let names = vec![
            "A".to_owned(),
            "B".to_owned(),
            "C".to_owned(),
            "D".to_owned(),
        ];
        let dist = vec![
            vec![0.00, 0.65, 0.60, 0.60],
            vec![0.65, 0.00, 1.15, 1.15],
            vec![0.60, 1.15, 0.00, 0.10],
            vec![0.60, 1.15, 0.10, 0.00],
        ];
        let tree = build_bionj_tree(&names, &dist).expect("bionj");
        let nwk = tree.to_newick();
        // Either ((A,B),(C,D)) or ((C,D),(A,B)) has the exact split {A,B} | {C,D}.
        let has_ab = nwk.contains("(A,B)") || nwk.contains("(B,A)");
        let has_cd = nwk.contains("(C,D)") || nwk.contains("(D,C)");
        assert!(
            has_ab && has_cd,
            "BioNJ must recover non-ultrametric split {{A,B}}|{{C,D}}, got {nwk}"
        );
    }

    #[test]
    fn parse_newick_rejects_flat_multifurcation_exceeding_max_depth() {
        let leaves: Vec<String> = (0..600).map(|i| format!("s{i}")).collect();
        let flat = format!("({});", leaves.join(","));
        assert!(
            parse_newick(&flat).is_err(),
            "flat multifurcation exceeding MAX_NEWICK_DEPTH must be rejected"
        );
    }

    #[test]
    fn bionj_handles_single_taxon_nonzero_diag_and_rejects_whitespace_names() {
        let one = build_bionj_tree(&["A".to_owned()], &[vec![0.0]]).expect("single taxon");
        assert_eq!(one.to_newick(), "A;");

        // Non-zero diagonal is zeroed defensively; whitespace names rejected.
        let names = vec!["A".to_owned(), "B".to_owned(), "C".to_owned()];
        let dist = vec![
            vec![0.5, 0.2, 0.8],
            vec![0.2, 0.5, 0.8],
            vec![0.8, 0.8, 0.5],
        ];
        let t = build_bionj_tree(&names, &dist).expect("bionj 3-taxon");
        assert_eq!(t.leaves().len(), 3);
        assert!(build_bionj_tree(&["  ".to_owned()], &[vec![0.0]]).is_err());
    }
}
