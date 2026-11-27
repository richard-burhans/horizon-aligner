//! Lossless conversions between [`crate::formats::HnexArchive`] and [`crate::tba::TbaBlock`].

use crate::formats::{HnexArchive, HnexBlock, HnexBubble, HnexRow, HnexSeqMeta};
use crate::tba::{TbaBlock, TbaRow};
use std::collections::BTreeMap;

impl HnexArchive {
    /// Construct an `.hnex` archive from reference-free [`TbaBlock`]s, optional guide tree,
    /// and explicit Horizon Bubbles.
    #[must_use]
    pub fn from_tba_blocks(
        tba_blocks: &[TbaBlock],
        guide_tree_newick: Option<String>,
        bubbles: Vec<HnexBubble>,
    ) -> Self {
        let mut seq_map: BTreeMap<(String, String), usize> = BTreeMap::new();
        let mut blocks = Vec::with_capacity(tba_blocks.len());
        for (block_id, blk) in tba_blocks.iter().enumerate() {
            let mut rows = Vec::with_capacity(blk.rows.len());
            for r in &blk.rows {
                let entry = seq_map
                    .entry((r.species.clone(), r.contig.clone()))
                    .or_insert(0);
                *entry = (*entry).max(r.src_size);
                rows.push(HnexRow {
                    species: r.species.clone(),
                    contig: r.contig.clone(),
                    start: r.start,
                    size: r.size,
                    strand: r.strand,
                    src_size: r.src_size,
                    seq: String::from_utf8_lossy(&r.seq).into_owned(),
                });
            }
            blocks.push(HnexBlock {
                block_id,
                score_bits: blk.score,
                rows,
                col_posteriors: Vec::new(),
            });
        }
        for bub in &bubbles {
            for sp in &bub.spans {
                let entry = seq_map
                    .entry((sp.species.clone(), sp.contig.clone()))
                    .or_insert(0);
                *entry = (*entry).max(sp.end);
            }
        }
        let sequences = seq_map
            .into_iter()
            .map(|((species, contig), length)| HnexSeqMeta {
                species,
                contig,
                length,
            })
            .collect();
        Self {
            version: "1.0".to_owned(),
            guide_tree_newick,
            sequences,
            blocks,
            bubbles,
        }
    }

    /// Convert aligned blocks in this archive losslessly back to [`TbaBlock`]s (for MAF/GFA export).
    #[must_use]
    pub fn to_tba_blocks(&self) -> Vec<TbaBlock> {
        self.blocks
            .iter()
            .map(|blk| TbaBlock {
                score: blk.score_bits,
                rows: blk
                    .rows
                    .iter()
                    .map(|r| TbaRow {
                        species: r.species.clone(),
                        contig: r.contig.clone(),
                        start: r.start,
                        size: r.size,
                        strand: r.strand,
                        src_size: r.src_size,
                        seq: r.seq.as_bytes().to_vec(),
                    })
                    .collect(),
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::formats::{HnexArchive, HnexBubble, HnexBubbleSpan};
    use crate::tba::{TbaBlock, TbaRow};

    #[test]
    fn hnex_archive_roundtrips_blocks_posteriors_bubbles_and_converts_losslessly_to_tba() {
        let mk_row = |sp: &str, start: usize, seq: &str| TbaRow {
            species: sp.to_owned(),
            contig: "chr1".to_owned(),
            start,
            size: seq.bytes().filter(|&b| b != b'-').count(),
            strand: '+',
            src_size: 200,
            seq: seq.as_bytes().to_vec(),
        };
        let tba_blocks = vec![TbaBlock {
            score: 42.5,
            rows: vec![
                mk_row("spA", 0, "ACGTACGT--ACGT"),
                mk_row("spB", 0, "ACGTACGTTTACGT"),
            ],
        }];
        let bubbles = vec![HnexBubble {
            bubble_id: 1,
            left_block: Some(0),
            right_block: None,
            delta_l_bits: -14.25,
            reason: "mdl_abstention_hypervariable_loop".to_owned(),
            spans: vec![
                HnexBubbleSpan {
                    species: "spA".to_owned(),
                    contig: "chr1".to_owned(),
                    start: 12,
                    end: 20,
                    seq: "ATATATAT".to_owned(),
                },
                HnexBubbleSpan {
                    species: "spB".to_owned(),
                    contig: "chr1".to_owned(),
                    start: 14,
                    end: 22,
                    seq: "GCGCGCGC".to_owned(),
                },
            ],
        }];
        let mut archive =
            HnexArchive::from_tba_blocks(&tba_blocks, Some("(spA,spB);".to_owned()), bubbles);
        archive.blocks[0].col_posteriors = vec![0.99; 14];

        let mut buf = Vec::new();
        archive.write_to(&mut buf).expect("write hnex");
        let parsed = HnexArchive::read_from(buf.as_slice()).expect("read hnex");

        assert_eq!(parsed.version, "1.0");
        assert_eq!(parsed.guide_tree_newick.as_deref(), Some("(spA,spB);"));
        assert_eq!(parsed.blocks.len(), 1);
        assert_eq!(parsed.blocks[0].col_posteriors.len(), 14);
        assert_eq!(parsed.bubbles.len(), 1);
        assert_eq!(parsed.bubbles[0].spans.len(), 2);
        assert!((parsed.bubbles[0].delta_l_bits - (-14.25)).abs() < 1e-4);

        let rt_tba = parsed.to_tba_blocks();
        assert_eq!(rt_tba.len(), 1);
        assert_eq!(rt_tba[0].rows[0].seq, b"ACGTACGT--ACGT");
        assert_eq!(rt_tba[0].rows[1].seq, b"ACGTACGTTTACGT");

        // Bubble-only contigs must also be tracked in archive.sequences.
        let bub_only = HnexArchive::from_tba_blocks(
            &[],
            None,
            vec![HnexBubble {
                bubble_id: 0,
                left_block: None,
                right_block: None,
                delta_l_bits: -1.0,
                reason: "ok".to_owned(),
                spans: vec![HnexBubbleSpan {
                    species: "spC".to_owned(),
                    contig: "chr9".to_owned(),
                    start: 0,
                    end: 4,
                    seq: "ACGT".to_owned(),
                }],
            }],
        );
        assert!(
            bub_only
                .sequences
                .iter()
                .any(|s| s.species == "spC" && s.contig == "chr9"),
            "bubble-only contigs must appear in archive.sequences"
        );
    }
}
