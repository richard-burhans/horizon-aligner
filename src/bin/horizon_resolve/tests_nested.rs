use super::*;

#[test]
#[allow(clippy::too_many_lines)]
fn nested_snarls_handles_inverted_interior_rev_ref_duplicated_ref_and_phased_diploid() -> Result<()>
{
    let dir = std::env::temp_dir().join(format!(
        "horizon_resolve_nested_edges_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir)?;

    // Case 1: Inverted interior segment (`>m1<inv>m2` where REF has `>m1>inv>m2`), phased diploid
    // sample `DIP` carrying alleles `1|2` (`inv_alt` on hap1, `<inv` on hap2) and `MISS` carrying
    // `.|.` (spanning deletion `>a>z` on both haps), plus a duplicated reference segment `dup_s`
    // that appears twice in REF (`REF_DUPLICATED`) and must NOT spawn bogus nested bubbles.
    let gfa_in = dir.join("nested_inv.gfa");
    // `inv` = "AACC" -> revcomp(`<inv`) = "GGTT" (distinct from "AACC" and "TGGT"!).
    std::fs::write(
        &gfa_in,
        "H\tVN:Z:1.1\n\
         S\ta\tACGTACGTAC\n\
         S\tdup_s\tGAGA\n\
         S\tm1\tTTTTT\n\
         S\tinv\tAACC\n\
         S\tinv_alt\tTGGT\n\
         S\tm2\tCCCCC\n\
         S\tz\tGTGTGTGTGT\n\
         L\ta\t+\tdup_s\t+\t0M\n\
         L\tdup_s\t+\tm1\t+\t0M\n\
         L\ta\t+\tz\t+\t0M\n\
         L\tm1\t+\tinv\t+\t0M\n\
         L\tm1\t+\tinv\t-\t0M\n\
         L\tm1\t+\tinv_alt\t+\t0M\n\
         L\tinv\t+\tm2\t+\t0M\n\
         L\tinv\t-\tm2\t+\t0M\n\
         L\tinv_alt\t+\tm2\t+\t0M\n\
         L\tm2\t+\tdup_s\t+\t0M\n\
         L\tdup_s\t+\tz\t+\t0M\n\
         W\tREF\t0\tchr1\t0\t42\t>a>dup_s>m1>inv>m2>dup_s>z\n\
         W\tDIP\t1\tchr1\t0\t42\t>a>dup_s>m1>inv_alt>m2>dup_s>z\n\
         W\tDIP\t2\tchr1\t0\t42\t>a>dup_s>m1<inv>m2>dup_s>z\n\
         W\tMISS\t1\tchr1\t0\t20\t>a>z\n\
         W\tMISS\t2\tchr1\t0\t20\t>a>z\n",
    )?;
    let cli = CliArgs {
        gfa_in: gfa_in.to_string_lossy().to_string(),
        out_path: None,
        vcf_path: Some("out.vcf".to_owned()),
        vcf_ref: Some("REF".to_owned()),
        nested_snarls: true,
        collapse_dups: false,
        tag_walks: false,
        max_bubble_steps: 16,
        min_anchor_bp: 100,
        divergence: None,
        threads: None,
    };
    let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
    let mut gfa_buf = Vec::new();
    let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
    let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;

    let mut vcf_buf = Vec::new();
    let vstats = export_snarl_vcf(
        &graph,
        &cap.seg_seqs,
        &cap.walks,
        Some("REF"),
        true,
        64,
        &mut vcf_buf,
    )?;
    let text = String::from_utf8(vcf_buf)?;
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    // Expect 1 LV=0 spanning deletion (a -> z) + 1 LV=1 nested bubble (m1 -> m2) with 3 alleles
    // (REF=AACC, ALT1=TGGT, ALT2=GGTT from `<inv`). `dup_s` is REF_DUPLICATED so it never anchors!
    assert_eq!(vstats.nested_sites, 1, "expected 1 LV=1 site:\n{text}");
    let lv1 = rows
        .iter()
        .find(|r| r[7].ends_with(";LV=1;PS=snarl_1"))
        .expect("LV=1 row with PS=snarl_1");
    // Samples sorted alphabetically: DIP(9), MISS(10), REF(11)
    // POS = len(a=10) + len(dup_s=4) + len(m1=5) + 1 = 20
    assert_eq!(lv1[1], "20", "lv1: {lv1:?}");
    assert_eq!(lv1[3], "AACC", "lv1: {lv1:?}");
    assert_eq!(lv1[4], "TGGT,GGTT", "lv1: {lv1:?}");
    assert_eq!(lv1[9], "1|2", "phased diploid DIP must be 1|2:\n{text}");
    assert_eq!(
        lv1[10], ".|.",
        "spanning-deletion MISS must be .|.:\n{text}"
    );
    assert_eq!(lv1[11], "0", "REF must be 0:\n{text}");

    // Case 2: Reference walk traverses reverse-oriented segments `<m1<r_snp<m2` inside a spanning
    // deletion. `m1`="AAAAA" -> `<m1`="TTTTT" (last base 'T'); `r_snp`="C" -> `<r_snp`="G";
    // `a_snp`="T" -> `<a_snp`="A". Verify exact revcomp REF/ALT and left-anchor base on deletion!
    let gfa_rev = dir.join("nested_rev_ref.gfa");
    std::fs::write(
        &gfa_rev,
        "H\tVN:Z:1.1\n\
         S\ta\tACGTACGTAC\n\
         S\tm1\tGAAAA\n\
         S\tr_snp\tC\n\
         S\ta_snp\tT\n\
         S\tm2\tCCCCC\n\
         S\tz\tGTGTGTGTGT\n\
         L\ta\t+\tm1\t-\t0M\n\
         L\ta\t+\tz\t+\t0M\n\
         L\tm1\t-\tr_snp\t-\t0M\n\
         L\tm1\t-\ta_snp\t-\t0M\n\
         L\tm1\t-\tm2\t-\t0M\n\
         L\tr_snp\t-\tm2\t-\t0M\n\
         L\ta_snp\t-\tm2\t-\t0M\n\
         L\tm2\t-\tz\t+\t0M\n\
         W\tREF\t0\tchr1\t0\t31\t>a<m1<r_snp<m2>z\n\
         W\tALT\t0\tchr1\t0\t31\t>a<m1<a_snp<m2>z\n\
         W\tIDEL\t0\tchr1\t0\t30\t>a<m1<m2>z\n\
         W\tSPAN\t0\tchr1\t0\t20\t>a>z\n",
    )?;
    let cli_rev = CliArgs {
        gfa_in: gfa_rev.to_string_lossy().to_string(),
        out_path: None,
        vcf_path: Some("out.vcf".to_owned()),
        vcf_ref: Some("REF".to_owned()),
        nested_snarls: true,
        collapse_dups: false,
        tag_walks: false,
        max_bubble_steps: 16,
        min_anchor_bp: 100,
        divergence: None,
        threads: None,
    };
    let graph_r = pass1_index(&cli_rev.gfa_in, cli_rev.min_anchor_bp)?;
    let mut gfa_buf_r = Vec::new();
    let mut cap_r = resolve::VcfCapture::with_segments(graph_r.num_segs());
    let _ = pass2_rewrite(
        &cli_rev.gfa_in,
        &mut gfa_buf_r,
        &graph_r,
        &cli_rev,
        Some(&mut cap_r),
    )?;
    let mut vcf_r = Vec::new();
    let vstats_r = export_snarl_vcf(
        &graph_r,
        &cap_r.seg_seqs,
        &cap_r.walks,
        Some("REF"),
        true,
        64,
        &mut vcf_r,
    )?;
    let _ = std::fs::remove_dir_all(&dir);
    let text_r = String::from_utf8(vcf_r)?;
    assert_eq!(
        vstats_r.nested_sites, 1,
        "expected 1 LV=1 on rev-ref:\n{text_r}"
    );
    let lv1_r: Vec<&str> = text_r
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .find(|r| r[7].ends_with(";LV=1;PS=snarl_1"))
        .expect("LV=1 rev-ref row with PS=snarl_1");
    // `<m1` = revcomp("GAAAA") = "TTTTC", so last oriented base of `<m1` is 'C'!
    // Since IDEL has empty interior, left-anchor 'C' is prepended at POS = 10 + 5 = 15:
    // REF = "C" + revcomp("C") = "CG", ALT = "CA,C" (ALT=<a_snp="A", IDEL="").
    assert_eq!(lv1_r[1], "15", "rev-ref anchored POS:\n{text_r}");
    assert_eq!(lv1_r[3], "CG", "rev-ref anchored REF:\n{text_r}");
    assert_eq!(lv1_r[4], "CA,C", "rev-ref anchored ALT:\n{text_r}");
    Ok(())
}

#[test]
fn per_contig_streaming_vcf_preserves_sorted_chrom_order_and_sequential_snarl_ids() -> Result<()> {
    let dir = std::env::temp_dir().join(format!(
        "horizon_resolve_multichrom_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir)?;
    let gfa_in = dir.join("mc.gfa");
    // Deliberately list chrB before chrA in the GFA; per-contig streaming must emit chrA first
    // (snarl_1) and chrB second (snarl_2), both with nested LV=1 SNPs inside spanning deletions.
    std::fs::write(
        &gfa_in,
        "H\tVN:Z:1.1\n\
         S\tb1\tACGTACGTAC\n\
         S\tbm1\tTTTTT\n\
         S\tbr\tA\n\
         S\tba\tG\n\
         S\tbm2\tCCCCC\n\
         S\tb2\tGTGTGTGTGT\n\
         S\ta1\tACGTACGTAC\n\
         S\tam1\tAAAAA\n\
         S\tar\tC\n\
         S\taa\tT\n\
         S\tam2\tGGGGG\n\
         S\ta2\tACACACACAC\n\
         W\tREF\t0\tchrB\t0\t31\t>b1>bm1>br>bm2>b2\n\
         W\tALT\t0\tchrB\t0\t31\t>b1>bm1>ba>bm2>b2\n\
         W\tDEL\t0\tchrB\t0\t20\t>b1>b2\n\
         W\tREF\t0\tchrA\t0\t31\t>a1>am1>ar>am2>a2\n\
         W\tALT\t0\tchrA\t0\t31\t>a1>am1>aa>am2>a2\n\
         W\tDEL\t0\tchrA\t0\t20\t>a1>a2\n",
    )?;
    let cli = CliArgs {
        gfa_in: gfa_in.to_string_lossy().to_string(),
        out_path: None,
        vcf_path: Some("out.vcf".to_owned()),
        vcf_ref: Some("REF".to_owned()),
        nested_snarls: true,
        collapse_dups: false,
        tag_walks: false,
        max_bubble_steps: 16,
        min_anchor_bp: 100,
        divergence: None,
        threads: Some(2),
    };
    let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
    let mut gfa_buf = Vec::new();
    let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
    let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;
    let mut vcf_buf = Vec::new();
    let vstats = export_snarl_vcf(
        &graph,
        &cap.seg_seqs,
        &cap.walks,
        Some("REF"),
        true,
        64,
        &mut vcf_buf,
    )?;
    let _ = std::fs::remove_dir_all(&dir);
    let text = String::from_utf8(vcf_buf)?;
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(
        vstats.sites_emitted, 4,
        "2 LV=0 + 2 LV=1 across chrA/chrB:\n{text}"
    );
    assert_eq!(vstats.nested_sites, 2);
    assert_eq!(rows[0][0], "chrA");
    assert_eq!(rows[1][0], "chrA");
    assert_eq!(rows[2][0], "chrB");
    assert_eq!(rows[3][0], "chrB");
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(row[2], format!("snarl_{}", idx + 1));
    }
    Ok(())
}

#[test]
#[allow(clippy::too_many_lines)]
fn vcf_export_sanitizes_iupac_bases_deduplicates_collapsed_alleles_and_suppresses_monomorphic()
-> Result<()> {
    let dir = std::env::temp_dir().join(format!(
        "horizon_resolve_iupac_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir)?;
    let gfa_in = dir.join("iupac.gfa");
    // Four snarls along chr1 exercising every IUPAC sanitization edge case:
    // 1. a -> b: REF=`CTTACTCTTCAT`, S1=`CTTACTYTTCAT`, S2=`CTTACTrTTCAT` (lowercase `r`),
    //    S3=`<rev_y` where `rev_y`=`ATGAAkAGTAAG` -> revcomp has `k` -> `N` = `CTTACTNTTCAT`!
    //    All three ALT haplotypes collapse to the single VCFv4.2 ALT `CTTACTNTTCAT` (AC=3, AN=4).
    // 2. b (ends in `Y`) -> c: deletion on S1/S2/S3 (empty interior) vs REF=`G` ->
    //    left-anchor base from `b` (`Y`) must sanitize to `N`, emitting REF=`NG`, ALT=`N`.
    // 3. c -> d: REF=`W` (-> `N`) vs S1=`Y` (-> `N`), S2=`R` (-> `N`), S3=`n` (-> `N`) ->
    //    all alleles collapse to `N`; site is monomorphic and MUST be suppressed!
    // 4. d -> e: REF=`M` (-> `N`) vs S1=`A`, S2=`C`, S3=`K` (-> `N` == REF) ->
    //    REF=`N`, ALT=`A,C`, S3 genotyped as `0` (REF), S1=`1`, S2=`2`.
    std::fs::write(
        &gfa_in,
        "H\tVN:Z:1.1\n\
         S\ta\tACGTACGTAC\n\
         S\tr1\tCTTACTCTTCAT\n\
         S\ty1\tCTTACTYTTCAT\n\
         S\tr_low\tCTTACTrTTCAT\n\
         S\trev_y\tATGAAkAGTAAG\n\
         S\tb\tACGTACGTAY\n\
         S\tg1\tG\n\
         S\tc\tTTTTTAAAAA\n\
         S\tw_ref\tW\n\
         S\ty_alt\tY\n\
         S\tr_alt\tR\n\
         S\tn_alt\tn\n\
         S\td\tCCCCCGGGGG\n\
         S\tm_ref\tM\n\
         S\ta_alt\tA\n\
         S\tc_alt\tC\n\
         S\tk_alt\tK\n\
         S\te\tACACACACAC\n\
         W\tREF\t0\tchr1\t0\t65\t>a>r1>b>g1>c>w_ref>d>m_ref>e\n\
         W\tS1\t0\tchr1\t0\t64\t>a>y1>b>c>y_alt>d>a_alt>e\n\
         W\tS2\t0\tchr1\t0\t64\t>a>r_low>b>c>r_alt>d>c_alt>e\n\
         W\tS3\t0\tchr1\t0\t64\t>a<rev_y>b>c>n_alt>d>k_alt>e\n",
    )?;
    let cli = CliArgs {
        gfa_in: gfa_in.to_string_lossy().to_string(),
        out_path: None,
        vcf_path: Some("out.vcf".to_owned()),
        vcf_ref: Some("REF".to_owned()),
        nested_snarls: false,
        collapse_dups: false,
        tag_walks: false,
        max_bubble_steps: 16,
        min_anchor_bp: 100,
        divergence: None,
        threads: None,
    };
    let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
    let mut gfa_buf = Vec::new();
    let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
    let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;

    let mut vcf_buf = Vec::new();
    let vstats = export_snarl_vcf(
        &graph,
        &cap.seg_seqs,
        &cap.walks,
        Some("REF"),
        false,
        64,
        &mut vcf_buf,
    )?;
    let _ = std::fs::remove_dir_all(&dir);
    let text = String::from_utf8(vcf_buf)?;
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();

    // Snarl 3 (W vs Y/R/n) collapsed to monomorphic N -> only 3 sites must be emitted.
    assert_eq!(
        vstats.sites_emitted, 3,
        "expected 3 polymorphic sites after IUPAC sanitization:\n{text}"
    );
    assert_eq!(rows.len(), 3, "vcf:\n{text}");

    // Every emitted REF and ALT character must belong strictly to {A, C, G, T, N, ','}.
    for row in &rows {
        assert!(
            row[3]
                .bytes()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')),
            "non-ACGTN REF allele in row: {row:?}"
        );
        assert!(
            row[4]
                .bytes()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b',')),
            "non-ACGTN ALT allele in row: {row:?}"
        );
    }

    // Samples sorted alphabetically: REF(9), S1(10), S2(11), S3(12).
    // Row 0 (a -> b): S1/S2/S3 all collapse to CTTACTNTTCAT (allele 1, AC=3, AN=4, NS=4).
    assert_eq!(rows[0][1], "11");
    assert_eq!(rows[0][3], "CTTACTCTTCAT");
    assert_eq!(rows[0][4], "CTTACTNTTCAT");
    assert!(
        rows[0][7].contains(";AC=3;AN=4;NS=4;LV=0"),
        "row0 info: {}",
        rows[0][7]
    );
    assert_eq!(rows[0][9], "0");
    assert_eq!(rows[0][10], "1");
    assert_eq!(rows[0][11], "1");
    assert_eq!(rows[0][12], "1");

    // Row 1 (b -> c): left-anchor `Y` sanitized to `N` -> REF=NG, ALT=N.
    assert_eq!(rows[1][1], "32");
    assert_eq!(rows[1][3], "NG");
    assert_eq!(rows[1][4], "N");
    assert!(
        rows[1][7].contains(";AC=3;AN=4;NS=4;LV=0"),
        "row1 info: {}",
        rows[1][7]
    );

    // Row 2 (d -> e): REF=`M`->`N`, S1=`A`(1), S2=`C`(2), S3=`K`->`N`(0, collapses to REF).
    assert_eq!(rows[2][1], "55");
    assert_eq!(rows[2][3], "N");
    assert_eq!(rows[2][4], "A,C");
    assert!(
        rows[2][7].contains(";AC=1,1;AN=4;NS=4;LV=0"),
        "row2 info: {}",
        rows[2][7]
    );
    assert_eq!(rows[2][9], "0");
    assert_eq!(rows[2][10], "1");
    assert_eq!(rows[2][11], "2");
    assert_eq!(rows[2][12], "0");
    // Row 0 (`CTTACTCTTCAT` -> `CTTACTNTTCAT`) is an MNV, Row 1 (`NG` -> `N`) is an indel,
    // Row 2 (`N` -> `A,C`) is a SNP.
    assert_eq!(vstats.snps, 1);
    assert_eq!(vstats.mnvs, 1);
    assert_eq!(vstats.indels, 1);
    assert_eq!(
        vstats.snps + vstats.mnvs + vstats.indels,
        vstats.sites_emitted
    );
    Ok(())
}

#[test]
#[allow(clippy::too_many_lines)]
fn nested_vcf_prevents_left_anchored_lv1_and_orphan_lv0_overlaps() -> Result<()> {
    let dir = std::env::temp_dir().join(format!(
        "horizon_resolve_lv1_overlap_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir)?;

    // Case 1: Left-anchored LV=1 deletion right at the start of an LV=0 snarl.
    // `a` (10 bp, ends in 'C') -> `d1` ("A") -> `m` ("TTT") -> `s1`/`s2` ("G"/"T") -> `z` (10 bp).
    // `LV=0` snarl `(a, z)` spans REF interior `d1+m+s1` = "ATTTG" at 1-based POS=11..15 (non-empty
    // alleles because S_ALT takes `alt_full` = "CCCCC", so NO left-anchor shift on LV=0 -> POS=11).
    // Meanwhile `S_DEL` takes `>a>m>s1>z`, creating an LV=1 bubble `(a, m)` with REF interior `d1`="A"
    // vs empty ALT -> 1-bp left-anchor shifts LV=1 to POS=10 (REF="CA", ALT="C", covering [10, 11]),
    // which overlaps LV=0 [11, 15] and sorts BEFORE POS=11! It MUST stay LV=1 linked to LV=0,
    // NOT be promoted to an overlapping LV=0!
    let gfa_left = dir.join("left_anchor.gfa");
    std::fs::write(
        &gfa_left,
        "H\tVN:Z:1.1\n\
         S\ta\tACGTACGTAC\n\
         S\td1\tA\n\
         S\tm\tTTT\n\
         S\ts1\tG\n\
         S\ts2\tT\n\
         S\talt_full\tCCCCC\n\
         S\tz\tGTGTGTGTGT\n\
         W\tREF\t0\tchr1\t0\t25\t>a>d1>m>s1>z\n\
         W\tS_DEL\t0\tchr1\t0\t24\t>a>m>s1>z\n\
         W\tS_SNP\t0\tchr1\t0\t25\t>a>d1>m>s2>z\n\
         W\tS_ALT\t0\tchr1\t0\t25\t>a>alt_full>z\n",
    )?;
    let cli = CliArgs {
        gfa_in: gfa_left.to_string_lossy().to_string(),
        out_path: None,
        vcf_path: Some("out.vcf".to_owned()),
        vcf_ref: Some("REF".to_owned()),
        nested_snarls: true,
        collapse_dups: false,
        tag_walks: false,
        max_bubble_steps: 16,
        min_anchor_bp: 100,
        divergence: None,
        threads: None,
    };
    let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
    let mut gfa_buf = Vec::new();
    let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
    let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;
    let mut vcf_buf = Vec::new();
    let _vstats = export_snarl_vcf(
        &graph,
        &cap.seg_seqs,
        &cap.walks,
        Some("REF"),
        true,
        64,
        &mut vcf_buf,
    )?;
    let text = String::from_utf8(vcf_buf)?;
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();

    // Assert NO two emitted LV=0 records overlap on REF coordinates!
    let mut lv0_spans: Vec<(usize, usize, &str)> = Vec::new();
    for r in &rows {
        if r[7].contains(";LV=0") {
            let pos: usize = r[1].parse()?;
            let end = pos + r[3].len();
            lv0_spans.push((pos, end, r[2]));
        }
    }
    assert_eq!(
        lv0_spans.len(),
        1,
        "expected exactly 1 top-level LV=0 snarl (left-anchored LV=1 must NOT promote to overlapping LV=0):\n{text}"
    );
    for w in lv0_spans.windows(2) {
        assert!(
            w[0].1 <= w[1].0,
            "overlapping LV=0 records {:?} and {:?}:\n{text}",
            w[0],
            w[1]
        );
    }
    // Verify the left-anchored deletion at POS=10 is still emitted as LV=1 linked to the LV=0 parent.
    let left_del = rows
        .iter()
        .find(|r| r[1] == "10" && r[3] == "CA" && r[4] == "C")
        .unwrap_or_else(|| panic!("missing left-anchored LV=1 deletion at POS=10:\n{text}"));
    assert!(
        left_del[7].contains(";LV=1;PS="),
        "left-anchored deletion must carry LV=1 and PS= parent ID, got:\n{text}"
    );

    // Case 2: Unemitted LV=0 (when anchor segment `z` has missing `*` sequence) produces multiple
    // overlapping orphaned LV=1 bubbles: outer `(a, b)` spanning [11, 17] and inner `(m1, m2)` at
    // POS=14. Only ONE non-overlapping outer bubble may be promoted to LV=0; the inner bubble MUST
    // attach to the promoted LV=0 via PS= instead of creating overlapping LV=0 records!
    let gfa_orphan = dir.join("orphan_overlap.gfa");
    std::fs::write(
        &gfa_orphan,
        "H\tVN:Z:1.1\n\
         S\ta\tACGTACGTAC\n\
         S\tm1\tAAA\n\
         S\ts1\tC\n\
         S\ts2\tG\n\
         S\tm2\tTTT\n\
         S\tb\tCCCCCCCCCC\n\
         S\tz\t*\tLN:i:10\n\
         W\tREF\t0\tchr1\t0\t37\t>a>m1>s1>m2>b>z\n\
         W\tH1\t0\tchr1\t0\t30\t>a>b>z\n\
         W\tH2\t0\tchr1\t0\t37\t>a>m1>s2>m2>b>z\n\
         W\tH3\t0\tchr1\t0\t20\t>a>z\n",
    )?;
    let cli2 = CliArgs {
        gfa_in: gfa_orphan.to_string_lossy().to_string(),
        out_path: None,
        vcf_path: Some("out.vcf".to_owned()),
        vcf_ref: Some("REF".to_owned()),
        nested_snarls: true,
        collapse_dups: false,
        tag_walks: false,
        max_bubble_steps: 16,
        min_anchor_bp: 100,
        divergence: None,
        threads: None,
    };
    let graph2 = pass1_index(&cli2.gfa_in, cli2.min_anchor_bp)?;
    let mut gfa_buf2 = Vec::new();
    let mut cap2 = resolve::VcfCapture::with_segments(graph2.num_segs());
    let _ = pass2_rewrite(&cli2.gfa_in, &mut gfa_buf2, &graph2, &cli2, Some(&mut cap2))?;
    let mut vcf_buf2 = Vec::new();
    let _vstats2 = export_snarl_vcf(
        &graph2,
        &cap2.seg_seqs,
        &cap2.walks,
        Some("REF"),
        true,
        64,
        &mut vcf_buf2,
    )?;
    let _ = std::fs::remove_dir_all(&dir);
    let text2 = String::from_utf8(vcf_buf2)?;
    let rows2: Vec<Vec<&str>> = text2
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    let mut lv0_spans2: Vec<(usize, usize, &str)> = Vec::new();
    for r in &rows2 {
        if r[7].contains(";LV=0") {
            let pos: usize = r[1].parse()?;
            lv0_spans2.push((pos, pos + r[3].len(), r[2]));
        }
    }
    assert_eq!(
        lv0_spans2.len(),
        1,
        "only outer orphan bubble may promote to LV=0; inner overlapping bubble must stay LV=1:\n{text2}"
    );
    assert!(
        rows2.iter().any(|r| r[7].contains(";LV=1;PS=")),
        "inner overlapping orphan bubble must attach to promoted LV=0 via PS=:\n{text2}"
    );
    Ok(())
}

#[test]
#[allow(clippy::too_many_lines)]
fn vcf_export_annotates_inversion_alleles_with_svtype_end_svlen() -> Result<()> {
    let dir = std::env::temp_dir().join(format!(
        "horizon_resolve_svtype_inv_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir)?;

    let inv_block =
        "ACGTCTAGGAATTCGGATCCATCGATCGTTAGGCCTAGCATGCTAGCTAACGTTAGGCCTTAGCATGCTAGCTAGGCCTA"; // 80 bp, non-palindromic
    let inv_block_rc: String = inv_block
        .bytes()
        .rev()
        .map(|b| horizon_aligner::dna::complement_base(b) as char)
        .collect();
    // Mutate 2 bases in the distinct-segment RC copy (97.5% RC identity) to test sequence-based INV detection.
    let mut rc_mut_bytes = inv_block_rc.as_bytes().to_vec();
    rc_mut_bytes[10] = if rc_mut_bytes[10] == b'A' { b'C' } else { b'A' };
    rc_mut_bytes[55] = if rc_mut_bytes[55] == b'G' { b'T' } else { b'G' };
    let rc_mut = String::from_utf8(rc_mut_bytes)?;
    let unrelated =
        "TTTTTTTTCCCCCCCCGGGGGGGGAAAAAAAATTTTTTTTCCCCCCCCGGGGGGGGAAAAAAAATTTTTTTTCCCCCCCC"; // 80 bp unrelated MNV

    let gfa_in = dir.join("inv.gfa");
    std::fs::write(
        &gfa_in,
        format!(
            "H\tVN:Z:1.1\n\
             S\ta\tACGTACGTAC\n\
             S\tinv1\t{inv_block}\n\
             S\tb\tGTGTGTGTGT\n\
             S\tinv2\t{inv_block}\n\
             S\tinv2_rc\t{rc_mut}\n\
             S\tc\tACACACACAC\n\
             S\tmnv_r\t{inv_block}\n\
             S\tmnv_a\t{unrelated}\n\
             S\td\tTATATATATA\n\
             S\tp1\tAA\n\
             S\tp2\tCC\n\
             S\tp3\tG\n\
             S\tinv3\t{inv_block}\n\
             S\te\tCGCGCGCGCG\n\
             W\tREF\t0\tchr1\t0\t374\t>a>inv1>b>inv2>c>mnv_r>d>p1>inv3>p2>e\n\
             W\tHAP\t0\tchr1\t0\t371\t>a<inv1>b>inv2_rc>c>mnv_a>d>p3<inv3>e\n"
        ),
    )?;

    let cli = CliArgs {
        gfa_in: gfa_in.to_string_lossy().to_string(),
        out_path: None,
        vcf_path: Some("out.vcf".to_owned()),
        vcf_ref: Some("REF".to_owned()),
        nested_snarls: false,
        collapse_dups: false,
        tag_walks: false,
        max_bubble_steps: 16,
        min_anchor_bp: 100,
        divergence: None,
        threads: None,
    };
    let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
    let mut gfa_buf = Vec::new();
    let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
    let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;
    let mut vcf_buf = Vec::new();
    let _vstats = export_snarl_vcf(
        &graph,
        &cap.seg_seqs,
        &cap.walks,
        Some("REF"),
        false,
        64,
        &mut vcf_buf,
    )?;
    let _ = std::fs::remove_dir_all(&dir);
    let text = String::from_utf8(vcf_buf)?;
    assert!(
        text.contains("##INFO=<ID=SVTYPE,"),
        "missing SVTYPE header:\n{text}"
    );
    assert!(
        text.contains("##INFO=<ID=END,"),
        "missing END header:\n{text}"
    );
    assert!(
        text.contains("##INFO=<ID=SVLEN,"),
        "missing SVLEN header:\n{text}"
    );

    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 4, "expected 4 snarl records:\n{text}");
    // Row 0: graph-oriented `<inv1` at POS=11, len=80 -> END=90, SVLEN=80
    assert!(
        rows[0][7].contains("SVTYPE=INV;END=90;SVLEN=80"),
        "row 0 missing SVTYPE=INV;END=90;SVLEN=80:\n{}",
        rows[0][7]
    );
    // Row 1: distinct-segment 97.5% RC `inv2_rc` at POS=101, len=80 -> END=180, SVLEN=80
    assert!(
        rows[1][7].contains("SVTYPE=INV;END=180;SVLEN=80"),
        "row 1 missing SVTYPE=INV;END=180;SVLEN=80:\n{}",
        rows[1][7]
    );
    // Row 2: unrelated 80 bp MNV must NOT be labeled SVTYPE=INV
    assert!(
        !rows[2][7].contains("SVTYPE=INV"),
        "row 2 non-inversion MNV falsely labeled SVTYPE=INV:\n{}",
        rows[2][7]
    );
    // Row 3: breakpoint-padded graph inversion (`>p1>inv3>p2` vs `>p3<inv3`, REF len=84) MUST be SVTYPE=INV
    assert!(
        rows[3][7].contains("SVTYPE=INV;END=364;SVLEN=84"),
        "row 3 breakpoint-padded graph inversion missing SVTYPE=INV:\n{}",
        rows[3][7]
    );
    Ok(())
}
