use super::*;

#[test]
fn mdl_gate_and_hmm_estimation_smoke() {
    let x = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let y = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let lcb: Vec<(usize, usize)> = (0..10).map(|i| (i * 4, i * 4)).collect();
    assert!(lcb_passes_mdl_gate(&lcb, x, y, y.len(), &|_| 4, 8.0));
    assert!(!lcb_passes_mdl_gate(
        &lcb,
        x,
        b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        40,
        &|_| 4,
        8.0
    ));

    let empty = std::collections::HashMap::new();
    let mut spans_fwd = empty.clone();
    let mut spans_rc = empty;
    let cfg = AlignConfig::default();
    let lcbs_fwd = vec![lcb];
    let y_rc = LazyRc::new(y);
    let ctx = AlignCtx {
        seq_x: x,
        seq_y: y,
        seq_y_rc: &y_rc,
        lcbs_fwd: &lcbs_fwd,
        lcbs_rc: &[],
        spans_fwd: &mut spans_fwd,
        spans_rc: &mut spans_rc,
        span: 4,
        effective_min_bits: 8.0,
        lcb_gate_bits: 8.0,
        cfg: &cfg,
    };
    let params = estimate_hmm_params(&ctx);
    assert!((0.60..=0.72).contains(&params.identity_prior));
}

#[test]
fn xdrop_respects_effective_min_bits_and_anchor_emission_skips_ns() {
    // Construct a 120-bp sequence where:
    // - [0..20) is a 20-bp matching flank (~28 bits: passes hardcoded 10.0, fails min_delta_bits=60.0)
    // - [20..40) is a 20-bp anchor with 'N' at positions 25..28 in BOTH x and y
    // - [40..120) is non-matching filler
    let mut x = b"AGTCATGCTAGTCGAACTGA".to_vec(); // 0..20 non-palindromic matching flank
    x.extend_from_slice(b"GCTAGNNNCTAGGCTAGCTA"); // 20..40 anchor with NNN
    x.extend_from_slice(&[b'A'; 40]); // 40..80
    let mut y = b"AGTCATGCTAGTCGAACTGA".to_vec();
    y.extend_from_slice(b"GCTAGNNNCTAGGCTAGCTA");
    y.extend_from_slice(&[b'C'; 40]);

    let lcb = vec![(20, 20)];
    let mut spans_fwd = std::collections::HashMap::new();
    spans_fwd.insert((20, 20), 20);
    let mut spans_rc = std::collections::HashMap::new();
    let cfg = AlignConfig {
        min_delta_bits: 60.0,
        max_ext: 40,
        ..AlignConfig::default()
    };
    let y_rc = LazyRc::new(&y);
    let lcbs_fwd = vec![lcb];
    let mut ctx = AlignCtx {
        seq_x: &x,
        seq_y: &y,
        seq_y_rc: &y_rc,
        lcbs_fwd: &lcbs_fwd,
        lcbs_rc: &[],
        spans_fwd: &mut spans_fwd,
        spans_rc: &mut spans_rc,
        span: 20,
        effective_min_bits: 60.0,
        lcb_gate_bits: 10.0,
        cfg: &cfg,
    };
    let ext = schedule_and_run_extension(&mut ctx).expect("extend");
    let pairs = runs_to_pairs(&ext.runs);
    // Anchor emission must NOT emit N-N pairs at 25..28.
    assert!(
        pairs.iter().all(|&(px, py)| x[px] != b'N' && y[py] != b'N'),
        "anchor emission must skip N-N wildcard positions, got {pairs:?}"
    );
    // X-drop flank [0..20) (~28 bits) must NOT be emitted when effective_min_bits = 60.0.
    assert!(
        pairs.iter().all(|&(px, _)| px >= 20),
        "xdrop flank bypassed effective_min_bits=60.0, got {pairs:?}"
    );
}

#[test]
fn lazy_rc_is_not_materialized_until_needed() {
    let seq = b"ACGTACGT";
    let lazy = LazyRc::new(seq);
    assert!(lazy.rc.get().is_none());
    assert_eq!(lazy.get(), b"ACGTACGT");
    assert!(lazy.rc.get().is_some());
}
