import sys

with open("src/bin/horizon_eval.rs", "r") as f:
    text = f.read()

# Replace hardcoded main
old_main = """fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        anyhow::bail!("usage: eval_pf_pv <pf.fa> <pv.fa> <rbest.maf>");
    }
    let pf = read_named_contig(Path::new(&args[1]), "CP131981.1")?;
    let pv = read_named_contig(Path::new(&args[2]), "CAJZCX010000015.1")?;
    let maf_all = parse_maf_pairs(Path::new(&args[3]), "CP131981.1", "CAJZCX010000015.1")?;

    let t0 = Instant::now();
    let aln = align_contig_pair(&pf, &pv)?;
    let elapsed = t0.elapsed().as_secs_f64();"""

new_main = """fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 6 {
        anyhow::bail!("usage: horizon_eval <ref_name> <target_name> <ref.fa> <target.fa> <truth.maf>");
    }
    let ref_name = &args[1];
    let target_name = &args[2];
    let pf = read_named_contig(Path::new(&args[3]), ref_name)?;
    let pv = read_named_contig(Path::new(&args[4]), target_name)?;
    let maf_all = parse_maf_pairs(Path::new(&args[5]), ref_name, target_name)?;

    let t0 = Instant::now();
    let aln = align_contig_pair(&pf, &pv)?;
    let elapsed = t0.elapsed().as_secs_f64();"""

text = text.replace(old_main, new_main)
text = text.replace("usage: eval_pf_pv", "usage: horizon_eval")
text = text.replace("eval_pf_pv", "horizon_eval")

with open("src/bin/horizon_eval.rs", "w") as f:
    f.write(text)
