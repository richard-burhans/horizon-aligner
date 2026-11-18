with open("src/formats.rs", "r") as f:
    text = f.read()

target = """        let (y_start, y_span, y_strand) = if is_rev {
            let y_lo = blk[blk.len() - 1].1;
            let y_hi = blk[0].1 + 1;
            (len_y - y_hi, y_hi - y_lo, '-')
        } else {
            let y_lo = blk[0].1;
            let y_hi = blk[blk.len() - 1].1 + 1;
            (y_lo, y_hi - y_lo, '+')
        };

        let _ = writeln!(out, "a score={}.0", blk.len());
        let _ = writeln!(out, "s {name_x} {x0} {} + {len_x} {ax}", x_end - x0);
        let _ = writeln!(
            out,
            "s {name_y} {y_start} {y_span} {y_strand} {len_y} {ay}\n"
        );"""

replacement = """        let y_lo = if is_rev { blk[blk.len() - 1].1 } else { blk[0].1 };
        let y_hi = if is_rev { blk[0].1 + 1 } else { blk[blk.len() - 1].1 + 1 };
        let y_iv = crate::coord::Interval::new(y_lo, y_hi, if is_rev { '-' } else { '+' }, len_y);

        let _ = writeln!(out, "a score={}.0", blk.len());
        let _ = writeln!(out, "s {name_x} {x0} {} + {len_x} {ax}", x_end - x0);
        let _ = writeln!(
            out,
            "s {name_y} {} {} {} {len_y} {ay}\n",
            y_iv.ucsc_start(), y_iv.span(), y_iv.strand
        );"""

if target in text:
    with open("src/formats.rs", "w") as f:
        f.write(text.replace(target, replacement))
    print("Patched formats.rs successfully.")
else:
    print("Could not find target in formats.rs.")
