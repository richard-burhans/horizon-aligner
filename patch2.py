with open("src/formats.rs", "r") as f:
    lines = f.readlines()

new_lines = []
skip = False
for line in lines:
    if "let (y_start, y_span, y_strand) = if is_rev {" in line:
        skip = True
        new_lines.append("        let y_lo = if is_rev { blk[blk.len() - 1].1 } else { blk[0].1 };\n")
        new_lines.append("        let y_hi = if is_rev { blk[0].1 + 1 } else { blk[blk.len() - 1].1 + 1 };\n")
        new_lines.append("        let y_iv = crate::coord::Interval::new(y_lo, y_hi, if is_rev { '-' } else { '+' }, len_y);\n")
        new_lines.append("\n")
        new_lines.append("        let _ = writeln!(out, \"a score={}.0\", blk.len());\n")
        new_lines.append("        let _ = writeln!(out, \"s {name_x} {x0} {} + {len_x} {ax}\", x_end - x0);\n")
        new_lines.append("        let _ = writeln!(\n")
        new_lines.append("            out,\n")
        new_lines.append("            \"s {name_y} {} {} {} {len_y} {ay}\\n\",\n")
        new_lines.append("            y_iv.ucsc_start(), y_iv.span(), y_iv.strand\n")
        new_lines.append("        );\n")
    elif skip and "        );" in line:
        skip = False
    elif not skip:
        new_lines.append(line)

with open("src/formats.rs", "w") as f:
    f.writelines(new_lines)
