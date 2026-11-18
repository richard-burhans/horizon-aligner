with open("src/egraph.rs", "r") as f:
    lines = f.readlines()

new_lines = []
in_fn = False
fn_lines = []
for i, line in enumerate(lines):
    if "/// Resolve LCB" in line or "#[must_use]" in line or "pub fn resolve_lcb_duplications" in line:
        fn_lines.append(line)
        in_fn = True
    elif in_fn:
        fn_lines.append(line)
        if line.startswith("}"):
            in_fn = False
    else:
        new_lines.append(line)

# find where #[cfg(test)] starts
test_start = -1
for i, line in enumerate(new_lines):
    if "#[cfg(test)]" in line:
        test_start = i
        break

if test_start != -1 and fn_lines:
    new_lines.insert(test_start, "".join(fn_lines) + "\n")
    with open("src/egraph.rs", "w") as f:
        f.write("".join(new_lines))
    print("Fixed egraph2.rs")
