with open("src/egraph.rs", "r") as f:
    content = f.read()

parts = content.split("mod tests {")
if len(parts) == 2:
    code = parts[0]
    tests = "mod tests {" + parts[1]
    
    # We want to extract the resolve_lcb_duplications that we added at the end of the file
    lines = tests.splitlines()
    fn_start = -1
    for i, line in enumerate(lines):
        if "pub fn resolve_lcb_duplications" in line:
            fn_start = i - 1  # include the docstring
            break
            
    if fn_start != -1:
        fn_code = "\n".join(lines[fn_start:])
        tests = "\n".join(lines[:fn_start])
        
        # Add #[must_use]
        fn_code = fn_code.replace("pub fn resolve_lcb_duplications", "#[must_use]\npub fn resolve_lcb_duplications")
        
        new_content = code + "\n" + fn_code + "\n\n" + tests
        
        with open("src/egraph.rs", "w") as f:
            f.write(new_content)
        print("Fixed egraph.rs")
