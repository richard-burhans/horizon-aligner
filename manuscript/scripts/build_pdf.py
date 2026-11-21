#!/usr/bin/env python3
"""Render manuscript/horizon_paper.md to manuscript/horizon_paper.pdf using markdown + weasyprint."""

from __future__ import annotations

import re
from pathlib import Path
import markdown
import weasyprint


CSS = """
@page {
  size: letter;
  margin: 0.78in 0.82in 0.82in 0.82in;
  @bottom-center {
    content: counter(page);
    font-family: "Liberation Serif", "Times New Roman", Georgia, serif;
    font-size: 9pt;
    color: #555;
  }
}
@page :first {
  @bottom-center { content: normal; }
}
body {
  font-family: "Liberation Serif", "Times New Roman", Georgia, serif;
  font-size: 10.5pt;
  line-height: 1.44;
  color: #1a1a1a;
  text-align: justify;
  hyphens: auto;
}
header.title-block {
  text-align: left;
  border-bottom: 1.2pt solid #222;
  padding-bottom: 10pt;
  margin-bottom: 14pt;
}
header.title-block h1.doc-title {
  font-family: "Liberation Sans", "Helvetica Neue", Arial, sans-serif;
  font-size: 17pt;
  line-height: 1.2;
  font-weight: 700;
  margin: 0 0 6pt 0;
  color: #111;
}
header.title-block .doc-meta {
  font-family: "Liberation Sans", "Helvetica Neue", Arial, sans-serif;
  font-size: 9.5pt;
  color: #333;
  margin: 0;
}
h1, h2, h3 {
  font-family: "Liberation Sans", "Helvetica Neue", Arial, sans-serif;
  color: #111;
  text-align: left;
  page-break-after: avoid;
  break-after: avoid;
}
h1 {
  font-size: 12.5pt;
  font-weight: 700;
  margin: 16pt 0 6pt 0;
  border-bottom: 0.5pt solid #ccc;
  padding-bottom: 2pt;
}
h2 {
  font-size: 10.8pt;
  font-weight: 700;
  margin: 12pt 0 4pt 0;
}
p {
  margin: 0 0 7pt 0;
}
code {
  font-family: "DejaVu Sans Mono", "Liberation Mono", monospace;
  font-size: 8.8pt;
  background-color: #f4f5f7;
  padding: 0.5pt 2.5pt;
  border-radius: 2pt;
}
img {
  display: block;
  max-width: 100%;
  height: auto;
  margin: 10pt auto 4pt auto;
  page-break-inside: avoid;
  break-inside: avoid;
}
p:has(img) {
  text-align: center;
  margin: 8pt 0 3pt 0;
  page-break-after: avoid;
  break-after: avoid;
}
table {
  width: 100%;
  border-collapse: collapse;
  margin: 10pt 0 12pt 0;
  font-size: 8.4pt;
  line-height: 1.3;
  text-align: left;
  page-break-inside: avoid;
  break-inside: avoid;
}
th, td {
  border-bottom: 0.5pt solid #d0d4da;
  padding: 4pt 5pt;
  vertical-align: top;
}
th {
  border-top: 1pt solid #222;
  border-bottom: 1pt solid #222;
  font-family: "Liberation Sans", "Helvetica Neue", Arial, sans-serif;
  font-weight: 700;
  background-color: #f7f8fa;
}
ul, ol {
  margin: 4pt 0 8pt 0;
  padding-left: 18pt;
}
li {
  margin-bottom: 3pt;
  text-align: left;
}
"""


def parse_frontmatter(text: str) -> tuple[dict[str, str], str]:
    m = re.match(r"^---\n(.*?)\n---\n(.*)$", text, flags=re.DOTALL)
    if not m:
        return {}, text
    fm_raw, body = m.group(1), m.group(2)
    meta: dict[str, str] = {}
    title_m = re.search(r'^title:\s*"?(.*?)"?\s*$', fm_raw, flags=re.MULTILINE)
    if title_m:
        meta["title"] = title_m.group(1)
    date_m = re.search(r'^date:\s*"?(.*?)"?\s*$', fm_raw, flags=re.MULTILINE)
    if date_m:
        meta["date"] = date_m.group(1)
    authors = re.findall(r'^\s*-\s*"?(.*?)"?\s*$', fm_raw, flags=re.MULTILINE)
    if authors:
        meta["author"] = ", ".join(authors)
    return meta, body


def render_md_to_pdf(md_path: Path, pdf_path: Path, manuscript_dir: Path) -> None:
    raw = md_path.read_text(encoding="utf-8")
    meta, body = parse_frontmatter(raw)
    # Strip pandoc-style {.unnumbered} attributes
    body = re.sub(r"\s*\{\.unnumbered\}", "", body)

    html_body = markdown.markdown(
        body,
        extensions=["tables", "fenced_code", "sane_lists"],
    )
    title = meta.get("title", "Horizon")
    author = meta.get("author", "")
    date = meta.get("date", "")
    meta_line = " | ".join(x for x in [author, date] if x)

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>{CSS}</style>
</head>
<body>
<header class="title-block">
  <h1 class="doc-title">{title}</h1>
  <p class="doc-meta">{meta_line}</p>
</header>
{html_body}
</body>
</html>
"""
    weasyprint.HTML(string=html_doc, base_url=str(manuscript_dir)).write_pdf(str(pdf_path))
    print(f"Wrote {pdf_path} ({pdf_path.stat().st_size} bytes)")


def main() -> None:
    manuscript_dir = Path(__file__).resolve().parents[1]
    for stem in ("horizon_paper", "horizon_companion"):
        render_md_to_pdf(
            manuscript_dir / f"{stem}.md",
            manuscript_dir / f"{stem}.pdf",
            manuscript_dir,
        )


if __name__ == "__main__":
    main()
