#!/usr/bin/env python3
"""Place each main figure (image + legend) and each main table (caption + table) directly after the paragraph that first
cites it, and drop the duplicated end-of-document legends and tables. Supplementary legends and the supplementary-table
list stay at the end. Reads the rendered manuscript; writes a Markdown file for md_to_docx.py (image lines: ![](path)).
Usage (repository root): python pipeline/11_manuscript/embed_figures_tables.py manuscript_mpe/manuscript.md manuscript_mpe/manuscript_embedded.md
"""
import re, sys

src, dst = sys.argv[1], sys.argv[2]
s = open(src, encoding="utf-8").read()
i_leg, i_supp, i_tab = s.index("# Figure legends"), s.index("# Supplementary tables"), s.index("# Tables")
main, legends, supp, tables = s[:i_leg], s[i_leg:i_supp], s[i_supp:i_tab], s[i_tab:]

# main-figure legends (Fig. 1-6) and supplementary-figure legends
paras = [p for p in legends.split("\n\n") if p.strip() and not p.startswith("# ")]
fig_leg = {int(re.match(r"\*\*Fig\. (\d+)\.", p).group(1)): p for p in paras if re.match(r"\*\*Fig\. \d+\.", p)}
supp_leg = [p for p in paras if re.match(r"\*\*Fig\. S\d+\.", p)]
# tables: caption paragraph followed by its pipe table
tab_blocks = {}
for m in re.finditer(r"(\*\*Table (\d+)\.\*\*[^\n]*)\n\n((?:\|[^\n]*\n?)+)", tables):
    tab_blocks[int(m.group(2))] = m.group(1).strip() + "\n\n" + m.group(3).strip()
assert sorted(fig_leg) == [1, 2, 3, 4, 5, 6] and sorted(tab_blocks) == [1, 2], (sorted(fig_leg), sorted(tab_blocks))

lines = main.split("\n")
start = next(k for k, l in enumerate(lines) if l.startswith("# 1. Introduction"))
cite_fig = lambda n: re.compile(r"\bFigs?\.\s+(?:[0-9a-z,–\- ]*?(?:and\s+)?)?" + str(n) + r"[a-z]?\b")
cite_tab = lambda n: re.compile(r"\bTables?\s+(?:[0-9 ,]*?(?:and\s+)?)?" + str(n) + r"\b(?!\.\d)")
inserts = {}                                        # line index -> blocks to insert after it, in citation order
for kind, items, pat in (("fig", fig_leg, cite_fig), ("tab", tab_blocks, cite_tab)):
    for n in sorted(items):
        k = next(k for k in range(start, len(lines)) if lines[k].strip() and not lines[k].startswith("#") and pat(n).search(lines[k]))
        block = (f"![](manuscript_mpe/figures/Fig{n}.png)\n\n" + fig_leg[n]) if kind == "fig" else tab_blocks[n]
        pos = pat(n).search(lines[k]).start()
        inserts.setdefault(k, []).append((pos, kind, n, block))
out = []
for k, l in enumerate(lines):
    out.append(l)
    for _, kind, n, block in sorted(inserts.get(k, [])):
        out += ["", block, ""]
        print(f"{'Fig.' if kind == 'fig' else 'Table'} {n} placed after line {k + 1}: {l[:70]}...")
text = "\n".join(out).rstrip() + "\n\n# Supplementary figure legends\n\n" + "\n\n".join(supp_leg) + "\n\n" + supp.strip() + "\n"
open(dst, "w", encoding="utf-8", newline="\n").write(text)
print("wrote", dst)
