#!/usr/bin/env python3
"""Fill {{key}} placeholders in the manuscript template with values from the numbers ledger.

Every number in the rendered text that depends on data comes from results/manuscript_numbers/ledger.tsv.
Unknown keys stop the run. Leading hyphen-minus signs in filled values become typographic minus signs.
Also lists every digit sequence typed directly in the template (outside placeholders), so each can be
checked by hand as a design constant, a citation year or a section number.
Usage (repository root): python pipeline/11_manuscript/render_manuscript.py manuscript_mpe/manuscript.md
  (concatenates manuscript_mpe/sections/*.md in name order)
"""
import re, sys
import pandas as pd

import glob
dst = sys.argv[1]
parts = sorted(glob.glob("manuscript_mpe/sections/*.md"))
L = pd.read_csv("results/manuscript_numbers/ledger.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("id")["value"]
t = "\n\n".join(open(f, encoding="utf-8").read().strip() for f in parts) + "\n"
# ---- citations: [@key; @key2] -> (Author, Year; Author et al., Year); reference list at {{REFERENCES}}
R = pd.read_csv("manuscript_mpe/references_verified.tsv", sep="\t", dtype=str, encoding="utf-8", keep_default_na=False).set_index("key")
bad = R[R.status != "VERIFIED"]
if len(bad):
    sys.exit("unverified references in table: " + ", ".join(bad.index))
cited = []
def surnames(a):
    return [x.split(", ")[0] for x in re.split(r"(?<=\.), ", a) if x]
def cite_one(k):
    if k not in R.index:
        sys.exit(f"citation key not in references_verified.tsv: {k}")
    if k not in cited:
        cited.append(k)
    sn, y = surnames(R.at[k, "authors_elsevier"]), R.at[k, "year"]
    who = sn[0] if len(sn) == 1 else f"{sn[0]} and {sn[1]}" if len(sn) == 2 else f"{sn[0]} et al."
    return f"{who}, {y}"
t = re.sub(r"\[(@[A-Za-z0-9]+(?:; ?@[A-Za-z0-9]+)*)\]",
           lambda m: "(" + "; ".join(cite_one(k.strip()[1:]) for k in m.group(1).split(";")) + ")", t)
def dotted(iso, full):
    # MPE style: a word of the NLM abbreviation takes a period unless it is a whole word of the full title
    words = set(re.findall(r"[a-z]+", full.lower()))
    return " ".join(w if w.lower() in words else w + "." for w in iso.split())


def ref_entry(k):
    r = R.loc[k]; title = r.title.rstrip(".")
    end = "" if title.endswith(("?", "!")) else "."
    pages = re.sub(r"(?<=\d)-(?=[A-Za-z]?\d)", "–", str(r.pages_or_article))
    loc = f"{dotted(r.journal_iso, r.journal_full)} {r.volume}, {pages}".strip().rstrip(",")
    return f"{r.authors_elsevier.rstrip('.')}., {r.year}. {title}{end} {loc}. https://doi.org/{r.doi}"
def ref_sort(k):
    # Elsevier name-year order: first author (surname, then initials); one author, then two (by second author), then et al.; then year
    auth = [x for x in re.split(r"(?<=\.), ", R.at[k, "authors_elsevier"]) if x]
    first = auth[0].lower()
    n = min(len(auth), 3)
    return (first.split(", ")[0], first, n, auth[1].lower() if n == 2 else "", R.at[k, "year"])
order = sorted(cited, key=ref_sort)
t = t.replace("{{REFERENCES}}", "\n\n".join(ref_entry(k) for k in order))
keys = re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", t)
missing = sorted({k for k in keys if k not in L.index})
if missing:
    sys.exit("unknown ledger keys: " + ", ".join(missing))


def fill(m):
    v = L[m.group(1)]
    if re.fullmatch(r"\d{5,}", v):              # thousands separators for counts (Elsevier style)
        v = f"{int(v):,}"
    return re.sub(r"^-(?=\d)", "−", v)


out = re.sub(r"\{\{([A-Za-z0-9_]+)\}\}", fill, t)
open(dst, "w", encoding="utf-8", newline="\n").write(out)
body = t.split("# References")[0]
typed = sorted({m for m in re.findall(r"(?<![A-Za-z_{])\d[\d.,]*", re.sub(r"\{\{[A-Za-z0-9_]+\}\}|\([^()]*\d{4}[^()]*\)", "", body))})
print(f"filled {len(keys)} placeholders ({len(set(keys))} distinct keys); {len(cited)} references cited -> {dst}")
print("numbers typed directly in the template (check each):", " ".join(typed))
