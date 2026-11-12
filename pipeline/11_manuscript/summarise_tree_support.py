#!/usr/bin/env python3
"""Summarise branch support for the species tree (results/tree_support): SH-aLRT and UFBoot from the registered LG+G4 analysis
and the partitioned ModelFinder analysis, gene concordance (gCF, 725 gene trees) and likelihood site concordance (sCFL) for each
branch, named by the clade it defines (tree rooted on A. thaliana). Reads IQ-TREE outputs only.
Usage (repository root): python pipeline/11_manuscript/summarise_tree_support.py
"""
import re
import pandas as pd

R = "results/tree_support"
NAMES = {
    frozenset({"Vrad", "Vung"}): "Vigna",
    frozenset({"Pvul", "Vrad", "Vung"}): "Phaseolus + Vigna",
    frozenset({"Gmax", "Pvul", "Vrad", "Vung"}): "Glycine + Phaseolus + Vigna",
    frozenset({"Ccaj", "Gmax", "Pvul", "Vrad", "Vung"}): "Millettioids",
    frozenset({"Mtru", "Tpra"}): "Medicago + Trifolium",
    frozenset({"Mtru", "Tpra", "Psat"}): "Medicago + Trifolium + Pisum",
    frozenset({"Mtru", "Tpra", "Psat", "Cari"}): "IRLC",
    frozenset({"Mtru", "Tpra", "Psat", "Cari", "Ljap"}): "Robinioid + IRLC",
    frozenset({"Mtru", "Tpra", "Psat", "Cari", "Ljap", "Ccaj", "Gmax", "Pvul", "Vrad", "Vung"}): "Millettioids + Robinioid + IRLC (excludes Lupinus)",
    frozenset({"Mtru", "Tpra", "Psat", "Cari", "Ljap", "Ccaj", "Gmax", "Pvul", "Vrad", "Vung", "Lang"}): "all legumes except Arachis",
}


def clades(newick):
    """Return {clade (frozenset of tips): label} for every internal node of a Newick string."""
    out, stack, pos = {}, [], 0
    tokens = re.findall(r"\(|\)[^,():;]*(?::[0-9.eE+-]+)?|[^,():;]+(?::[0-9.eE+-]+)?|,|;", newick)
    for t in tokens:
        if t == "(":
            stack.append(set())
        elif t.startswith(")"):
            members = stack.pop()
            label = t[1:].split(":")[0]
            out[frozenset(members)] = label
            if stack:
                stack[-1] |= members
        elif t not in (",", ";"):
            name = t.split(":")[0]
            if stack:
                stack[-1].add(name)
    return out


rows = []
lab_lg = clades(open(f"{R}/lg_boot.treefile").read())
lab_part = clades(open(f"{R}/part.treefile").read())
lab_g = clades(open(f"{R}/concord_gcf.cf.tree").read())
lab_s = clades(open(f"{R}/concord_scfl.cf.tree").read())
gstat = pd.read_csv(f"{R}/concord_gcf.cf.stat", sep="\t", comment="#")
sstat = pd.read_csv(f"{R}/concord_scfl.cf.stat", sep="\t", comment="#")
for c, name in NAMES.items():
    alrt, ufb = lab_lg[c].split("/")
    palrt, pufb = lab_part[c].split("/")
    gcf = float(lab_g[c].split("/")[2]); scf = float(lab_s[c].split("/")[2])
    g = gstat.iloc[(gstat.gCF - gcf).abs().argmin()]; s = sstat.iloc[(sstat.sCF - scf).abs().argmin()]
    rows.append(dict(clade=name, tips=",".join(sorted(c)), shalrt_LG=float(alrt), ufboot_LG=float(ufb), shalrt_partitioned=float(palrt), ufboot_partitioned=float(pufb),
                     gCF=gcf, gDF1=g.gDF1, gDF2=g.gDF2, gDFP=g.gDFP, gN=int(g.gN), sCF=scf, sDF1=s.sDF1, sDF2=s.sDF2))
T = pd.DataFrame(rows)
T["concordant_beats_alternatives"] = (T.gCF > T[["gDF1", "gDF2"]].max(axis=1)) & (T.sCF > T[["sDF1", "sDF2"]].max(axis=1))
T.to_csv(f"{R}/support_summary.tsv", sep="\t", index=False)
strip = lambda s: re.sub(r"\)[^,():;]*", ")", re.sub(r":[0-9.eE+-]+", "", s))
same = strip(open(f"{R}/lg_boot.treefile").read()) == strip(open(f"{R}/part.treefile").read()) == strip(open("results/server_logs/concat.fa.treefile").read())
pd.DataFrame([dict(topology_identical_registered_LG_partitioned=same, n_branches=len(T))]).to_csv(f"{R}/support_topology.tsv", sep="\t", index=False)
print(T.to_string(index=False)); print("identical topology (registered, LG+bootstrap, partitioned):", same)
