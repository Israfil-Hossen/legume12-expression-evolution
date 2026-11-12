#!/usr/bin/env python3
"""EXPLORATORY (post hoc, prompted by a pre-submission review; not in any specification).

R1  Relative dosage (Step 5 Part 3) by expression level: log2(x + 1) compresses a doubling at low expression,
    which would bias the registered statistic towards 'partial'. Reported per bin of reference expression
    (mean log2 scaled TPM of the single-copy species), with calibration orthogroups from the same bin, and next to
    the value full dosage would give under the same pseudocount.
R2  Absolute-scale copy-number gradient: median divergence of multi-copy minus single-copy orthogroups
    (log2 units) and the ratio of the two medians, within studies and between species, the between-species excess
    over the two species' within-study values, and its change with divergence time (55 pairs without C. cajan,
    species bootstrap with the Step 2 rules).

Expression loading and the dosage definition are copied from step5_analysis.py; the script first reproduces the
registered dosage estimates and stops if they differ.
Usage (repository root): python pipeline/12_step5/explore_review_checks.py <DATA_dir>
"""
import os, sys, itertools
import numpy as np
import pandas as pd

DATA = sys.argv[1]
RES, OUT = "results", "results/exploratory_review"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261009)
NBOOT = 2000
SPECIES = ["Gmax", "Pvul", "Vrad", "Vung", "Ccaj", "Mtru", "Cari", "Psat", "Tpra", "Ljap", "Lang", "Adur"]
OFR = os.path.join(DATA, "proteomes", "OrthoFinder", "Results_Sep23")
LOG = open(os.path.join(OUT, "report_review_checks.txt"), "w", encoding="utf-8")


def say(s=""):
    print(s); LOG.write(s + "\n"); LOG.flush()


# ------------------------------------------------------------------ expression (as step5_analysis.py)
qc1 = pd.read_csv(f"{RES}/results_step2/qc_samples.tsv", sep="\t"); qc1 = qc1[qc1.qc_pass]
qc2 = pd.read_csv(f"{RES}/results_step4/qc_samples_step4.tsv", sep="\t"); qc2 = qc2[qc2.qc_pass]
og = pd.read_csv(os.path.join(OFR, "Orthogroups", "Orthogroups.tsv"), sep="\t", dtype=str).fillna("")
gene2og, cn = {}, {}
for sp in SPECIES:
    g2o = {}
    for o, v in zip(og.Orthogroup, og[sp]):
        if v:
            gs = [g[len(sp) + 1:] for g in v.split(", ")]
            for g in gs:
                g2o[g] = o
            cn[(o, sp)] = len(gs)
    gene2og[sp] = g2o
CN = pd.Series(cn).unstack(fill_value=0).reindex(columns=SPECIES, fill_value=0)
SC = CN.index[(CN == 1).all(axis=1)]
tx2g = {sp: pd.read_csv(os.path.join(DATA, "tables", f"{sp}_tx2gene.tsv"), sep="\t", header=None, names=["Name", "gene"]) for sp in SPECIES}
OGT = {}
for b, sp in list(zip(qc1.biorep, qc1.species)) + list(zip(qc2.biorep, qc2.species)):
    q = pd.read_csv(os.path.join(DATA, "quant", b, "quant.sf"), sep="\t", usecols=["Name", "TPM"])
    gt = q.merge(tx2g[sp], on="Name").groupby("gene").TPM.sum()
    gt = gt[gt.index.map(lambda g: g in gene2og[sp])]
    OGT[b] = gt.groupby(gt.index.map(gene2og[sp])).sum().reindex(CN.index[CN[sp] > 0], fill_value=0.0)
sf = {b: np.median(OGT[b].reindex(SC).values) for b in OGT}
C = np.median([sf[b] for b in qc1.biorep])
QC1 = {sp: list(qc1[qc1.species == sp].biorep) for sp in SPECIES}


def og_mean(bioreps):
    M = pd.concat([OGT[b] for b in bioreps], axis=1)
    return np.log2(M / pd.Series(sf)[bioreps].values * C + 1).mean(axis=1)


SPM = {sp: og_mean(QC1[sp]) for sp in SPECIES}
others = [s for s in SPECIES if s != "Gmax"]
m = (CN["Gmax"] == 2) & (CN[others].isin([0, 1]).all(axis=1)) & ((CN[others] == 1).sum(axis=1) >= 1)
pkG = pd.read_csv(f"{RES}/step3/Gmax/pairs_ks.tsv", sep="\t")
dG = {frozenset((x[5:], y[5:])): (c, k) for x, y, c, k in zip(pkG.gene1, pkG.gene2, pkG.collinear, pkG.ks)}
genesG = {o: [g[5:] for g in v.split(", ")] for o, v in zip(og.Orthogroup, og["Gmax"]) if v}
wgd = [o for o in CN.index[m] if (lambda t: t is not None and t[0] == 1 and np.isfinite(t[1]) and t[1] < 0.3)(dG.get(frozenset(genesG[o])))]
cal = list(SC)
REF = {o: np.mean([SPM[s][o] for s in others if CN.at[o, s] == 1]) for o in wgd + cal}   # reference level (log2 scale)

# ------------------------------------------------------------------ R1 dosage by expression bin
say("== R1. Relative dosage by reference expression level (EXPLORATORY, post hoc)")
reg = pd.read_csv(f"{RES}/results_step5/tests_step5.tsv", sep="\t").set_index("analysis").estimate
rows = []
for lab, bio, key in (("study 1", QC1["Gmax"], "3_dosage_Gmax"), ("study 2", list(qc2[qc2.species == "Gmax"].biorep), "3_dosage_Gmax_study2")):
    Gm = og_mean(bio)
    D = pd.Series({o: Gm[o] - REF[o] for o in wgd + cal})
    stat = np.median(D[wgd]) - np.median(D[cal])
    if abs(stat - reg[key]) > 1e-9:
        sys.exit(f"registered dosage not reproduced for {lab}: {stat} vs {reg[key]}")
    say(f"  {lab}: registered estimate reproduced ({stat:+.3f}); WGD {len(wgd)}, calibration {len(cal)}")
    ref = pd.Series(REF)
    edges = np.quantile(ref[wgd], [0, .25, .5, .75, 1.0]); edges[0], edges[-1] = -np.inf, np.inf   # quartiles of the WGD set
    for i in range(4):
        inb = (ref > edges[i]) & (ref <= edges[i + 1])
        w = [o for o in wgd if inb[o]]; c = [o for o in cal if inb[o]]
        if len(w) < 20 or len(c) < 20:
            continue
        Dw, Dc = D[w].values, D[c].values
        est = np.median(Dw) - np.median(Dc)
        bt = [np.median(rng.choice(Dw, len(Dw))) - np.median(rng.choice(Dc, len(Dc))) for _ in range(NBOOT)]
        x = 2 ** ref[w].values - 1                                    # reference per-copy level, scaled TPM
        full = np.median(np.log2(2 * x + 1) - np.log2(x + 1))         # full dosage under the same pseudocount
        rows.append(dict(study=lab, bin=i + 1, ref_tpm_lo=round(2 ** max(edges[i], 0) - 1, 2) if np.isfinite(edges[i]) else 0.0,
                         ref_tpm_hi=round(2 ** edges[i + 1] - 1, 2) if np.isfinite(edges[i + 1]) else np.inf,
                         n_wgd=len(w), n_cal=len(c), rel_dosage=est, lo=np.quantile(bt, .025), hi=np.quantile(bt, .975),
                         full_dosage_expected=full, fraction_of_full=est / full))
R1 = pd.DataFrame(rows); R1.to_csv(f"{OUT}/R1_dosage_by_expression.tsv", sep="\t", index=False)
say(R1.round(3).to_string(index=False))

# ------------------------------------------------------------------ R2 absolute-scale gradient
say("\n== R2. Absolute-scale copy-number gradient (EXPLORATORY, post hoc)")
def grad(f):
    d = pd.read_csv(f, sep="\t")
    s, mlt = d.loc[d.copies == 1, "div"], d.loc[d.copies >= 2, "div"]
    return np.median(mlt) - np.median(s), np.log2(np.median(mlt) / np.median(s)), len(s), len(mlt)
W = {}
for sp in SPECIES:
    f = f"{RES}/results_step2/div_within_primary_{sp}.tsv.gz"
    if os.path.exists(f):
        W[sp] = grad(f)
TT = pd.read_csv("metadata/step2_divergence_times.tsv", sep="\t")
rows = []
for r in TT.itertuples():
    a, b = sorted([r.sp_a, r.sp_b], key=SPECIES.index)
    dif, lr, ns, nm = grad(f"{RES}/results_step2/div_between_primary_{a}_{b}.tsv.gz")
    wd = [W[s][0] for s in (a, b) if s in W]; wl = [W[s][1] for s in (a, b) if s in W]
    rows.append(dict(sp_a=a, sp_b=b, time=r.median_mya, n_single=ns, n_multi=nm, diff_between=dif, logratio_between=lr,
                     diff_within_mean=np.mean(wd), logratio_within_mean=np.mean(wl),
                     excess_diff=dif - np.mean(wd), excess_logratio=lr - np.mean(wl)))
R2 = pd.DataFrame(rows); R2.to_csv(f"{OUT}/R2_absolute_gradient_pairs.tsv", sep="\t", index=False)
pd.DataFrame([dict(species=k, diff_within=v[0], logratio_within=v[1], n_single=v[2], n_multi=v[3]) for k, v in W.items()]).to_csv(
    f"{OUT}/R2_absolute_gradient_within.tsv", sep="\t", index=False)
say("  within studies: diff " + ", ".join(f"{k} {v[0]:+.3f}" for k, v in W.items()))
say(f"  within, median diff {np.median([v[0] for v in W.values()]):+.3f}, median log2 ratio {np.median([v[1] for v in W.values()]):+.3f}")
say(f"  between, median diff {R2.diff_between.median():+.3f}, median log2 ratio {R2.logratio_between.median():+.3f}")
say(f"  between minus within: median excess diff {R2.excess_diff.median():+.3f} (positive in {(R2.excess_diff > 0).sum()}/{len(R2)}); "
    f"median excess log2 ratio {R2.excess_logratio.median():+.3f} (positive in {(R2.excess_logratio > 0).sum()}/{len(R2)})")

P = R2[(R2.sp_a != "Ccaj") & (R2.sp_b != "Ccaj")].copy()
P["key"] = P.sp_a + " " + P.sp_b
rngT = P.time.max() - P.time.min()
sp_list = sorted(set(P.sp_a) | set(P.sp_b), key=SPECIES.index)
def change(d, col):
    return np.polyfit(d.time, d[col], 1)[0] * rngT if len(d) >= 3 and d.time.var() > 0 else np.nan
def boot(col):
    idx = P.set_index("key")
    out = []
    for _ in range(NBOOT):
        s = list(rng.choice(sp_list, len(sp_list), replace=True))
        keys = [" ".join(sorted([s[i], s[j]], key=SPECIES.index)) for i, j in itertools.combinations(range(len(s)), 2) if s[i] != s[j]]
        d = idx.loc[[k for k in keys if k in idx.index]]
        out.append(change(d, col))
    return np.nanquantile(out, [.025, .975])
trend = []
for col, lab in (("diff_between", "between-species gradient, difference"), ("excess_diff", "excess over within, difference"),
                 ("logratio_between", "between-species gradient, log2 ratio"), ("excess_logratio", "excess over within, log2 ratio")):
    est = change(P, col); lo, hi = boot(col)
    trend.append(dict(measure=col, n_pairs=len(P), change=est, lo=lo, hi=hi))
    say(f"  change over {P.time.min():.1f}-{P.time.max():.1f} Myr ({len(P)} pairs), {lab}: {est:+.3f} [{lo:+.3f}, {hi:+.3f}]")
pd.DataFrame(trend).to_csv(f"{OUT}/R2_time_trend.tsv", sep="	", index=False)

# ------------------------------------------------------------------ R3 absolute gradient: two studies of one species vs species, identical orthogroups
say("\n== R3. Absolute-scale gradient between studies vs between species on identical orthogroups (EXPLORATORY, post hoc)")
def diff_ratio(div, copies):
    s, mlt = div[copies == 1], div[copies >= 2]
    return np.median(mlt) - np.median(s), np.log2(np.median(mlt) / np.median(s))
rows = []
for A in ("Gmax", "Pvul", "Mtru"):
    st = pd.read_csv(f"{RES}/results_step4/div_study_{A}.tsv.gz", sep="\t").set_index("Orthogroup")
    for B in SPECIES:
        if B == A:
            continue
        a, b = sorted([A, B], key=SPECIES.index)
        d = pd.read_csv(f"{RES}/results_step2/div_between_primary_{a}_{b}.tsv.gz", sep="\t").set_index("Orthogroup")
        d = d[d.index.isin(st.index)]
        ds, rs = diff_ratio(d["div"].values, d.copies.values)
        dt, rt = diff_ratio(st.loc[d.index, "div_study"].values, d.copies.values)
        rows.append(dict(focal=A, other=B, n_og=len(d), diff_species=ds, diff_study=dt, logratio_species=rs, logratio_study=rt))
R3 = pd.DataFrame(rows); R3.to_csv(f"{OUT}/R3_absolute_study_vs_species.tsv", sep="\t", index=False)
for A, g in R3.groupby("focal", sort=False):
    say(f"  {A}: median diff study {g.diff_study.median():+.3f} vs species {g.diff_species.median():+.3f} (species > study in {(g.diff_species > g.diff_study).sum()}/{len(g)}); "
        f"log2 ratio study {g.logratio_study.median():+.3f} vs species {g.logratio_species.median():+.3f}")
say(f"  all 33: diff species minus study, median {(R3.diff_species - R3.diff_study).median():+.3f} (positive in {(R3.diff_species > R3.diff_study).sum()}/33); "
    f"log2 ratio species minus study, median {(R3.logratio_species - R3.logratio_study).median():+.3f} (positive in {(R3.logratio_species > R3.logratio_study).sum()}/33)")
