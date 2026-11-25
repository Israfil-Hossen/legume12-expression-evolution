#!/usr/bin/env python3
"""EXPLORATORY (added 2026-11-24, after all tests specified in advance were known; prompted by a pre-submission review).

Is the conservation of ohnolog expression asymmetry (2B), and its decline with divergence time, specific to duplicates?
Null: pairs of NON-paralogous genes, i.e. two different orthogroups that are single-copy in all 12 legumes. Each of the
572 ohnolog orthogroups is matched to a random gene pair with (i) a similar mean expression of the two genes, (ii) a
similar mean absolute log2 ratio between them, and (iii) the same set of species in which the pair is available.
For every matched null set the 2B statistics are recomputed exactly as for the ohnologs: the per-species mean log2
ratio, the Spearman correlation of mean ratios between every two species (rho_r) and the change in rho_r over the
observed divergence-time range (OLS slope x range; 55 species pairs without Cajanus cajan).

Expression is computed as in pipeline/12_step5/step5_analysis.py (gene TPM = sum of transcript TPM; samples scaled so
that their median over the orthogroups single-copy in all 12 species equals the grand median over Step 2 samples;
between-species study only, samples that passed quality control). The script first reproduces the registered 2B
estimate (-0.224) from the same data, as a check.

Usage (repository root): python pipeline/12_step5/explore_2B_random_pair_null.py <DATA_dir with quant/, tables/, proteomes/>
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

DATA = sys.argv[1]
RES, OUT = "results", "results/exploratory_review"
NSETS = 1000
rng = np.random.default_rng(20261124)
SPECIES = ["Gmax", "Pvul", "Vrad", "Vung", "Ccaj", "Mtru", "Cari", "Psat", "Tpra", "Ljap", "Lang", "Adur"]
OFR = os.path.join(DATA, "proteomes", "OrthoFinder", "Results_Sep23")
LOG = open(os.path.join(OUT, "report_R7_random_pair_null.txt"), "w", encoding="utf-8")


def say(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


# ------------------------------------------------------------------ samples, orthogroups, expression (as step5_analysis.py)
qc1 = pd.read_csv(f"{RES}/results_step2/qc_samples.tsv", sep="\t")
qc1 = qc1[qc1.qc_pass]
QC1 = {sp: list(qc1[qc1.species == sp].biorep) for sp in SPECIES}
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
say(f"orthogroups single-copy in all 12 legumes: {len(SC)}")
tx2g = {sp: pd.read_csv(os.path.join(DATA, "tables", f"{sp}_tx2gene.tsv"), sep="\t", header=None, names=["Name", "gene"])
        for sp in SPECIES}
GENE, SF = {}, {}
for sp in SPECIES:
    for b in QC1[sp]:
        q = pd.read_csv(os.path.join(DATA, "quant", b, "quant.sf"), sep="\t", usecols=["Name", "TPM"])
        gt = q.merge(tx2g[sp], on="Name").groupby("gene").TPM.sum()
        gt = gt[gt.index.map(lambda g: g in gene2og[sp])]
        GENE[b] = gt
        ot = gt.groupby(gt.index.map(gene2og[sp])).sum()
        SF[b] = np.median(ot.reindex(SC, fill_value=0.0).values)
C = np.median([SF[b] for sp in SPECIES for b in QC1[sp]])
say(f"scaling constant C = {C:.4f}  (Step 4/5 report: 8.3075)")

# per species: mean scaled TPM and mean log2(scaled TPM + 1) of every gene
MT, ML = {}, {}
for sp in SPECIES:
    M = pd.concat([GENE[b] / SF[b] * C for b in QC1[sp]], axis=1).fillna(0.0)
    MT[sp] = M.mean(axis=1)
    ML[sp] = np.log2(M + 1).mean(axis=1)

TT = pd.read_csv("metadata/step2_divergence_times.tsv", sep="\t")
TIME = {(r.sp_a, r.sp_b): r.median_mya for r in TT.itertuples()}
TIME.update({(b, a): t for (a, b), t in list(TIME.items())})
SPP = [(a, b) for i, a in enumerate(SPECIES) for b in SPECIES[i + 1:]]
PRIM = [k for k, (a, b) in enumerate(SPP) if "Ccaj" not in (a, b)]
tvec = np.array([TIME[p] for p in SPP])[PRIM]
trange = tvec.max() - tvec.min()


FITS = []


def stats_2b(R, keep_fit=False):
    """R: DataFrame pairs x species of mean log2 ratios (NaN = not available). Returns rho_r per species pair and trend."""
    rr = []
    for a, b in SPP:
        ok = R[a].notna() & R[b].notna()
        rr.append(stats.spearmanr(R.loc[ok, a], R.loc[ok, b]).correlation if ok.sum() > 2 else np.nan)
    rr = np.array(rr)
    y = rr[PRIM]
    slope, icpt = np.polyfit(tvec, y, 1)
    if keep_fit:
        FITS.append((slope, icpt))
    return rr, slope * trange


# ------------------------------------------------------------------ ohnolog pairs exactly as step5_analysis.py
h = pd.read_csv(os.path.join(OFR, "Phylogenetic_Hierarchical_Orthogroups", "N1.tsv"), sep="\t", dtype=str).fillna("")
ncopy = h[SPECIES].apply(lambda c: c.map(lambda v: 0 if v == "" else len(v.split(", "))))
h["n_sc"] = (ncopy == 1).sum(axis=1)
cnt = h.groupby("OG").size()
h2 = h[h.OG.isin(cnt[cnt == 2].index)]
keep = h2.groupby("OG").n_sc.min()
keep = keep[keep >= 10].index
h2 = h2[h2.OG.isin(keep)].sort_values(["OG", "HOG"])
rows = []
for o, sub in h2.groupby("OG"):
    hA, hB = sub.iloc[0], sub.iloc[1]
    for sp in SPECIES:
        a, b = hA[sp], hB[sp]
        if a and b and "," not in a and "," not in b:
            gA, gB = a[len(sp) + 1:], b[len(sp) + 1:]
            tA, tB = MT[sp].get(gA, 0.0), MT[sp].get(gB, 0.0)
            if tA >= 1 or tB >= 1:          # expression filter of Step 5
                rows.append(dict(OG=o, species=sp, r=ML[sp].get(gA, 0.0) - ML[sp].get(gB, 0.0),
                                 lev=(ML[sp].get(gA, 0.0) + ML[sp].get(gB, 0.0)) / 2))
OH = pd.DataFrame(rows)
R_oh = OH.pivot(index="OG", columns="species", values="r").reindex(columns=SPECIES)
rr_oh, trend_oh = stats_2b(R_oh)
say(f"ohnolog orthogroups: {R_oh.shape[0]}; species-level pairs: {int(R_oh.notna().sum().sum())}")
say(f"REPRODUCTION CHECK: median rho_r {np.nanmedian(rr_oh):+.3f} (registered +0.562); change over {trange:.1f} Myr "
    f"{trend_oh:+.4f} (registered -0.2235)")
ref = pd.read_csv(f"{RES}/results_step5/2B_species_pairs.tsv", sep="\t")
dev = np.nanmax(np.abs(rr_oh - ref.rho_r.values))
say(f"  max |rho_r - committed rho_r| over the 66 species pairs: {dev:.2e}")

# features used for matching: mean expression level and mean |ratio| over the species where the pair is available
feat_oh = pd.DataFrame(dict(lev=OH.groupby("OG").lev.mean(), absr=OH.groupby("OG").r.apply(lambda v: np.mean(np.abs(v)))))
avail_oh = R_oh.notna()

# ------------------------------------------------------------------ candidate null pairs: two different single-copy orthogroups
SCG = {sp: pd.Series({o: g for g, o in gene2og[sp].items() if o in set(SC)}) for sp in SPECIES}
Lsc = pd.DataFrame({sp: ML[sp].reindex(SCG[sp].reindex(SC).values).values for sp in SPECIES}, index=SC)
Tsc = pd.DataFrame({sp: MT[sp].reindex(SCG[sp].reindex(SC).values).values for sp in SPECIES}, index=SC).fillna(0.0)
Lsc = Lsc.fillna(0.0)
NPOOL = 200000
i1 = rng.integers(0, len(SC), NPOOL)
i2 = rng.integers(0, len(SC), NPOOL)
ok = i1 != i2
i1, i2 = i1[ok], i2[ok]
L1, L2 = Lsc.values[i1], Lsc.values[i2]
T1, T2 = Tsc.values[i1], Tsc.values[i2]
expressed = (T1 >= 1) | (T2 >= 1)                        # same expression filter, per species
Rpool = np.where(expressed, L1 - L2, np.nan)
LEVpool = np.where(expressed, (L1 + L2) / 2, np.nan)
pool_lev = np.nanmean(LEVpool, axis=1)
pool_abs = np.nanmean(np.abs(Rpool), axis=1)
say(f"null pool: {len(i1)} random pairs of single-copy orthogroups")

# 2-D matching: 10 x 10 quantile bins of (level, |ratio|) defined on the ohnolog pairs
qlev = np.quantile(feat_oh.lev, np.linspace(0, 1, 11)); qabs = np.quantile(feat_oh.absr, np.linspace(0, 1, 11))
qlev[0], qabs[0], qlev[-1], qabs[-1] = -np.inf, -np.inf, np.inf, np.inf
bin_pool = np.digitize(pool_lev, qlev[1:-1]) * 10 + np.digitize(pool_abs, qabs[1:-1])
bin_oh = np.digitize(feat_oh.lev.values, qlev[1:-1]) * 10 + np.digitize(feat_oh.absr.values, qabs[1:-1])
pool_by_bin = {b: np.where(bin_pool == b)[0] for b in np.unique(bin_oh)}
small = {b: len(v) for b, v in pool_by_bin.items() if len(v) < 50}
say(f"matching bins used: {len(pool_by_bin)}; bins with < 50 candidate null pairs: {small if small else 'none'}")

trends, medians = [], []
spp_mat = np.array([[SPECIES.index(a), SPECIES.index(b)] for a, b in SPP])
for s in range(NSETS):
    pick = np.array([rng.choice(pool_by_bin[b]) for b in bin_oh])
    Rn = pd.DataFrame(Rpool[pick], columns=SPECIES, index=feat_oh.index)
    Rn = Rn.where(avail_oh.reindex(index=feat_oh.index).values)     # same species availability as the ohnolog it replaces
    rr, tr = stats_2b(Rn, keep_fit=True)
    trends.append(tr)
    medians.append(np.nanmedian(rr))
trends, medians = np.array(trends), np.array(medians)

# unmatched null (random pairs, availability only), for reference
trends_u, medians_u = [], []
for s in range(200):
    pick = rng.integers(0, len(i1), len(feat_oh))
    Rn = pd.DataFrame(Rpool[pick], columns=SPECIES, index=feat_oh.index).where(avail_oh.reindex(index=feat_oh.index).values)
    rr, tr = stats_2b(Rn)
    trends_u.append(tr)
    medians_u.append(np.nanmedian(rr))

p_trend = (1 + np.sum(trends <= trend_oh)) / (1 + NSETS)
p_med = (1 + np.sum(medians >= np.nanmedian(rr_oh))) / (1 + NSETS)
say("")
say("RESULT (matched null: non-paralogous gene pairs matched on expression level, |ratio| and species availability)")
say(f"  ohnologs:   median rho_r {np.nanmedian(rr_oh):+.3f}; change over {trange:.1f} Myr {trend_oh:+.3f}")
say(f"  null sets ({NSETS}): median rho_r {np.median(medians):+.3f} (2.5-97.5%: {np.quantile(medians, .025):+.3f} to "
    f"{np.quantile(medians, .975):+.3f}); change {np.median(trends):+.3f} ({np.quantile(trends, .025):+.3f} to "
    f"{np.quantile(trends, .975):+.3f})")
say(f"  one-sided empirical P: ohnolog conservation higher than null {p_med:.4f}; ohnolog decline steeper than null {p_trend:.4f}")
say(f"  unmatched random pairs (200 sets): median rho_r {np.median(medians_u):+.3f}; change {np.median(trends_u):+.3f} "
    f"({np.quantile(trends_u, .025):+.3f} to {np.quantile(trends_u, .975):+.3f})")
pd.DataFrame(dict(set=np.arange(NSETS), median_rho_r=medians, change=trends)).to_csv(
    f"{OUT}/R7_random_pair_null.tsv", sep="\t", index=False)
pd.DataFrame([dict(statistic="median_rho_r", ohnolog=np.nanmedian(rr_oh), null_median=np.median(medians),
                   null_lo=np.quantile(medians, .025), null_hi=np.quantile(medians, .975), p_one_sided=p_med),
              dict(statistic="change_over_range", ohnolog=trend_oh, null_median=np.median(trends),
                   null_lo=np.quantile(trends, .025), null_hi=np.quantile(trends, .975), p_one_sided=p_trend)]).to_csv(
    f"{OUT}/R7_random_pair_null_summary.tsv", sep="\t", index=False)
# fitted null lines over the observed time range (for Fig. 5c): median and 2.5-97.5% of the fitted values
grid = np.linspace(tvec.min(), tvec.max(), 50)
fv = np.array([s_ * grid + i_ for s_, i_ in FITS])
pd.DataFrame(dict(time=grid, median=np.median(fv, axis=0), lo=np.quantile(fv, .025, axis=0), hi=np.quantile(fv, .975, axis=0))).to_csv(
    f"{OUT}/R7_random_pair_null_lines.tsv", sep="	", index=False)
say("written: results/exploratory_review/R7_random_pair_null.tsv, R7_random_pair_null_summary.tsv, R7_random_pair_null_lines.tsv")
