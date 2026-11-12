#!/usr/bin/env python3
"""Step 5 analyses (docs/step5_spec.md): 1A laboratory share of expression distance, 1B expression-level
control, 1C read ambiguity, 2A-2C ohnolog within-sample ratios, 3 relative dosage.

Inputs: committed results (results/results_step2, results_step4), the deposit archive contents
(quant.sf, ambig_info.tsv, tables, OrthoFinder tables) unpacked to DATA, metadata/, the IQ-TREE tree.
Details the specification leaves open are marked IMPLEMENTATION.
Usage: python pipeline/12_step5/step5_analysis.py <DATA_dir>   (from the repository root)
"""
import os, re, sys, math
import numpy as np
import pandas as pd
from scipy import stats

DATA = sys.argv[1]
RES, OUT = "results", "results/results_step5"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261009)
NDRAW, NBOOT, MATERIAL, NEGLIGIBLE = 50, 2000, 0.050, 0.025
SPECIES = ["Gmax", "Pvul", "Vrad", "Vung", "Ccaj", "Mtru", "Cari", "Psat", "Tpra", "Ljap", "Lang", "Adur"]
FOCAL = ["Gmax", "Pvul", "Mtru"]
OFR = os.path.join(DATA, "proteomes", "OrthoFinder", "Results_Sep23")
LOG = open(os.path.join(OUT, "report_step5.txt"), "w", encoding="utf-8")


def say(s=""):
    print(s); LOG.write(s + "\n"); LOG.flush()


def srho(x, y):
    return stats.spearmanr(x, y).correlation


def partial_spearman(x, y, covs):
    """Spearman partial correlation: Pearson correlation of rank residuals after regressing on covariate ranks."""
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    Z = np.column_stack([np.ones(len(rx))] + [stats.rankdata(c) for c in covs])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None); by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    return np.corrcoef(rx - Z @ bx, ry - Z @ by)[0, 1]


def boot_median(v):
    v = np.asarray([x for x in v if np.isfinite(x)])
    bt = [np.median(rng.choice(v, len(v), replace=True)) for _ in range(NBOOT)]
    return np.median(v), np.quantile(bt, .025), np.quantile(bt, .975)


say(f"Step 5 analysis   data: {DATA}")

# ============================================================ samples
qc1 = pd.read_csv(f"{RES}/results_step2/qc_samples.tsv", sep="\t"); qc1 = qc1[qc1.qc_pass].copy(); qc1["study"] = 1
qc2 = pd.read_csv(f"{RES}/results_step4/qc_samples_step4.tsv", sep="\t"); qc2 = qc2[qc2.qc_pass].copy(); qc2["study"] = 2
S = pd.concat([qc1[["species", "group", "biorep", "study"]], qc2[["species", "group", "biorep", "study"]]], ignore_index=True)
say(f"samples: study 1 {len(qc1)}, study 2 {len(qc2)}")

# ============================================================ orthogroups (Step 2 run)
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
say(f"orthogroups with a legume copy: {len(CN)}; single-copy in all 12: {len(SC)}")

# ============================================================ expression (gene TPM, orthogroup TPM, scaling as Step 2/4)
tx2g = {sp: pd.read_csv(os.path.join(DATA, "tables", f"{sp}_tx2gene.tsv"), sep="\t", header=None, names=["Name", "gene"])
        for sp in SPECIES}
GENE, OGT, AMB = {}, {}, {}
for _, s in S.iterrows():
    b, sp = s.biorep, s.species
    q = pd.read_csv(os.path.join(DATA, "quant", b, "quant.sf"), sep="\t", usecols=["Name", "TPM"])
    a = pd.read_csv(os.path.join(DATA, "quant", b, "aux_info", "ambig_info.tsv"), sep="\t")
    # IMPLEMENTATION: ambig_info rows follow quant.sf order; the extra trailing rows are genome decoys (zero counts, checked)
    a = a.iloc[:len(q)].reset_index(drop=True)
    q["U"], q["A"] = a.UniqueCount.values, a.AmbigCount.values
    m = q.merge(tx2g[sp], on="Name")
    gt = m.groupby("gene")[["TPM", "U", "A"]].sum()
    gt["og"] = gt.index.map(gene2og[sp])
    gt = gt[gt.og.notna()]
    GENE[b] = gt.TPM
    ot = gt.groupby("og")[["TPM", "U", "A"]].sum()
    OGT[b] = ot.TPM.reindex(CN.index[CN[sp] > 0], fill_value=0.0)
    AMB[b] = (ot.A / (ot.U + ot.A)).replace([np.inf, -np.inf], np.nan)
sf = {b: np.median(OGT[b].reindex(SC).values) for b in S.biorep}
C = np.median([sf[b] for b in qc1.biorep])
say(f"scaling constant C = {C:.4f}  (Step 4 report: 8.3075)")


def og_matrix(sp, bioreps):
    M = pd.concat([OGT[b] for b in bioreps], axis=1); M.columns = bioreps
    return np.log2(M / pd.Series(sf)[bioreps] * C + 1)


def gene_matrix(genes, bioreps):
    M = pd.concat([GENE[b].reindex(genes, fill_value=0.0) for b in bioreps], axis=1); M.columns = bioreps
    return M / pd.Series(sf)[bioreps] * C        # scaled TPM (not logged)


# ============================================================ draws (identical rules to Step 2)
def split_draw(groups, k):
    gs = {}
    for b, g in groups.items():
        gs.setdefault(g, []).append(b)
    names = sorted(gs); cap = {g: len(gs[g]) // 2 for g in names}
    if sum(cap.values()) < k:
        return None
    al = {g: 0 for g in names}
    for _ in range(k):
        av = [g for g in names if cap[g] - al[g] > 0]; pk = av[rng.integers(len(av))]; al[pk] += 1
    h1, h2 = [], []
    for g in names:
        if al[g]:
            s = list(rng.choice(gs[g], 2 * al[g], replace=False)); h1 += s[:al[g]]; h2 += s[al[g]:]
    return h1, h2


def draws_within(groups, k=2):
    out = []
    for _ in range(NDRAW):
        h = split_draw(groups, k)
        if h is None:
            return None
        out.append(h)
    return out


def draws_between(cols_a, cols_b, k=2):
    return [(list(rng.choice(cols_a, k, replace=False)), list(rng.choice(cols_b, k, replace=False))) for _ in range(NDRAW)]


def div_from_draws(X, draws):
    """X: DataFrame rows x samples; mean over draws of |mean half1 - mean half2| / sqrt(2)."""
    acc = 0
    for h1, h2 in draws:
        acc = acc + (X[h1].mean(axis=1) - X[h2].mean(axis=1)).abs() / math.sqrt(2)
    return acc / len(draws)


def pair_file(a, b):
    s = sorted([a, b], key=SPECIES.index)
    return f"{RES}/results_step2/div_between_primary_{s[0]}_{s[1]}.tsv.gz"


def within_file(sp):
    return f"{RES}/results_step2/div_within_primary_{sp}.tsv.gz"


TT = pd.read_csv("metadata/step2_divergence_times.tsv", sep="\t")
TIME = {(r.sp_a, r.sp_b): r.median_mya for r in TT.itertuples()}
TIME.update({(b, a): t for (a, b), t in list(TIME.items())})
B2 = pd.read_csv(f"{RES}/results_step2/rho_between_primary.tsv", sep="\t")
W2 = pd.read_csv(f"{RES}/results_step2/rho_within_primary.tsv", sep="\t")


def species_boot_trend(df, xcol, ycol):
    """Step 2 primary-test logic: OLS slope x observed range; species bootstrap with pairs among resampled species."""
    rng_x = df[xcol].max() - df[xcol].min()
    def fit(d):
        if len(d) < 3 or d[xcol].var() == 0:
            return np.nan
        return np.polyfit(d[xcol], d[ycol], 1)[0] * rng_x
    est = fit(df); sp = sorted(set(df.sp_a) | set(df.sp_b))
    key = {frozenset((a, b)): i for i, (a, b) in enumerate(zip(df.sp_a, df.sp_b))}
    bt = []
    for _ in range(NBOOT):
        s = list(rng.choice(sp, len(sp), replace=True)); idx = []
        for i in range(len(s)):
            for j in range(i + 1, len(s)):
                if s[i] != s[j] and frozenset((s[i], s[j])) in key:
                    idx.append(key[frozenset((s[i], s[j]))])
        bt.append(fit(df.iloc[idx]) if idx else np.nan)
    return est, np.nanquantile(bt, .025), np.nanquantile(bt, .975), rng_x


TESTS = []

# ============================================================ 1A laboratory share of expression distance
say("\n== 1A. Expression distance within study, between studies, between species (Step 2 pair sets)")
study = {A: pd.read_csv(f"{RES}/results_step4/div_study_{A}.tsv.gz", sep="\t").set_index("Orthogroup") for A in FOCAL}
win = {sp: pd.read_csv(within_file(sp), sep="\t").set_index("Orthogroup") for sp in SPECIES if os.path.exists(within_file(sp))}
rows, boots = [], {}
for A in FOCAL:
    for Bsp in [s for s in SPECIES if s != A]:
        d = pd.read_csv(pair_file(A, Bsp), sep="\t").set_index("Orthogroup")
        ds, dst, dw = d["div"].values, study[A].loc[d.index, "div_study"].values, win[A].loc[d.index, "div"].values
        f = lambda i: (np.median(dst[i]) - np.median(dw[i])) / (np.median(ds[i]) - np.median(dw[i]))
        n = len(ds); allidx = np.arange(n)
        boots[(A, Bsp)] = np.array([f(rng.integers(0, n, n)) for _ in range(NBOOT)])
        rows.append(dict(focal=A, other=Bsp, n_og=n, d_within=np.median(dw), d_study=np.median(dst), d_species=np.median(ds), F=f(allidx)))
P1A = pd.DataFrame(rows); P1A.to_csv(f"{OUT}/1A_pairs.tsv", sep="\t", index=False)
per = {A: np.mean([boots[(A, b)] for b in P1A[P1A.focal == A].other], axis=0) for A in FOCAL}
Fall = np.mean([per[A] for A in FOCAL], axis=0)
Fest = P1A.groupby("focal").F.mean().reindex(FOCAL).mean()
lo, hi = np.quantile(Fall, [.025, .975])
v = ("laboratory accounts for at least half of the between-species excess" if Fest >= 0.5 and lo > 0.25 else
     "evolution dominates the between-species excess" if Fest <= 0.25 else "partial")
for A in FOCAL:
    r = P1A[P1A.focal == A]
    say(f"  {A}: median distance within {r.d_within.iloc[0]:.3f} | between studies {r.d_study.median():.3f} | "
        f"between species median {r.d_species.median():.3f};  F mean {r.F.mean():+.3f} "
        f"[{np.quantile(per[A], .025):+.3f}, {np.quantile(per[A], .975):+.3f}]")
say(f"  F (all 33 pairs) {Fest:+.3f} [{lo:+.3f}, {hi:+.3f}]  -> {v}   (PRIMARY 1A)")
TESTS.append(dict(analysis="1A_F", estimate=Fest, lo=lo, hi=hi, verdict=v))
# secondary: noise-corrected between-species distance against divergence time
rows = []
for r in B2.itertuples():
    a, b = r.sp_a, r.sp_b
    if "Ccaj" in (a, b) or a not in win or b not in win:
        continue
    d = pd.read_csv(pair_file(a, b), sep="\t").set_index("Orthogroup")
    corr = np.median(d["div"]) - (np.median(win[a].loc[d.index, "div"]) + np.median(win[b].loc[d.index, "div"])) / 2
    rows.append(dict(sp_a=a, sp_b=b, time=TIME[(a, b)], d_species=np.median(d["div"]), d_corrected=corr))
T1A = pd.DataFrame(rows); T1A.to_csv(f"{OUT}/1A_time.tsv", sep="\t", index=False)
est, lo, hi, rx = species_boot_trend(T1A, "time", "d_corrected")
say(f"  secondary: noise-corrected distance ~ time ({len(T1A)} pairs, range {rx:.1f} Myr): change {est:+.4f} log2 [{lo:+.4f}, {hi:+.4f}]  (reported)")
TESTS.append(dict(analysis="1A_time_secondary", estimate=est, lo=lo, hi=hi, verdict="reported"))

# ============================================================ 1B expression-level control
say("\n== 1B. Copy-number gradient with expression level held constant (partial Spearman)")
QC1 = qc1.groupby("species").biorep.apply(list).to_dict()
MEANEXP, PERCOPY = {}, {}
for sp in SPECIES:
    M = og_matrix(sp, QC1[sp])
    MEANEXP[sp] = M.mean(axis=1)
    copies = CN.loc[M.index, sp]
    PERCOPY[sp] = np.log2((2 ** M - 1).div(copies, axis=0) + 1).mean(axis=1)
rows = []
for sp, w in win.items():
    e, pc = MEANEXP[sp].loc[w.index].values, PERCOPY[sp].loc[w.index].values
    raw = srho(w.copies, w["div"])
    reg = float(W2.loc[W2.species == sp, "rho"].iloc[0])
    assert abs(raw - reg) < 1e-6, (sp, raw, reg)       # same values as the registered Step 2 gradient
    dec = pd.qcut(e, 10, labels=False, duplicates="drop")
    rows.append(dict(level="within", label=sp, rho_raw=raw, rho_partial_expr=partial_spearman(w.copies, w["div"], [e]),
                     rho_partial_percopy=partial_spearman(w.copies, w["div"], [pc]),
                     rho_decile_median=np.median([srho(w.copies.values[dec == k], w["div"].values[dec == k]) for k in np.unique(dec)])))
for r in B2.itertuples():
    d = pd.read_csv(pair_file(r.sp_a, r.sp_b), sep="\t").set_index("Orthogroup")
    e = (MEANEXP[r.sp_a].loc[d.index].values + MEANEXP[r.sp_b].loc[d.index].values) / 2
    pc = (PERCOPY[r.sp_a].loc[d.index].values + PERCOPY[r.sp_b].loc[d.index].values) / 2
    raw = srho(d.copies, d["div"]); assert abs(raw - r.rho) < 1e-6
    rows.append(dict(level="between", label=f"{r.sp_a}-{r.sp_b}", rho_raw=raw, rho_partial_expr=partial_spearman(d.copies, d["div"], [e]),
                     rho_partial_percopy=partial_spearman(d.copies, d["div"], [pc]), rho_decile_median=np.nan))
P1B = pd.DataFrame(rows); P1B.to_csv(f"{OUT}/1B_partial.tsv", sep="\t", index=False)
w = P1B[P1B.level == "within"]
est, lo, hi = boot_median(w.rho_partial_expr)
v = ("the gradient persists after control for expression level" if est >= MATERIAL and lo > 0 else
     "the gradient is explained by expression level" if est <= NEGLIGIBLE else "partial")
say(f"  within species ({len(w)}): raw median {w.rho_raw.median():+.4f}; partial | mean expression {est:+.4f} [{lo:+.4f}, {hi:+.4f}]  -> {v}   (PRIMARY 1B)")
say(f"    partial | per-copy expression: median {w.rho_partial_percopy.median():+.4f};  median within-decile rho: {w.rho_decile_median.median():+.4f}")
bt = P1B[P1B.level == "between"]
say(f"  between species ({len(bt)}): raw median {bt.rho_raw.median():+.4f}; partial | mean expression {bt.rho_partial_expr.median():+.4f}; | per-copy {bt.rho_partial_percopy.median():+.4f}")
for r in w.itertuples():
    say(f"    {r.label:5s} raw {r.rho_raw:+.4f}  partial(expr) {r.rho_partial_expr:+.4f}  partial(per-copy) {r.rho_partial_percopy:+.4f}  decile median {r.rho_decile_median:+.4f}")
TESTS.append(dict(analysis="1B_partial_expr", estimate=est, lo=lo, hi=hi, verdict=v))

# ============================================================ 1C read ambiguity
say("\n== 1C. Read ambiguity (salmon AmbigCount share), per orthogroup, median over QC-passed study-1 samples")
rows = []
for sp, w in win.items():
    amb = pd.concat([AMB[b] for b in QC1[sp]], axis=1).median(axis=1).reindex(w.index)
    ok = amb.notna().values
    ww, e, am = w[ok], MEANEXP[sp].loc[w.index].values[ok], amb.values[ok]
    rows.append(dict(species=sp, n_og=int(ok.sum()), rho_copies_ambiguity=srho(ww.copies, am),
                     rho_partial_amb_expr=partial_spearman(ww.copies, ww["div"], [am, e]),
                     rho_ambiguity_div_given_copies=partial_spearman(am, ww["div"], [ww.copies])))
P1C = pd.DataFrame(rows); P1C.to_csv(f"{OUT}/1C_ambiguity.tsv", sep="\t", index=False)
est, lo, hi = boot_median(P1C.rho_partial_amb_expr)
v = ("the gradient persists after control for expression level and read ambiguity" if est >= MATERIAL and lo > 0 else
     "the gradient is explained by expression level and read ambiguity" if est <= NEGLIGIBLE else "partial")
for r in P1C.itertuples():
    say(f"    {r.species:5s} n {r.n_og:5d}  rho(copies, ambiguity) {r.rho_copies_ambiguity:+.3f}  partial gradient {r.rho_partial_amb_expr:+.4f}  "
        f"rho(ambiguity, div | copies) {r.rho_ambiguity_div_given_copies:+.3f}")
say(f"  median partial gradient | ambiguity + expression {est:+.4f} [{lo:+.4f}, {hi:+.4f}]  -> {v}   (PRIMARY 1C)")
TESTS.append(dict(analysis="1C_partial_amb_expr", estimate=est, lo=lo, hi=hi, verdict=v))

# ============================================================ Part 2: ohnolog pairs
say("\n== Part 2. Ohnolog pairs (orthogroups split into exactly two N1 HOGs, both single-copy in >= 10 species)")
h = pd.read_csv(os.path.join(OFR, "Phylogenetic_Hierarchical_Orthogroups", "N1.tsv"), sep="\t", dtype=str).fillna("")
ncopy = h[SPECIES].apply(lambda c: c.map(lambda v: 0 if v == "" else len(v.split(", "))))
h["n_sc"] = (ncopy == 1).sum(axis=1)
cnt = h.groupby("OG").size()
h2 = h[h.OG.isin(cnt[cnt == 2].index)]
keep = h2.groupby("OG").n_sc.min(); keep = keep[keep >= 10].index
h2 = h2[h2.OG.isin(keep)].sort_values(["OG", "HOG"])
PAIRS = []
for o, sub in h2.groupby("OG"):
    hA, hB = sub.iloc[0], sub.iloc[1]            # copy A = smaller HOG identifier
    for sp in SPECIES:
        a, b = hA[sp], hB[sp]
        if a and b and "," not in a and "," not in b:
            PAIRS.append(dict(OG=o, HOG_A=hA.HOG, HOG_B=hB.HOG, species=sp, gA=a[len(sp) + 1:], gB=b[len(sp) + 1:]))
PAIRS = pd.DataFrame(PAIRS)
# collinearity and Ks from the Step 3 reciprocal-best-hit table (pairs absent from it: collinear unknown -> treated as not collinear)
col = []
for sp in SPECIES:
    pk = pd.read_csv(f"{RES}/step3/{sp}/pairs_ks.tsv", sep="\t")
    d = {frozenset((x[len(sp) + 1:], y[len(sp) + 1:])): (c, k) for x, y, c, k in zip(pk.gene1, pk.gene2, pk.collinear, pk.ks)}
    for r in PAIRS[PAIRS.species == sp].itertuples():
        c, k = d.get(frozenset((r.gA, r.gB)), (0, np.nan)); col.append((r.Index, c, k))
col = pd.DataFrame(col, columns=["i", "collinear", "ks"]).set_index("i")
PAIRS = PAIRS.join(col)
# expression filter (study 1): mean scaled TPM of at least one copy >= 1
RATIO1, keep_rows = {}, []
for sp in SPECIES:
    pp = PAIRS[PAIRS.species == sp]
    X = gene_matrix(pd.unique(pd.concat([pp.gA, pp.gB])), QC1[sp])
    ok = (X.loc[pp.gA].mean(axis=1).values >= 1) | (X.loc[pp.gB].mean(axis=1).values >= 1)
    pp = pp[ok]; keep_rows += list(pp.index)
    R = pd.DataFrame(np.log2(X.loc[pp.gA].values + 1) - np.log2(X.loc[pp.gB].values + 1), index=pp.OG.values, columns=QC1[sp])
    RATIO1[sp] = R
PAIRS = PAIRS.loc[keep_rows]
PAIRS.to_csv(f"{OUT}/2_pairs.tsv", sep="\t", index=False)
say(f"  ohnolog orthogroups: {PAIRS.OG.nunique()}; species-pairs used: {len(PAIRS)}; per species: "
    + ", ".join(f"{sp} {len(RATIO1[sp])}" for sp in SPECIES))

# ---------------- 2A laboratory protection
say("\n== 2A. Laboratory excess for within-sample ratios versus absolute expression (same genes, same draws)")
rows, bootE = [], {}
for A in FOCAL:
    pp = PAIRS[PAIRS.species == A]
    b1 = QC1[A]; b2 = list(qc2[qc2.species == A].biorep)
    X = gene_matrix(pd.unique(pd.concat([pp.gA, pp.gB])), b1 + b2)
    L = np.log2(X + 1)
    R = pd.DataFrame(L.loc[pp.gA].values - L.loc[pp.gB].values, index=pp.OG.values, columns=b1 + b2)
    LA = pd.DataFrame(L.loc[pp.gA].values, index=pp.OG.values, columns=b1 + b2)
    LB = pd.DataFrame(L.loc[pp.gB].values, index=pp.OG.values, columns=b1 + b2)
    grp = dict(zip(qc1[qc1.species == A].biorep, qc1[qc1.species == A].group))
    dw, db = draws_within(grp), draws_between(b1, b2)
    rw, rb = div_from_draws(R, dw).values, div_from_draws(R, db).values
    aw = ((div_from_draws(LA, dw) + div_from_draws(LB, dw)) / 2).values
    ab = ((div_from_draws(LA, db) + div_from_draws(LB, db)) / 2).values
    E = lambda i: (np.median(rb[i]) - np.median(rw[i])) / (np.median(ab[i]) - np.median(aw[i]))
    n = len(rw); bootE[A] = np.array([E(rng.integers(0, n, n)) for _ in range(NBOOT)])
    rows.append(dict(species=A, n_pairs=n, ratio_within=np.median(rw), ratio_study=np.median(rb), abs_within=np.median(aw),
                     abs_study=np.median(ab), X_ratio=np.median(rb) - np.median(rw), X_abs=np.median(ab) - np.median(aw), E=E(np.arange(n))))
    pd.DataFrame(dict(OG=pp.OG.values, ratio_within=rw, ratio_study=rb, abs_within=aw, abs_study=ab)).to_csv(f"{OUT}/2A_divergence_{A}.tsv.gz", sep="\t", index=False)
P2A = pd.DataFrame(rows); P2A.to_csv(f"{OUT}/2A_summary.tsv", sep="\t", index=False)
Eall = np.mean([bootE[A] for A in FOCAL], axis=0); Eest = P2A.E.mean(); lo, hi = np.quantile(Eall, [.025, .975])
v = ("within-sample ratios largely remove the laboratory effect" if Eest <= 0.5 and hi < 1 else
     "within-sample ratios do not remove the laboratory effect" if Eest >= 0.9 else "partial")
for r in P2A.itertuples():
    say(f"  {r.species}: {r.n_pairs} pairs | ratio: within {r.ratio_within:.3f}, between studies {r.ratio_study:.3f} (excess {r.X_ratio:+.3f}) | "
        f"absolute: within {r.abs_within:.3f}, between studies {r.abs_study:.3f} (excess {r.X_abs:+.3f}) | E {r.E:+.3f}")
say(f"  E (mean of 3) {Eest:+.3f} [{lo:+.3f}, {hi:+.3f}]  -> {v}   (PRIMARY 2A)")
TESTS.append(dict(analysis="2A_E", estimate=Eest, lo=lo, hi=hi, verdict=v))

# ---------------- 2B conservation of asymmetry and its erosion with time
say("\n== 2B. Conservation of ohnolog asymmetry between species, and change with divergence time")
MR = {sp: RATIO1[sp].mean(axis=1) for sp in SPECIES}
pd.DataFrame(MR).rename_axis("OG").to_csv(f"{OUT}/2B_mean_ratios.tsv", sep="\t")   # per-species mean ratio (added 2026-10-09 for the figure)
COLL = {sp: set(PAIRS[(PAIRS.species == sp) & (PAIRS.collinear == 1)].OG) for sp in SPECIES}
rows = []
for i, a in enumerate(SPECIES):
    for b in SPECIES[i + 1:]:
        common = MR[a].index.intersection(MR[b].index)
        cc = [o for o in common if o in COLL[a] and o in COLL[b]]
        rows.append(dict(sp_a=a, sp_b=b, time=TIME[(a, b)], n=len(common), rho_r=srho(MR[a][common], MR[b][common]),
                         n_col=len(cc), rho_r_col=srho(MR[a][cc], MR[b][cc]) if len(cc) > 10 else np.nan))
P2B = pd.DataFrame(rows)
# patristic distances on the IQ-TREE tree (A. thaliana pruned)
def parse_newick(s):
    s = s.strip().rstrip(";"); pos = 0
    def node():
        nonlocal pos
        n = {"name": "", "len": 0.0, "children": []}
        if s[pos] == "(":
            pos += 1
            while True:
                n["children"].append(node())
                if s[pos] == ",": pos += 1; continue
                if s[pos] == ")": pos += 1; break
        m = re.match(r"[^:,()]*", s[pos:]); n["name"] = m.group(0); pos += len(m.group(0))
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":([0-9.eE+-]+)", s[pos:]); n["len"] = float(m.group(1)); pos += len(m.group(0))
        return n
    return node()
tree = parse_newick(open(f"{RES}/server_logs/concat.fa.treefile").read())
tree["children"] = [c for c in tree["children"] if c["name"] != "Atha"]
depth, path = {}, {}
def walk(n, d, anc):
    d2 = d + n["len"]; anc2 = anc + [id(n)]
    if not n["children"]:
        depth[n["name"]] = d2; path[n["name"]] = anc2; return
    for c in n["children"]:
        walk(c, d2, anc2)
nd = {}
def collect(n, d):
    nd[id(n)] = d + n["len"]
    for c in n["children"]:
        collect(c, d + n["len"])
walk(tree, 0.0, []); collect(tree, 0.0)
def patristic(a, b):
    common = [x for x, y in zip(path[a], path[b]) if x == y][-1]
    return depth[a] + depth[b] - 2 * nd[common]
P2B["patristic"] = [patristic(a, b) for a, b in zip(P2B.sp_a, P2B.sp_b)]
P2B.to_csv(f"{OUT}/2B_species_pairs.tsv", sep="\t", index=False)
prim = P2B[(P2B.sp_a != "Ccaj") & (P2B.sp_b != "Ccaj")].reset_index(drop=True)
est, lo, hi, rx = species_boot_trend(prim, "time", "rho_r")
v = ("asymmetry erodes with divergence time (evolution detectable)" if est <= -MATERIAL and hi < 0 else
     "asymmetry does not change with divergence time" if abs(est) <= NEGLIGIBLE else "partial")
say(f"  median rho_r between species: {P2B.rho_r.median():+.3f} (range {P2B.rho_r.min():+.3f} to {P2B.rho_r.max():+.3f}); pairs of ohnologs per comparison {P2B.n.min()}-{P2B.n.max()}")
say(f"  rho_r ~ divergence time ({len(prim)} pairs, range {rx:.1f} Myr): change {est:+.4f} [{lo:+.4f}, {hi:+.4f}]  -> {v}   (PRIMARY 2B)")
TESTS.append(dict(analysis="2B_trend", estimate=est, lo=lo, hi=hi, verdict=v))
pc = prim.dropna(subset=["rho_r_col"]).reset_index(drop=True)
e2, l2, h2_ = species_boot_trend(pc, "time", "rho_r_col")[:3]
say(f"  collinear-only: median rho_r {P2B.rho_r_col.median():+.3f}; trend {e2:+.4f} [{l2:+.4f}, {h2_:+.4f}]  (reported)")
e3, l3, h3, r3 = species_boot_trend(P2B, "patristic", "rho_r")
say(f"  rho_r ~ patristic distance (66 pairs, Cajanus included): change {e3:+.4f} [{l3:+.4f}, {h3:+.4f}]  (reported)")
TESTS += [dict(analysis="2B_trend_collinear", estimate=e2, lo=l2, hi=h2_, verdict="reported"),
          dict(analysis="2B_trend_patristic", estimate=e3, lo=l3, hi=h3, verdict="reported")]
# ceilings
rows = []
for sp in SPECIES:
    grp = dict(zip(qc1[qc1.species == sp].biorep, qc1[qc1.species == sp].group))
    dw = draws_within(grp)
    if dw is None:
        continue
    R = RATIO1[sp]
    rows.append(dict(species=sp, ceiling="replicate", rho=np.mean([srho(R[h1].mean(axis=1), R[h2].mean(axis=1)) for h1, h2 in dw])))
for A in FOCAL:
    pp = PAIRS[PAIRS.species == A]; b2 = list(qc2[qc2.species == A].biorep)
    X = gene_matrix(pd.unique(pd.concat([pp.gA, pp.gB])), b2)
    R2 = pd.Series((np.log2(X.loc[pp.gA].values + 1) - np.log2(X.loc[pp.gB].values + 1)).mean(axis=1), index=pp.OG.values)
    rows.append(dict(species=A, ceiling="laboratory", rho=srho(MR[A][R2.index], R2)))
CEIL = pd.DataFrame(rows); CEIL.to_csv(f"{OUT}/2B_ceilings.tsv", sep="\t", index=False)
say(f"  replicate ceiling (half vs half, 11 species): median {CEIL[CEIL.ceiling == 'replicate'].rho.median():+.3f}; "
    f"laboratory ceiling (study 1 vs study 2): " + ", ".join(f"{r.species} {r.rho:+.3f}" for r in CEIL[CEIL.ceiling == 'laboratory'].itertuples()))
# direction consistency per ohnolog orthogroup
mr = pd.DataFrame(MR)
npos, npres = (mr > 0).sum(axis=1), mr.notna().sum(axis=1)
pv = np.array([stats.binomtest(int(k), int(n), 0.5).pvalue for k, n in zip(npos, npres)])
order = np.argsort(pv); q = np.empty_like(pv); m_ = len(pv)
q[order] = np.minimum.accumulate((pv[order] * m_ / np.arange(1, m_ + 1))[::-1])[::-1]
DIR = pd.DataFrame(dict(OG=mr.index, n_species=npres.values, n_A_higher=npos.values, p=pv, q=q)); DIR.to_csv(f"{OUT}/2B_direction.tsv", sep="\t", index=False)
say(f"  direction: {len(DIR)} ohnolog orthogroups; same copy higher in all species where present: {(DIR.n_A_higher.isin([0]) | (DIR.n_A_higher == DIR.n_species)).sum()}; "
    f"binomial FDR 5% significant: {(DIR.q < 0.05).sum()} ({100 * (DIR.q < 0.05).mean():.1f}%)")

# ---------------- 2C asymmetric versus balanced pairs (descriptive)
say("\n== 2C. Consistently asymmetric versus balanced ohnolog pairs (descriptive)")
absr = mr.abs(); sgn = np.sign(mr)
frac_same = pd.concat([(sgn > 0).sum(axis=1), (sgn < 0).sum(axis=1)], axis=1).max(axis=1) / npres
cls = np.where((frac_same >= 0.8) & (absr.mean(axis=1) >= 1), "asymmetric",
               np.where(((absr < 0.5).sum(axis=1) / npres) >= 0.8, "balanced", "other"))
CL = pd.DataFrame(dict(OG=mr.index, cls=cls, mean_abs_r=absr.mean(axis=1).values, frac_same_direction=frac_same.values))
kc = PAIRS.groupby("OG").agg(ks_median=("ks", "median"), frac_collinear=("collinear", "mean"))
CL = CL.merge(kc, left_on="OG", right_index=True, how="left")
hdr = os.path.join(RES, "atha_protein_headers.txt")
if os.path.exists(hdr):
    prod = {}
    for line in open(hdr, encoding="utf-8"):
        m = re.match(r">(\S+) (.*) \[Arabidopsis thaliana\]", line.strip())
        if m:
            prod[m.group(1)] = m.group(2)
    p2g = pd.read_csv(os.path.join(DATA, "tables", "Atha_prot2gene.tsv"), sep="\t", header=None, names=["prot", "gene"])
    gprod = p2g.assign(name=p2g.prot.map(prod)).dropna().groupby("gene").name.first().to_dict()
    ogat = dict(zip(og.Orthogroup, og["Atha"]))
    TF = re.compile(r"transcription factor|MYB|WRKY|bHLH|NAC domain|AP2|homeobox|MADS|bZIP|GRAS|TCP", re.I)
    def category(o):
        names = [gprod.get(g[len("Atha_"):], "") for g in ogat.get(o, "").split(", ") if g]
        txt = " ".join(names)
        if not txt.strip():
            return "no A. thaliana annotation"
        if re.search(r"ribosomal protein", txt, re.I):
            return "ribosomal protein"
        if TF.search(txt):
            return "transcription factor"
        return "other"
    CL["category"] = CL.OG.map(category)
else:
    CL["category"] = "headers not available"
CL.to_csv(f"{OUT}/2C_classes.tsv", sep="\t", index=False)
say("  " + CL.groupby("cls").agg(n=("OG", "size"), ks_median=("ks_median", "median"), frac_collinear=("frac_collinear", "mean")).round(3).to_string().replace("\n", "\n  "))
say("  " + pd.crosstab(CL.cls, CL.category).to_string().replace("\n", "\n  "))

# ============================================================ Part 3 relative dosage
say("\n== 3. Relative dosage after lineage-specific duplication (relative to each library's transcript pool)")
SPM = {sp: og_matrix(sp, QC1[sp]).mean(axis=1) for sp in SPECIES}
def dosage(X, ogs, Xmean):
    return np.array([Xmean[o] - np.mean([SPM[s][o] for s in SPECIES if s != X and CN.at[o, s] == 1]) for o in ogs])
one_extra = {}
for X in SPECIES:
    others = [s for s in SPECIES if s != X]
    m = (CN[X] == 2) & (CN[others].isin([0, 1]).all(axis=1)) & ((CN[others] == 1).sum(axis=1) >= 1)
    one_extra[X] = list(CN.index[m])
pkG = pd.read_csv(f"{RES}/step3/Gmax/pairs_ks.tsv", sep="\t")
dG = {frozenset((x[5:], y[5:])): (c, k, t1, t2) for x, y, c, k, t1, t2 in zip(pkG.gene1, pkG.gene2, pkG.collinear, pkG.ks, pkG.type1, pkG.type2)}
genesG = {o: [g[5:] for g in v.split(", ")] for o, v in zip(og.Orthogroup, og["Gmax"]) if v}
wgd = [o for o in one_extra["Gmax"] if (lambda t: t is not None and t[0] == 1 and np.isfinite(t[1]) and t[1] < 0.3)(dG.get(frozenset(genesG[o])))]
cal = list(SC)
res3 = {}
for lab, bio in (("study 1", QC1["Gmax"]), ("study 2", list(qc2[qc2.species == "Gmax"].biorep))):
    Gm = og_matrix("Gmax", bio).mean(axis=1)
    Dw, Dc = dosage("Gmax", wgd, Gm), dosage("Gmax", cal, Gm)
    stat = np.median(Dw) - np.median(Dc)
    bt = [np.median(rng.choice(Dw, len(Dw))) - np.median(rng.choice(Dc, len(Dc))) for _ in range(NBOOT)]
    res3[lab] = (stat, *np.quantile(bt, [.025, .975]))
    pd.DataFrame(dict(OG=wgd + cal, set=["glycine_wgd"] * len(wgd) + ["calibration"] * len(cal), D=np.r_[Dw, Dc])).to_csv(
        f"{OUT}/3_dosage_Gmax_{lab.replace(' ', '')}.tsv.gz", sep="\t", index=False)
stat, lo, hi = res3["study 1"]
v = ("relative family expression rises with copy number (dosage retained)" if stat >= 0.75 else
     "relative family expression stays constant (dosage sharing)" if stat <= 0.25 else "partial")
say(f"  Glycine-WGD orthogroups (Gmax 2 copies, collinear, Ks < 0.3; others 1): {len(wgd)}; calibration (single-copy in all 12): {len(cal)}")
say(f"  relative dosage, study 1: {stat:+.3f} log2 [{lo:+.3f}, {hi:+.3f}]  -> {v}   (PRIMARY 3)")
s2 = res3["study 2"]; say(f"  replication with Gmax study 2: {s2[0]:+.3f} [{s2[1]:+.3f}, {s2[2]:+.3f}]")
TESTS += [dict(analysis="3_dosage_Gmax", estimate=stat, lo=lo, hi=hi, verdict=v),
          dict(analysis="3_dosage_Gmax_study2", estimate=s2[0], lo=s2[1], hi=s2[2], verdict="replication")]
# descriptive: tandem/proximal lineage-specific pairs in every species
rows = []
for X in SPECIES:
    pk = pd.read_csv(f"{RES}/step3/{X}/pairs_ks.tsv", sep="\t")
    dX = {frozenset((x[len(X) + 1:], y[len(X) + 1:])): (t1, t2) for x, y, t1, t2 in zip(pk.gene1, pk.gene2, pk.type1, pk.type2)}
    gX = {o: [g[len(X) + 1:] for g in v.split(", ")] for o, v in zip(og.Orthogroup, og[X]) if v}
    tand = [o for o in one_extra[X] if (lambda t: t is not None and t[0] in (2, 3) and t[1] in (2, 3))(dX.get(frozenset(gX[o])))]
    if len(tand) < 10:
        rows.append(dict(species=X, n_tandem=len(tand), rel_dosage=np.nan)); continue
    Dt, Dc = dosage(X, tand, SPM[X]), dosage(X, cal, SPM[X])
    rows.append(dict(species=X, n_tandem=len(tand), rel_dosage=np.median(Dt) - np.median(Dc)))
P3 = pd.DataFrame(rows); P3.to_csv(f"{OUT}/3_tandem_descriptive.tsv", sep="\t", index=False)
say("  tandem/proximal lineage-specific duplicates (descriptive): " + ", ".join(f"{r.species} n={r.n_tandem} {r.rel_dosage:+.2f}" for r in P3.itertuples()))

pd.DataFrame(TESTS).to_csv(f"{OUT}/tests_step5.tsv", sep="\t", index=False)
say(f"\nwrote {OUT}")
LOG.close()
