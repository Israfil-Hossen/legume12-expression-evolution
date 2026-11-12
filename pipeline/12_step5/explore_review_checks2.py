#!/usr/bin/env python3
"""EXPLORATORY (post hoc, prompted by a pre-submission review; not in any specification). Committed results only.

R4  Secondary test with matched orthogroup sets: rho_within recomputed on each species pair's own (equal-copy)
    orthogroup set, compared with rho_between on that set; mean over pairs, species bootstrap (Step 2 rules).
R5  Variance-based study share F2 = (d_study^2 - d_within^2) / (d_species^2 - d_within^2), next to the registered
    linear F, per comparison and per species.
R6  Permutation test of the ohnolog-asymmetry trend (2B): species labels of the divergence-time matrix permuted
    (Mantel-type, 10,000 permutations, C. cajan excluded as in the registered test), for all pairs and collinear pairs.
Usage (repository root): python pipeline/12_step5/explore_review_checks2.py
"""
import os, itertools
import numpy as np
import pandas as pd
from scipy import stats

RES, OUT = "results", "results/exploratory_review"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261009)
NBOOT, NPERM = 2000, 10000
SPECIES = ["Gmax", "Pvul", "Vrad", "Vung", "Ccaj", "Mtru", "Cari", "Psat", "Tpra", "Ljap", "Lang", "Adur"]
LOG = open(os.path.join(OUT, "report_review_checks2.txt"), "w", encoding="utf-8")


def say(s=""):
    print(s); LOG.write(s + "\n"); LOG.flush()


def key(a, b):
    return " ".join(sorted([a, b], key=SPECIES.index))


def boot_pairs(sp_list):
    s = list(rng.choice(sp_list, len(sp_list), replace=True))
    return [key(s[i], s[j]) for i, j in itertools.combinations(range(len(s)), 2) if s[i] != s[j]]


# ------------------------------------------------------------------ R4
say("== R4. Secondary test on matched orthogroup sets (EXPLORATORY, post hoc)")
B = pd.read_csv(f"{RES}/results_step2/rho_between_primary.tsv", sep="\t")
Wreg = pd.read_csv(f"{RES}/results_step2/rho_within_primary.tsv", sep="\t").set_index("species").rho
win = {sp: pd.read_csv(f"{RES}/results_step2/div_within_primary_{sp}.tsv.gz", sep="\t").set_index("Orthogroup")
       for sp in SPECIES if os.path.exists(f"{RES}/results_step2/div_within_primary_{sp}.tsv.gz")}
rows = []
for r in B.itertuples():
    a, b = sorted([r.sp_a, r.sp_b], key=SPECIES.index)
    d = pd.read_csv(f"{RES}/results_step2/div_between_primary_{a}_{b}.tsv.gz", sep="\t").set_index("Orthogroup")
    rw = {}
    for s in (a, b):
        if s in win:
            w = win[s].loc[d.index]
            rw[s] = stats.spearmanr(w.copies, w["div"]).correlation
    rows.append(dict(key=key(a, b), sp_a=a, sp_b=b, n_og=len(d), rho_between=r.rho, rho_within_matched=np.mean(list(rw.values())),
                     rho_within_registered=np.mean([Wreg[s] for s in (a, b) if s in Wreg.index])))
R4 = pd.DataFrame(rows)
R4["diff_matched"] = R4.rho_between - R4.rho_within_matched
R4.to_csv(f"{OUT}/R4_matched_within.tsv", sep="\t", index=False)
idx = R4.set_index("key"); spl = SPECIES
bt = [idx.loc[[k for k in boot_pairs(spl) if k in idx.index]].diff_matched.mean() for _ in range(NBOOT)]
est = R4.diff_matched.mean()
pd.DataFrame([dict(statistic="mean_rho_between_minus_matched_within", n_pairs=len(R4), estimate=est, lo=np.quantile(bt, .025),
                   hi=np.quantile(bt, .975), n_negative=int((R4.diff_matched < 0).sum()), median_matched_within=R4.rho_within_matched.median())]).to_csv(
    f"{OUT}/R4_summary.tsv", sep="	", index=False)
say(f"  rho_within on matched sets: median {R4.rho_within_matched.median():+.3f} (registered, all orthogroups: median {Wreg.median():+.3f})")
say(f"  mean rho_between - matched rho_within over {len(R4)} pairs: {est:+.3f} [{np.quantile(bt, .025):+.3f}, {np.quantile(bt, .975):+.3f}]; "
    f"negative in {(R4.diff_matched < 0).sum()}/{len(R4)} pairs")

# ------------------------------------------------------------------ R5
say("\n== R5. Variance-based study share F2 (EXPLORATORY, post hoc)")
A = pd.read_csv(f"{RES}/results_step5/1A_pairs.tsv", sep="\t")
A["F2"] = (A.d_study ** 2 - A.d_within ** 2) / (A.d_species ** 2 - A.d_within ** 2)
A.to_csv(f"{OUT}/R5_F2.tsv", sep="\t", index=False)
for sp, g in A.groupby("focal", sort=False):
    say(f"  {sp}: F (registered, linear) {g.F.mean():.2f}; F2 (variance) {g.F2.mean():.2f}")
say(f"  all 33: F {A.groupby('focal').F.mean().mean():.2f}; F2 {A.groupby('focal').F2.mean().mean():.2f}")

# ------------------------------------------------------------------ R6
say("\n== R6. Permutation test of the ohnolog-asymmetry trend (EXPLORATORY, post hoc)")
P = pd.read_csv(f"{RES}/results_step5/2B_species_pairs.tsv", sep="\t")
P = P[(P.sp_a != "Ccaj") & (P.sp_b != "Ccaj")].copy()
sp11 = [s for s in SPECIES if s != "Ccaj"]
T = {key(r.sp_a, r.sp_b): r.time for r in P.itertuples()}
perm_rows = []
for col, lab in (("rho_r", "all ohnolog pairs"), ("rho_r_col", "collinear pairs")):
    d = P[np.isfinite(P[col])]
    rngT = d.time.max() - d.time.min()
    obs = np.polyfit(d.time, d[col], 1)[0] * rngT
    perm = []
    for _ in range(NPERM):
        mp = dict(zip(sp11, rng.permutation(sp11)))
        tt = np.array([T[key(mp[a], mp[b])] for a, b in zip(d.sp_a, d.sp_b)])
        perm.append(np.polyfit(tt, d[col], 1)[0] * rngT)
    p = (np.sum(np.abs(perm) >= abs(obs)) + 1) / (NPERM + 1)
    say(f"  {lab} ({len(d)} species pairs): change {obs:+.3f}; permutation P = {p:.4f} (two-sided, {NPERM} species-label permutations)")
    perm_rows.append(dict(measure=col, n_pairs=len(d), change=obs, p_perm=p, n_perm=NPERM))
pd.DataFrame(perm_rows).to_csv(f"{OUT}/R6_permutation.tsv", sep="	", index=False)
