#!/usr/bin/env python3
"""EXPLORATORY (post hoc, not in docs/step5_spec.md): robustness of the 2B result (erosion of ohnolog
asymmetry with divergence time). Reads results/results_step5; writes results/exploratory_step5/.
  a. drop G. max and L. angustifolius (lineage-specific polyploidy; fewest ohnolog pairs)
  b. drop comparisons with fewer than 100 shared ohnolog pairs
  c. attenuation-corrected conservation: rho_r / sqrt(replicate ceiling A x replicate ceiling B)
  d. joint regression: rho_r ~ time + mean replicate ceiling of the two species
"""
import os, numpy as np, pandas as pd
from scipy import stats
OUT = "results/exploratory_step5"; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261010); NBOOT = 2000
P = pd.read_csv("results/results_step5/2B_species_pairs.tsv", sep="\t")
CE = pd.read_csv("results/results_step5/2B_ceilings.tsv", sep="\t")
rep = CE[CE.ceiling == "replicate"].set_index("species").rho
P = P[(P.sp_a != "Ccaj") & (P.sp_b != "Ccaj")].reset_index(drop=True)


def trend(df, y):
    rx = df.time.max() - df.time.min()
    fit = lambda d: np.polyfit(d.time, d[y], 1)[0] * rx if len(d) >= 3 and d.time.var() > 0 else np.nan
    sp = sorted(set(df.sp_a) | set(df.sp_b)); key = {frozenset(p): i for i, p in enumerate(zip(df.sp_a, df.sp_b))}
    bt = []
    for _ in range(NBOOT):
        s = list(rng.choice(sp, len(sp))); idx = [key[frozenset((s[i], s[j]))] for i in range(len(s)) for j in range(i + 1, len(s))
                                                if s[i] != s[j] and frozenset((s[i], s[j])) in key]
        bt.append(fit(df.iloc[idx]) if idx else np.nan)
    return fit(df), np.nanquantile(bt, .025), np.nanquantile(bt, .975), len(df)


out = []
def show(label, r):
    out.append(dict(check=label, estimate=r[0], lo=r[1], hi=r[2], n_pairs=r[3]))
    print(f"  {label:58s} change {r[0]:+.3f} [{r[1]:+.3f}, {r[2]:+.3f}]  ({r[3]} comparisons)")

print("EXPLORATORY robustness of 2B (registered estimate -0.224 [-0.281, -0.147])")
show("a. without G. max and L. angustifolius", trend(P[~P.sp_a.isin(["Gmax", "Lang"]) & ~P.sp_b.isin(["Gmax", "Lang"])].reset_index(drop=True), "rho_r"))
show("b. only comparisons with >= 100 shared ohnolog pairs", trend(P[P.n >= 100].reset_index(drop=True), "rho_r"))
Q = P[P.sp_a.isin(rep.index) & P.sp_b.isin(rep.index)].copy()
Q["rho_r_corr"] = Q.rho_r / np.sqrt(rep[Q.sp_a].values * rep[Q.sp_b].values)
show("c. attenuation-corrected rho_r (Lotus has no ceiling)", trend(Q.reset_index(drop=True), "rho_r_corr"))
Q["mean_ceiling"] = (rep[Q.sp_a].values + rep[Q.sp_b].values) / 2
X = np.column_stack([np.ones(len(Q)), Q.time, Q.mean_ceiling]); beta, *_ = np.linalg.lstsq(X, Q.rho_r, rcond=None)
rx = Q.time.max() - Q.time.min()
print(f"  d. joint OLS rho_r ~ time + mean ceiling: time effect over {rx:.1f} Myr {beta[1] * rx:+.3f}; ceiling coefficient {beta[2]:+.3f}"
      f";  Spearman(rho_r, mean ceiling) {stats.spearmanr(Q.rho_r, Q.mean_ceiling).correlation:+.3f}")
out.append(dict(check="d. joint OLS time effect (ceiling held constant)", estimate=beta[1] * rx, lo=np.nan, hi=np.nan, n_pairs=len(Q)))
print(f"  closest pair: {P.loc[P.time.idxmin(), 'sp_a']}-{P.loc[P.time.idxmin(), 'sp_b']} {P.time.min():.1f} Myr rho_r {P.loc[P.time.idxmin(), 'rho_r']:+.3f}; "
      f"laboratory ceilings {', '.join(f'{r.species} {r.rho:+.3f}' for r in CE[CE.ceiling == 'laboratory'].itertuples())}")
pd.DataFrame(out).to_csv(f"{OUT}/2B_robustness.tsv", sep="\t", index=False)
