# Step 3 specification — addendum: duplicate age, the paralog-level design, expression-evolution models, orthology robustness

**Status: FINAL, 2026-09-23.** Frozen for registration; the author's registration of
this file on OSF is its approval. Written before any Step 2 result existed: no
expression values from the 12-species data had been analysed, and only sample-level
QC (salmon mapping rates of two test samples) had been inspected. This file is not
edited from here on; later changes are recorded as dated deviations at the end.

Uses the Step 2 data and processing unchanged (`docs/step2_spec.md`, sections 3–5):
same 12 species, same QC-passed samples, same salmon quantification, same
single-copy-orthogroup scaling. Magnitude floors are the manuscript's:
**material ≥ 0.050, negligible ≤ 0.025** on the Spearman scale.

---

## A. The paralog-level design: does expression divergence between paralogs grow with duplicate age?

**Why.** The manuscript argues that orthogroup-level cross-species comparisons
cannot show whether duplication relaxes constraint. The constructive question is
whether the valid design (paralogs compared within one species, in the same
samples, at gene level) can. Gene-level comparison has its own technical hazard:
near-identical young paralogs share reads, so their individual estimates are
noisy. This analysis measures that noise directly instead of assuming it away.

**Paralog pairs and age (per species).**
- Proteins: longest isoform per gene (the OrthoFinder input). CDS: RefSeq `cds_from_genomic`.
- DIAMOND all-vs-all within species (e-value ≤ 1e-10, up to 5 hits per gene).
- MCScanX collinearity with default parameters; its `duplicate_gene_classifier`
  assigns each gene to singleton, dispersed, proximal, tandem or WGD/segmental.
- Paralog pairs: reciprocal best within-species DIAMOND hits.
- Ks: protein alignment (MAFFT) back-translated to codons (pal2nal), then `yn00` (PAML).
  Pairs with Ks > 3 are excluded (saturation). Low-Ks pairs are kept; they are the
  pairs the ambiguity hazard concerns.

**Expression.** Gene TPM per sample (salmon + transcript→gene map), scaled as in
Step 2, log2(x + 1). Pairs in which both genes have mean TPM < 1 are excluded
(not expressed in leaf).

**Measures per pair.**
- **D**, expression divergence: |mean log2 expression of gene 1 − of gene 2| over all QC-passed samples.
- **N**, noise floor of that contrast: for 50 stratified half-splits (2 vs 2, as in
  Step 2), |(gene1 − gene2)<sub>half 1</sub> − (gene1 − gene2)<sub>half 2</sub>| / √2, averaged.
- **D − N**, divergence in excess of measurement noise.

**Predictions and decision rules.** Spearman ρ against Ks within each species; the
summary is the median over species, with a 95% interval from resampling species
2,000 times.

| | Prediction | Decision |
|---|---|---|
| A1 | ρ(D, Ks) > 0 (classic duplicate divergence) | reported |
| A2 | ρ(N, Ks) < 0 (young paralogs noisier: read-assignment ambiguity) | median ρ ≤ −0.050 with interval excluding 0 → **ambiguity noise confirmed** |
| **A3 (primary)** | ρ(D − N, Ks) > 0 | median ρ ≥ 0.050, interval excluding 0 → **evolutionary divergence detectable at gene level**; ≤ 0.025 → **not detectable**; otherwise **partial** |

**Exploratory, no decision rule:** WGD/segmental versus other duplicate classes,
compared within Ks bins.

## B. Expression-evolution model with and without within-species variance

**Why.** Phylogenetic models of expression evolution estimate an evolutionary rate
from species means. If part of the between-species variance is measurement noise
that scales with family size, the fitted rate will too. Models that include
within-species variance (e.g. EVE) should remove that part.

**Data.** Orthogroups present in all 12 legumes. Per species: mean per-copy
expression (orthogroup sum ÷ copies, log2) over QC-passed samples, and its squared
standard error (within-species variance ÷ n). Tree: the Step 2 IQ-TREE tree,
*A. thaliana* pruned. Predictor: mean log2 copy number over the 12 species.

**Model.** Brownian motion fitted by maximum likelihood per orthogroup, twice:
(i) σ²<sub>naive</sub>, species means only; (ii) σ²<sub>me</sub>, with each species'
squared standard error added to the diagonal of the covariance matrix.

| | Prediction | Decision |
|---|---|---|
| B1 | ρ(σ²<sub>naive</sub>, copies) > 0 (the phenomenon, in model form) | reported |
| **B2 (primary)** | Δρ = ρ(σ²<sub>me</sub>, copies) − ρ(σ²<sub>naive</sub>, copies), 95% interval from 2,000 bootstraps over orthogroups | Δρ ≤ −0.050 → **within-species variance accounts for a material part of the gradient** (technical); \|Δρ\| ≤ 0.025 → **no change** (evolutionary); otherwise **partial** |

## C. Orthology robustness

Repeat the Step 2 primary test with OrthoFinder's hierarchical orthogroups at the
root node (`Phylogenetic_Hierarchical_Orthogroups/N0.tsv`) in place of orthogroups.
Everything else unchanged. Reported as agreement or disagreement with the Step 2
primary verdict; no separate decision rule.

## D. Reporting

Every result is reported whatever it shows. Analyses A–C are confirmatory only in
the sense defined here; anything else is labelled exploratory in the manuscript.

---

## Deviations

### D2 — 2026-09-24: A3 statistic replaced (the registered one is biased)

Recorded before any Step 3 data were analysed: paralog pairs and Ks had not been
computed, and the analysis code had been run only on synthetic data.

While validating the code on synthetic data, the registered A3 statistic (Spearman
ρ between D − N and Ks) turned out to be biased. D is averaged over all samples and N
over 2 + 2 half-splits, so they carry different amounts of noise. When noise differs
between pairs, which is exactly what A2 predicts, D − N varies with noise alone. In
simulations with no true divergence it gave ρ = +0.44, a spurious "evolutionary"
result (pipeline/09_step3/validation/a3_bias_sim.R, a3_fix_sim.R).

Change: for each pair and each of the same 50 stratified half-splits, the product of
the two half-means of the paralog difference, m₁·m₂, is an unbiased estimate of the
squared true divergence δ², whatever the noise level. Its average over draws is
related to Ks by **Pearson** correlation, because rank correlation is itself biased by
noise that differs between pairs. In the same simulations: no divergence r = +0.003;
equal divergence at all ages r = +0.000; divergence growing with Ks r = +0.891.

A3 (primary) is now: Pearson r between mean m₁·m₂ and Ks, per species; median over
species with the species bootstrap; the same floors and decision wording as before.
The registered statistic is still computed and reported, labelled as biased. A1 and
A2 are unchanged.

### D3 — 2026-09-25: analysis C uses the legume crown node, not N0

Recorded before analysis C was run on the real data (its only attempt stopped at
start-up because the input file was missing). The paralog jobs for analysis A had
run, but no A, B or C statistic had been computed.

Section C assumed OrthoFinder 2 output, where `Orthogroups.tsv` (MCL clusters) and
`Phylogenetic_Hierarchical_Orthogroups/N0.tsv` (root HOGs) are two different
orthology definitions. The Step 2 run used OrthoFinder 3.1.5 (log of `Results_Sep23`).
From version 3.1 OrthoFinder no longer writes `N0.tsv`; `Orthogroups/Orthogroups.tsv`
holds the N0 HOGs instead (OrthoFinder release notes). The Step 2 primary analysis
therefore already used root HOGs, and C as registered would repeat it exactly.

Change: C uses the HOGs of the legume crown node, the node below the root that
contains all 12 legumes and not *A. thaliana*. The script picks the one node file in
`Phylogenetic_Hierarchical_Orthogroups/` whose species columns are exactly the 12
legumes and stops if there is not exactly one. With *A. thaliana* as outgroup, this is
the finest HOG level at which every species pair in the study is still compared
within one group, so it tests whether the verdict depends on the orthology level.
Everything else in C is unchanged. The registered N0 version is reported as identical
to the Step 2 primary analysis by construction. Checked on synthetic data laid out
like OrthoFinder 3.1.5 output: default mode is identical to before, and C gives the
same result as the earlier N0-based test.

Implementation note (same day, before C produced any result): the first version of the
rule looked for a node file whose species columns were exactly the 12 legumes. In
OrthoFinder 3.1.5, every node file keeps all 13 species columns (empty outside the
clade), so the run stopped at start-up with "found 0". The node is now taken from
OrthoFinder's rooted species tree (`Species_Tree/SpeciesTree_rooted_node_labels.txt`):
the one node whose descendants are exactly the 12 legumes. In this run the root N0
separates *A. thaliana* from the legumes, so the node is N1. As a check, the report
gives the number of HOGs in that file that contain *A. thaliana* genes; it should be 0.
