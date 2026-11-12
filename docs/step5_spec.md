# Step 5 specification — laboratory versus evolution, and the ancient WGD seen from within samples

**Status: FINAL, 2026-10-09.** Written and committed before any of the analyses below
were run. Later changes are recorded as dated deviations at the end; the text above
them is not edited.

## What had been seen when this was written

Stated so that a reader can judge what could have influenced these choices.

- All registered results of Steps 2–4 (`results/results_step2`, `results_step3`,
  `results_step3C`, `results_step4`) and the post hoc exploratory analyses E1–E4
  (`results/exploratory_step2`). E1 included the median between-species divergence
  against time (Spearman +0.063).
- The Step 3 exploratory table (median D − N by Ks bin and collinearity), and on
  2026-10-09 a quick per-species version of it: in all 11 species with a half-split,
  collinear paralog pairs at Ks 0.1–1 had lower expression divergence D than
  non-collinear pairs.
- On 2026-10-09, the structure of ancient duplicate pairs only: 587 orthogroups split
  into exactly two legume-crown (N1) hierarchical orthogroups that are both single-copy
  in at least 10 species; per species about 550 such pairs, 60–80% collinear, median Ks
  0.64–1.15. **No expression value of these pairs was looked at.**
- Not seen: any expression ratio between ohnologs, any expression-level-controlled
  gradient, any read-ambiguity value, any dosage computation, any between-study
  divergence other than the copy-number gradient of Step 4.

## Data and shared definitions

Same data and processing as Steps 2 and 4: the 108 QC-passed Step 2 samples (study 1)
and the 32 QC-passed Step 4 samples (study 2 of *G. max*, *P. vulgaris*, *M.
truncatula*); the Step 2 OrthoFinder run; per-sample scaling to the grand median of
the Step 2 single-copy-orthogroup medians (constant C, as in Step 4); log2(x + 1).
Divergence between two sets of samples = |difference of means| / √2 with 2 samples per
side and 50 draws (half-splits stratified by group, as in Step 2). Intervals: 2,000
bootstrap resamples of the stated unit. Floors on the Spearman scale as before:
material ≥ 0.050, negligible ≤ 0.025. Cajanus is excluded from divergence-time
regressions (deviation D1 of Step 2); species bootstrap as in Step 2.

---

## Part 1 — How much of cross-species expression divergence is laboratory?

### 1A. Expression distance: within study, between studies, between species

For each focal species A ∈ {Gmax, Pvul, Mtru} and each other species B, on the Step 2
orthogroup set of pair (A, B) (equal copy number), take the median over orthogroups of:
d_within(A) (half-split within study 1), d_study(A) (study 1 vs study 2) and
d_species(A, B) (the Step 2 between-species divergence).

**Primary statistic:** the laboratory share
F(A, B) = [d_study(A) − d_within(A)] / [d_species(A, B) − d_within(A)],
the fraction of the between-species excess over replicate noise that is reproduced by
changing laboratory within species. F = mean over B, then over the three A. Interval:
2,000 bootstraps of orthogroups within each pair (divergences held fixed), as Step 4.

| F | Verdict |
|---|---|
| ≥ 0.50 and lower bound > 0.25 | laboratory accounts for at least half of the between-species excess |
| ≤ 0.25 | evolution dominates the between-species excess |
| otherwise | partial |

**Secondary:** across the 55 pairs (Cajanus excluded), regress noise-corrected
between-species distance, d_species(A, B) − [d_within(A) + d_within(B)]/2, on TimeTree
divergence time; report the predicted change over the observed range in log2 units,
species bootstrap. Lotus has no half-split and is excluded from this regression.
Reported, no decision rule.

### 1B. Does the copy-number gradient survive control for expression level?

For every within-species half-split (11 species) and every between-species pair (66),
recompute the gradient as a **partial Spearman** correlation between copy number and
divergence, controlling for the orthogroup's mean log2 expression in the compared
samples. Secondary covariate: mean per-copy expression (log2 of sum ÷ copies).

**Primary statistic:** median partial ρ over the 11 within-species comparisons, with a
species-bootstrap interval.

| Median partial ρ | Verdict |
|---|---|
| ≥ 0.050 and interval excluding 0 | the gradient persists after control for expression level |
| ≤ 0.025 | the gradient is explained by expression level |
| otherwise | partial |

Also reported: the same for the 66 between-species pairs, and the median gradient
within expression deciles.

### 1C. Does read ambiguity explain the gradient?

From salmon's `aux_info/ambig_info.tsv` (UniqueCount, AmbigCount per transcript),
ambiguity of an orthogroup in a sample = Σ AmbigCount / Σ (UniqueCount + AmbigCount)
over the transcripts of its genes; per species, the median over that species'
QC-passed samples. Reported: ρ(copy number, ambiguity) per species.

**Primary statistic:** median over the 11 within-species comparisons of the partial
Spearman ρ between copy number and divergence controlling for **both** ambiguity and
mean log2 expression; same verdict table as 1B (wording: "persists after control for
expression level and read ambiguity" / "is explained by expression level and read
ambiguity" / partial). Also reported: partial ρ(ambiguity, divergence | copy number).
If the ambiguity files cannot be obtained, 1C is reported as not done.

---

## Part 2 — The papilionoid whole-genome duplication seen from within samples

### Pairs

Ohnolog orthogroups: orthogroups that OrthoFinder splits into exactly two N1
(legume-crown) hierarchical orthogroups, both single-copy in at least 10 of the 12
species. In each species, a pair is used if both HOGs have exactly one gene there and
the mean scaled TPM of at least one gene is ≥ 1. Copy A is the HOG with the smaller
identifier; the labelling is therefore the same in every species. Primary set: all
such pairs; sensitivity: pairs collinear (MCScanX) in the species concerned.

Ratio in one sample: r = log2(x_A + 1) − log2(x_B + 1), x = scaled TPM.

### 2A. Are within-sample ratios protected from the laboratory effect?

For A ∈ {Gmax, Pvul, Mtru}: for ohnolog ratios and, on the same genes, for absolute
log2 expression (mean of the two genes' divergences), compute median d_within (study
1 half-split) and median d_study (study 1 vs study 2). Laboratory excess
X = median d_study − median d_within.

**Primary statistic:** E = X_ratio / X_absolute, mean over the three species; interval
from 2,000 bootstraps of ohnolog pairs.

| E | Verdict |
|---|---|
| ≤ 0.50 and upper bound < 1 | within-sample ratios largely remove the laboratory effect |
| ≥ 0.90 | within-sample ratios do not remove the laboratory effect |
| otherwise | partial |

### 2B. Is ohnolog expression asymmetry conserved, and does it erode with time?

For each species pair, conservation of asymmetry ρ_r(A, B) = Spearman correlation,
over ohnolog pairs present in both species, between the species' mean ratios. Same
quantity between the two halves of one study (replicate ceiling) and between the two
studies of Gmax, Pvul and Mtru (laboratory ceiling).

**Primary test:** regression of ρ_r(A, B) on TimeTree divergence time over the 55 pairs
without Cajanus; statistic = predicted change in ρ_r over the observed range; species
bootstrap.

| Change | Verdict |
|---|---|
| ≤ −0.050 and interval excluding 0 | asymmetry erodes with divergence time (evolution detectable) |
| absolute change ≤ 0.025 | asymmetry does not change with divergence time |
| otherwise | partial |

**Also reported (no decision rule):** median ρ_r between species versus the replicate
and laboratory ceilings; the same regression with the collinear-only set and on
patristic distance (IQ-TREE tree, Cajanus included); for each ohnolog pair, the number
of species in which copy A is higher (binomial test against 0.5 per pair, with the
proportion of pairs significant at FDR 5%).

### 2C. Which ohnolog pairs are asymmetric? (descriptive)

Pairs classified as consistently asymmetric (same copy higher in ≥ 80% of species where
present, mean |r| ≥ 1) or balanced (mean |r| < 0.5 in ≥ 80% of species). Compared on
Ks, collinearity and, if *A. thaliana* product names can be obtained for the
orthogroup, functional categories (ribosomal proteins, transcription factors, other),
descriptive only.

---

## Part 3 — Relative dosage after lineage-specific duplication

Orthogroups in which exactly one species X has two copies and every other species
present has one copy. For X ∈ {Gmax (Glycine WGD: pair collinear and Ks < 0.3), and,
descriptively, every species with tandem pairs (MCScanX tandem/proximal)}: relative
dosage D_X = scaled log2 family total in X minus the mean of the same quantity in the
single-copy species, per orthogroup. Calibration (no copy-number difference, same
studies): the same D_Gmax computed on orthogroups single-copy in all 12 species. This
removes the overall offset between the Gmax study and the other studies.

**Primary statistic:** median D_Gmax over Glycine-WGD orthogroups minus the median of
the calibration distribution; interval from 2,000 bootstraps of orthogroups (both sets
resampled). Replication: the same statistic with Gmax study 2 in place of study 1.

| Value (log2 units) | Verdict |
|---|---|
| ≥ 0.75 | relative family expression rises with copy number (dosage retained) |
| ≤ 0.25 | relative family expression stays constant (dosage sharing) |
| otherwise | partial |

Limitation stated in advance: TPM-based measures are relative to each library's
transcript pool; they cannot detect a change in absolute transcriptome size.

---

## Reporting

Every result is reported whatever it shows, with the verdict wording above. Any
departure is recorded below as a dated deviation. Analyses not listed here are labelled
exploratory.

---

## Deviations

*(none yet)*
