# Step 4 specification — addendum: two studies of one species (the "between-study" rung)

**Status: FINAL, 2026-09-24.** Frozen for registration; the author's registration of
this file on OSF is its approval. This file is not edited from here on; later
changes are recorded as dated deviations at the end.

**Timing, stated plainly.** This addendum was written **after** the Step 2 results and
the exploratory analyses E1–E4 had been seen. Its motivation is therefore post hoc:
exploratory analysis E2 showed that between-species ρ is about 3.8 times what
within-species replicate noise alone predicts. The analysis below is specified
**before** any Step 4 data are downloaded or quantified.

---

## 1. Question

In the Step 2 design every species comes from a different study, so a difference
between two species is also a difference between two laboratories, growth conditions
and protocols. Is the copy-number gradient between species larger than the gradient
between two independent studies of the *same* species, where no evolutionary
divergence exists?

## 2. Predictions

| Account | Prediction |
|---|---|
| **Study (technical)** | Between-species ρ is no larger than between-study ρ on the same orthogroups |
| **Species (evolutionary not excluded)** | Between-species ρ exceeds between-study ρ by a material amount |

## 3. Data

A second leaf study for *Glycine max*, *Phaseolus vulgaris* and *Medicago truncatula*
(the species with the most public leaf data). Study 1 is the Step 2 study.
Exact runs: `metadata/step4_runs.tsv` (33 biological replicates, 33 runs).

| Species | Study 2 | Control groups × biological replicates | Tissue (BioSample) |
|---|---|---|---|
| Gmax | PRJNA877350 | 3 mock groups × 3 | leaves, R3–R5 |
| Pvul | PRJNA741786, genotype "Resistant" | 4 mock time points × 3 | leaf |
| Mtru | PRJNA975868, accession A17 | 4 mock time points × 3 (library series SM) | leaf |

**Selection rule, applied to all verified candidates:** leaf in the BioSample record;
unstressed or mock controls; at least 4 biological replicates with distinct libraries;
at most 4 control groups within one genotype (if a study has several genotypes, the
first by name is used, as for *Lupinus* in Step 2); among those, the study with the
most control replicates, ties broken by read depth. Candidates rejected on
verification: PRJNA877710 (BioSample tissue given as "Anther; Leaf samples"),
PRJNA210619 and PRJNA471093 (no replicates), PRJNA900512 (design not recoverable).
*Medicago*: the three 72 hpi libraries of the "RM" series are excluded because their
prefix places them outside the SM mock series used at every other time point.

## 4. Processing

Identical to Step 2 (`docs/step2_spec.md` section 5): same fastp and salmon settings,
same indices, same QC thresholds, same orthogroups and copy numbers (the Step 2
OrthoFinder run; these species' proteomes are unchanged), same scaling and log2(x + 1).
2 samples per side throughout.

## 5. Measurements

For focal species A (Gmax, Pvul, Mtru) and every other species B (11 per focal species):
- **Orthogroup set:** the Step 2 set for pair (A, B): orthogroups with equal copy number in A and B.
- **ρ_species(A, B):** the Step 2 between-species value on that set.
- **ρ_study(A | set):** Spearman ρ between copy number and between-study divergence
  of A on the same set; divergence = |mean of 2 study-1 samples − mean of 2 study-2
  samples| / √2, averaged over 50 draws.
- **Δ(A, B) = ρ_species(A, B) − ρ_study(A | set).** Identical orthogroups, identical statistic.
- **Δ** = mean over B for each A, then mean over the three A.

95% intervals: 2,000 bootstrap resamples of orthogroups within each pair, the
per-orthogroup divergences held fixed.

## 6. Decision rule

Floors as before: material ≥ 0.050, negligible ≤ 0.025.
- Δ ≥ 0.050 and interval excluding 0 → **between-species gradient exceeds between-study** (species-level component; evolutionary account not excluded)
- |Δ| ≤ 0.025 → **between-species gradient no larger than between-study** (study account)
- otherwise → **partial**; report the magnitude

Per-species Δ values are reported alongside.

## 7. Also reported (descriptive, no decision rule)

- ρ_study on all orthogroups, against the half-split ρ of study 1 (Step 2) and of study 2
- The attenuation reference (observed/expected, the manuscript's method) for the study rung, beside the same quantity for the species rung

## 8. Reporting

Every result above is reported whatever it shows. Any departure is recorded below as a
dated deviation.

---

## Deviations

*(none yet)*
