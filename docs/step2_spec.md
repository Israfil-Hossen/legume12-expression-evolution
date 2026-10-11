# Step 2 specification — copy number and expression divergence across 12 legumes

**Status: APPROVED by the author, 2026-09-22.** Written and committed before any
expression data for the new species had been downloaded or quantified. This file is
not edited from here on; any later change is recorded as a dated deviation at the
end, never by rewriting the text above it.

---

## 1. Question

The manuscript shows, in three species, that the copy-number gradient in expression
divergence is already fully present between two halves of one experiment, where
nothing evolutionary differs. The objection it cannot answer with three species is
whether the gradient *grows with evolutionary distance* once measurement noise is
held fixed. This step asks exactly that, across 12 legume species spanning the five
main Papilionoideae clades.

## 2. Predictions (fixed now)

| Account | Prediction |
|---|---|
| **Evolutionary** — duplication relaxes constraint on expression | The copy-number gradient between two species increases with their divergence time, and exceeds the within-species half-split gradient |
| **Technical** — the gradient is a property of measurement | The between-species gradient does not change with divergence time and is close to the within-species half-split gradient |

## 3. Species and genomes

NCBI RefSeq annotations only, the rule used in the manuscript. Diploids only.

| Code | Species | Clade | Accession |
|---|---|---|---|
| Gmax | *Glycine max* | Millettioid | GCF_000004515.6 |
| Pvul | *Phaseolus vulgaris* | Millettioid | GCF_000499845.2 |
| Vrad | *Vigna radiata* | Millettioid | GCF_000741045.1 |
| Vung | *Vigna unguiculata* | Millettioid | GCF_004118075.2 |
| Ccaj | *Cajanus cajan* | Millettioid | GCF_000230855.2 |
| Mtru | *Medicago truncatula* | IRLC | GCF_003473485.1 |
| Cari | *Cicer arietinum* | IRLC | GCF_000331145.2 |
| Psat | *Pisum sativum* | IRLC | GCF_024323335.1 |
| Tpra | *Trifolium pratense* | IRLC | GCF_020283565.1 |
| Ljap | *Lotus japonicus* | Robinioid | GCF_012489685.1 |
| Lang | *Lupinus angustifolius* | Genistoid | GCF_001865875.1 |
| Adur | *Arachis duranensis* | Dalbergioid | GCF_000817695.3 |
| Atha | *Arabidopsis thaliana* | outgroup (proteins only) | GCF_000001735.4 |

Excluded, with reason: *Lens culinaris* (no annotated assembly), *Lupinus albus* and
*Vicia faba* (GenBank annotation only), *Arachis hypogaea* (allotetraploid).

## 4. Datasets

**Design: control (unstressed) leaf expression**, one study per species. The exact
runs are listed in `metadata/step2_runs.tsv` (109 biological replicates, 118 runs).
All 12 tissues are confirmed as leaf in the BioSample records.

| Species | Study | Groups × biological replicates | Median reads |
|---|---|---|---|
| Gmax | PRJNA186522 | 4 lines × 3 | 75M |
| Pvul | PRJNA911667 | 3 time points × 3 | 49M |
| Mtru | PRJNA701852 | 4 time points × 3 | 26M |
| Vrad | PRJNA764584 | 4 groups × 3 | 22M |
| Vung | PRJNA743032 | 2 time points × 5 | 46M |
| Ccaj | PRJNA1250142 | 2 genotypes × 3 | 18M |
| Cari | PRJNA436616 | 2 genotypes × 4 | 50M |
| Psat | PRJNA850002 | 4 groups × 3 | **6M** (flagged) |
| Tpra | PRJNA1091176 | 2 months × 5 | 23M |
| Lang | PRJEB105685 | 2 lines × 5 | 85M |
| Adur | PRJNA726938 | 1 × 5 | 34M |
| Ljap | PRJDB9534 | 1 × **3** (flagged) | 31M |

**Selection rules.**
- Controls must be biological replicates of one genotype within a group. Panels of
  different accessions are not replicates.
- A biological replicate is a separate library (distinct `library_name`) or GEO
  sample. Runs from the same library are lanes and are merged before quantification.
- Standard mRNA-seq only. Degradome, small-RNA and methylation libraries are excluded.
- *Lupinus*: lines 83A-476 and 95826 were chosen by name order, before any data were
  seen, to match the two-group design of the other species.

## 5. Processing — identical to the manuscript

- fastp: `-q 20 -u 30 -l 36 --trim_poly_g`, plus `--detect_adapter_for_pe` for paired-end
- salmon: decoy-aware index (genome as decoy), `k = 31`; quant `-l A --seqBias --gcBias --numBootstraps 30`
- tximport, `countsFromAbundance = "no"`; transcript→gene from the GFF, keyed on NCBI GeneID
- OrthoFinder rerun on all 12 legumes + *A. thaliana*, longest isoform per gene.
  The new orthogroups replace the manuscript's for this analysis only.

**Sample QC (fixed thresholds).** Exclude a sample if salmon mapping rate < 0.70, if
fewer than 5 million reads survive fastp, or if more than 50% of its reads are
assigned to rRNA genes (the chickpea lesson). If exclusions leave a species below
4 replicates it drops out of the half-split. Below 2, it drops out entirely.

## 6. Measurements

**Orthogroup expression.** Per sample, sum gene TPM over paralogs within the
orthogroup. Scale each sample so that the median over orthogroups single-copy in all
12 species is equal across samples. Use log2(x + 1).

**Divergence** between two sets of samples = |difference in mean log2 expression| / √2,
the manuscript's statistic.

**Fixed sample size.** Every comparison uses **2 samples per side**, so that no
comparison is noisier than another by construction:
- *Within species (half-split):* 2 vs 2 samples of one species, drawn stratified by
  group so the halves have the same group composition. 50 draws, averaged per
  orthogroup. Needs ≥ 4 replicates, so 11 species (not *Lotus*).
- *Between species:* 2 samples of species A vs 2 of species B. 50 draws. All 66 pairs.

**Copy-number confound in expression levels.** A sum over paralogs rises with copy
number, so between-species level differences inflate wherever two species differ in
copy number. **Primary:** each pairwise comparison uses only orthogroups with the
same copy number in both species, and that shared number is the predictor.
**Sensitivity:** per-copy expression (orthogroup sum ÷ copies) on all orthogroups.

**Gradient** = Spearman ρ between copy number and divergence, per comparison: 11
within-species values (ρ_within) and 66 between-species values (ρ_between).

**Divergence time.** Primary: TimeTree median estimates for each species pair, to be
fetched and frozen in an appendix to this file **before** expression data are
analysed. Sensitivity: patristic distance on the concatenated single-copy-ortholog
tree (IQ-TREE).

## 7. Tests and decision rules

Magnitude floors are the manuscript's: **material ≥ 0.050, negligible ≤ 0.025** on
the Spearman scale.

**Primary test.** Regress ρ_between on divergence time across the 66 pairs. The
quantity is the predicted change in ρ between the shortest and longest divergence
times observed. Pairs share species, so the 95% interval comes from resampling
species (not pairs) with replacement, 2,000 times.
- Predicted change ≥ 0.050 and interval excludes zero → **supports evolutionary account**
- Predicted change ≤ 0.025 → **supports technical account**
- In between → **partial**; report the magnitude and do not claim either

**Secondary test.** Mean ρ_between minus mean ρ_within, with the same species-level
bootstrap and the same floors.

## 8. Sensitivity analyses (declared now; all reported)

1. Per-copy expression instead of the equal-copy restriction
2. 3 samples per side, using only species with ≥ 6 replicates
3. Excluding *Pisum* (low read depth)
4. Reads thinned to a common depth of 5 million per sample
5. Patristic distance instead of TimeTree divergence time

## 9. Reporting

Every result above is reported whatever it shows. No analysis is dropped because of
its outcome. Any departure from this file is added below as a dated deviation with
its reason.

## 10. Limitations known in advance

- Species is confounded with laboratory, growth conditions and developmental stage
  (leaf ages range from 2 weeks to 6 months), as in the manuscript.
- *Pisum* has low read depth. *Lotus* has only 3 replicates and has no half-split.
- Orthogroups are rebuilt with 13 proteomes, so copy numbers differ from the
  manuscript's three-species orthogroups.

---

## Deviations

### D1 — 2026-09-23: *Cajanus* removed from the primary divergence-time regression

Recorded before any expression values were analysed. The only quantification output
inspected so far is the salmon mapping rate of two test samples, as a pipeline check.

TimeTree places *Cajanus cajan* outside the Millettioid clade. Its median age to
*G. max* (49.5 Myr) equals its age to *M. truncatula* and *L. japonicus*, and 25 of
the 220 species triples are inconsistent with any single tree, every one of them
involving *Cajanus*. TimeTree's adjusted (synthesis-tree) ages are tree-consistent but
place *Cajanus* the same way (53.9 Myr to every Millettioid and Hologalegina
species). *Cajanus* belongs to Phaseoleae with *Glycine*, *Phaseolus* and *Vigna*, so
both sources misplace it.

Change: the primary regression (section 7) uses the 55 pairs among the other 11
species. *Cajanus* remains in the within-species half-split, in the secondary test,
and in sensitivity analysis 5 (patristic distances on our own species tree), where
its position is taken from the data. No divergence time was chosen by hand.


---

## Appendix A — divergence times (frozen 2026-09-23)

Fetched from the TimeTree pairwise API on 2026-09-23: the precomputed median age
across published studies, with its interval and number of studies, plus the
synthesis-tree adjusted age. All 66 pairs are in `metadata/step2_divergence_times.tsv`,
committed together with this appendix. Median ages range from 9.7 Myr
(*V. radiata*–*V. unguiculata*) to 56.0 Myr (Genistoid and Dalbergioid species against
all others).
