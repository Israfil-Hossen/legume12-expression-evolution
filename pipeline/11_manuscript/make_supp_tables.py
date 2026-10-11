#!/usr/bin/env python3
"""Supplementary Tables S1-S12 as one Excel workbook, assembled only from committed result files
(no statistic is computed here except joins, reshaping and counting rows).
Usage (repository root): python pipeline/11_manuscript/make_supp_tables.py
"""
import glob, os, re
import pandas as pd

RES, OUT = "results", "manuscript_mpe/Supplementary_Tables.xlsx"
SPECIES = ["Gmax", "Pvul", "Vrad", "Vung", "Ccaj", "Mtru", "Cari", "Psat", "Tpra", "Ljap", "Lang", "Adur"]
LATIN = dict(Gmax="Glycine max", Pvul="Phaseolus vulgaris", Vrad="Vigna radiata", Vung="Vigna unguiculata",
             Ccaj="Cajanus cajan", Mtru="Medicago truncatula", Cari="Cicer arietinum", Psat="Pisum sativum",
             Tpra="Trifolium pratense", Ljap="Lotus japonicus", Lang="Lupinus angustifolius", Adur="Arachis duranensis")
rd = lambda p, **k: pd.read_csv(p, sep="\t", **k)
T = {}

# Table S2: samples ---
runs = pd.concat([rd("metadata/step2_runs.tsv").assign(set="between-species (Step 2)"),
                  rd("metadata/step4_runs.tsv").assign(set="second study (Step 4)")])
agg = runs.groupby("biorep").agg(set=("set", "first"), species=("species", "first"), study=("study", "first"),
                                  group=("group", "first"), runs=("run", lambda x: ";".join(x)), layout=("layout", "first"),
                                  read_count_ENA=("read_count", "sum")).reset_index()
qc = pd.concat([rd(f"{RES}/results_step2/qc_samples.tsv"), rd(f"{RES}/results_step4/qc_samples_step4.tsv")])
s1 = agg.merge(qc[["biorep", "pct_mapped", "reads_after", "rrna_frac", "mapped_reads", "qc_reason", "qc_pass"]], on="biorep", how="left")
s1.insert(2, "species_name", s1.species.map(LATIN))
rm = pd.read_csv("metadata/study_metadata/run_level_metadata.tsv", sep="\t")
rmb = rm.groupby("biorep").agg(genotype_recorded=("genotype", "first"), age_recorded=("age", "first"), instrument=("instrument_model", "first"),
                               read_length=("read_len", lambda x: "/".join(str(int(v)) for v in sorted(set(x)))), library_selection=("library_selection", "first"),
                               biosample=("sample_accession", "first")).reset_index()
s1 = s1.merge(rmb, on="biorep", how="left")
T["S2 Samples"] = (s1.sort_values(["set", "species", "biorep"]),
                   "All biological replicates. read_count_ENA, sequenced fragments (read pairs counted once); reads_after, after fastp (pairs counted once); "
                   "pct_mapped, salmon mapping rate; rrna_frac, fraction of reads assigned to nuclear rRNA transcripts; qc_pass, passed all thresholds; "
                   "genotype_recorded, age_recorded, as given in the BioSample record (retrieved 9 October 2026); read_length, bases per read.")
T["S2b Studies"] = (pd.read_csv("metadata/study_metadata/study_characteristics.tsv", sep="\t"),
                    "The 15 studies: genotypes used, plant age or stage, platform, layout, read length and library selection, as recorded in BioSample/ENA.")

# Table S1: divergence times ---
tt = rd("metadata/step2_divergence_times.tsv")
sp2 = rd(f"{RES}/results_step5/2B_species_pairs.tsv")[["sp_a", "sp_b", "patristic"]]
s2 = tt.merge(pd.concat([sp2, sp2.rename(columns={"sp_a": "sp_b", "sp_b": "sp_a"})]), on=["sp_a", "sp_b"], how="left")
T["S1b Tree support"] = (pd.read_csv(f"{RES}/tree_support/support_summary.tsv", sep="	"),
                         "Species tree: SH-aLRT and UFBoot (1,000 replicates each) under LG+G4 and under a partitioned model with ModelFinder (one partition "
                         "per orthogroup), gene concordance (gCF, 725 gene trees) and likelihood site concordance (sCF) per branch, with the discordance "
                         "factors for the two alternative resolutions (DF1, DF2) and polyphyly (gDFP); IQ-TREE 3.1.3.")
T["S1 Divergence times"] = (s2, "TimeTree pairwise ages retrieved and frozen on 23 September 2026 (median, interval, number of studies, synthesis-tree "
                                "adjusted age) and patristic distance on the IQ-TREE species tree (substitutions/site).")

# Table S6: within-study gradients ---
w = rd(f"{RES}/results_step2/rho_within_primary.tsv").rename(columns={"rho": "rho_raw"})
b1 = rd(f"{RES}/results_step5/1B_partial.tsv"); b1 = b1[b1.level == "within"].rename(columns={"label": "species"}).drop(columns=["level", "rho_raw"])
c1 = rd(f"{RES}/results_step5/1C_ambiguity.tsv")
dep = rd(f"{RES}/exploratory_step2/within_rho_vs_depth.tsv")[["species", "median_reads_M"]]
s3 = w.merge(b1, on="species").merge(c1, on="species", suffixes=("", "_amb")).merge(dep, on="species")
T["S6 Within-study gradient"] = (s3, "Copy-number gradient between replicate halves of one study (Spearman rho), raw and as partial correlations controlling for "
                                     "mean expression, per-copy expression, and expression plus read ambiguity; median rho within expression deciles; rho between "
                                     "copy number and read ambiguity; median reads (millions).")

# Table S7: between-species gradients ---
keys = ["sp_a", "sp_b"]
s4 = rd(f"{RES}/results_step2/rho_between_primary.tsv").rename(columns={"n_og": "n_og_primary", "rho": "rho_primary"})
for v in ("S1_per_copy", "S2_three_per_side", "S3_no_Pisum", "S4_thinned_5M"):
    x = rd(f"{RES}/results_step2/rho_between_{v}.tsv").rename(columns={"n_og": f"n_og_{v}", "rho": f"rho_{v}"})
    s4 = s4.merge(x, on=keys, how="left")
c = rd(f"{RES}/results_step3C/rho_between_primary.tsv").rename(columns={"n_og": "n_og_C_HOG", "rho": "rho_C_HOG"})
s4 = s4.merge(c, on=keys, how="left").merge(tt[["sp_a", "sp_b", "median_mya"]], on=keys, how="left")
pb = rd(f"{RES}/results_step5/1B_partial.tsv"); pb = pb[pb.level == "between"].copy()
pb[["sp_a", "sp_b"]] = pb.label.str.split("-", expand=True)
s4 = s4.merge(pb[["sp_a", "sp_b", "rho_partial_expr", "rho_partial_percopy"]], on=keys, how="left")
r4 = rd(f"{RES}/exploratory_review/R4_matched_within.tsv")[["sp_a", "sp_b", "rho_within_matched", "diff_matched"]]
r2 = rd(f"{RES}/exploratory_review/R2_absolute_gradient_pairs.tsv").drop(columns=["time"])
s4 = s4.merge(r4, on=keys, how="left").merge(r2.add_prefix("abs_").rename(columns={"abs_sp_a": "sp_a", "abs_sp_b": "sp_b"}), on=keys, how="left")
T["S7 Between-species gradient"] = (s4, "Between-species copy-number gradient for all 66 species pairs: primary analysis, sensitivity analyses S1-S4 and hierarchical "
                                        "orthogroups (C), with the number of orthogroups compared and the TimeTree age; rho_partial_expr, partial rho controlling for mean expression (1B); "
                                        "EXPLORATORY: rho_within_matched, mean within-study rho of the two species on the pair's orthogroups, and diff_matched, rho_primary minus it; "
                                        "abs_*, absolute-scale gradient (median divergence of multi-copy minus single-copy orthogroups, log2 units, and log2 ratio of the two medians).")

# Table S8: between-study comparison ---
p4 = rd(f"{RES}/results_step4/pairs_step4.tsv").rename(columns={"lo": "delta_lo", "hi": "delta_hi"})
a1 = rd(f"{RES}/results_step5/1A_pairs.tsv").drop(columns=["n_og"])
r5 = rd(f"{RES}/exploratory_review/R5_F2.tsv")[["focal", "other", "F2"]]
r3 = rd(f"{RES}/exploratory_review/R3_absolute_study_vs_species.tsv").drop(columns=["n_og"]).add_prefix("abs_").rename(columns={"abs_focal": "focal", "abs_other": "other"})
T["S8 Between-study"] = (p4.merge(a1, on=["focal", "other"]).merge(r5, on=["focal", "other"]).merge(r3, on=["focal", "other"]), "For each species with two studies (focal) and each other species: gradient between species and between "
                         "the two studies on identical orthogroups, Delta with 95% orthogroup-bootstrap interval, median expression distances within study, "
                         "between studies and between species, and the study share F ('laboratory share' in the specification); EXPLORATORY: F2, variance-based share; "
                         "abs_*, absolute-scale gradients between studies and between species on the same orthogroups.")

# Table S9: paralogs ---
rows = []
for f in sorted(glob.glob(f"{RES}/server_logs/paralogs/*.log")):
    t = open(f).read()
    m = re.search(r"(\w+): genes with longest protein (\d+); kept with matching CDS (\d+)", t)
    n = re.search(r"reciprocal-best-hit pairs (\d+); with Ks (\d+); collinear (\d+)", t)
    rows.append(dict(species=m.group(1), genes=int(m.group(2)), genes_with_cds=int(m.group(3)), rbh_pairs=int(n.group(1)),
                     pairs_with_ks=int(n.group(2)), collinear_pairs=int(n.group(3))))
s6 = pd.DataFrame(rows).merge(rd(f"{RES}/results_step3/A_per_species.tsv"), on="species", how="left")
T["S9 Paralogs"] = (s6, "Paralog pairs per genome (reciprocal best DIAMOND hits), with Ks and MCScanX collinearity; n_pairs, pairs used (Ks <= 3, expressed); "
                        "correlations with Ks of D (A1), the noise floor N (A2), D - N (registered, biased) and the mean product m1*m2 (A3 after deviation D2).")

# Table S10: ohnologs ---
pr = rd(f"{RES}/results_step5/2_pairs.tsv")
mr = rd(f"{RES}/results_step5/2B_mean_ratios.tsv").set_index("OG")
mr.columns = [f"mean_ratio_{c}" for c in mr.columns]
info = pr.groupby("OG").agg(HOG_A=("HOG_A", "first"), HOG_B=("HOG_B", "first"), n_species=("species", "size"),
                            frac_collinear=("collinear", "mean"), ks_median=("ks", "median"))
s7 = info.join(mr).join(rd(f"{RES}/results_step5/2B_direction.tsv").set_index("OG")[["n_A_higher", "p", "q"]]) \
         .join(rd(f"{RES}/results_step5/2C_classes.tsv").set_index("OG")[["cls", "mean_abs_r", "frac_same_direction", "category"]]).reset_index()
T["S10 Ohnologs"] = (s7, "Ohnolog orthogroups from the papilionoid whole-genome duplication: the two legume-crown hierarchical orthogroups (copy A, copy B), "
                        "species with a usable pair, collinear fraction, median Ks, mean log2(A/B) in each species, number of species with copy A higher and "
                        "binomial test (p, BH q), descriptive class and A. thaliana annotation category.")

# Table S11: dosage summary ---
t5 = rd(f"{RES}/results_step5/tests_step5.tsv")
s8a = t5[t5.analysis.str.startswith("3_")].assign(set="Glycine whole-genome duplicates")
s8b = rd(f"{RES}/results_step5/3_tandem_descriptive.tsv").assign(set="tandem/proximal lineage-specific duplicates (descriptive)")
s8c = rd(f"{RES}/exploratory_review/R1_dosage_by_expression.tsv").assign(set="EXPLORATORY: Glycine whole-genome duplicates by quartile of reference expression")
T["S11 Dosage summary"] = (pd.concat([s8a, s8b, s8c], ignore_index=True), "Relative dosage (log2): registered estimate for Glycine whole-genome duplicates "
                          "(both G. max studies) and, descriptively, for tandem/proximal duplicates present in one species only; EXPLORATORY: by quartile of reference "
                          "expression (scaled TPM of the single-copy species), with the value full dosage would give under log2(x + 1) and the fraction of it observed.")

# Table S12: per-orthogroup dosage ---
s9 = rd(f"{RES}/results_step5/3_dosage_Gmax_study1.tsv.gz").rename(columns={"D": "D_study1"}).merge(
     rd(f"{RES}/results_step5/3_dosage_Gmax_study2.tsv.gz").rename(columns={"D": "D_study2"}), on=["OG", "set"])
T["S12 Dosage per orthogroup"] = (s9, "Per-orthogroup relative family expression of G. max minus the mean of single-copy species (log2), for Glycine whole-genome "
                                     "duplicate orthogroups and the calibration set (single-copy in all 12 species), with each G. max study.")

# Table S4: all registered tests ---
t2 = rd(f"{RES}/results_step2/tests.tsv").assign(step="2")
t3 = rd(f"{RES}/results_step3/tests_step3AB.tsv").rename(columns={"analysis": "test"}).assign(step="3")
t3c = rd(f"{RES}/results_step3C/tests.tsv").assign(step="3C")
t4 = rd(f"{RES}/results_step4/tests_step4.tsv").rename(columns={"level": "test", "delta": "estimate"}).assign(step="4")
t5r = t5.rename(columns={"analysis": "test"}).assign(step="5")
s10 = pd.concat([t2, t3, t3c, t4, t5r], ignore_index=True)[["step", "test", "variant", "n", "estimate", "lo", "hi", "verdict", "flag"]]
T["S4 All tests"] = (s10, "All pre-registered tests at full precision, with the verdict in the wording of each specification.")

# Table S5: software and parameters ---
sv = pd.concat([rd("metadata/software_versions_server.tsv"), rd("metadata/software_versions_local.tsv")], ignore_index=True)
led = rd(f"{RES}/manuscript_numbers/ledger.tsv", dtype=str)
extra = led[led.id.str.startswith("ver_") & led.id.str.endswith("dupenv") | led.id.isin(["ver_salmon", "ver_orthofinder", "ver_iqtree"])]
sv = pd.concat([sv, pd.DataFrame(dict(tool=extra.id.str.replace("ver_", "").str.replace("_dupenv", ""), version=extra.value,
                                      environment=extra.note.fillna(""), **{"how recorded": extra.source}))], ignore_index=True)
params = led[led.id.str.startswith("des_")][["id", "value", "source"]].rename(columns={"id": "parameter"})
T["S5 Software"] = (sv, "Software versions and how each was recorded.")
T["S5b Parameters"] = (params, "Design constants fixed in the specifications.")

# Table S3: deviations and exploratory analyses ---
s12 = pd.DataFrame([
    ("D1", "deviation", "Step 2", "Cajanus cajan removed from divergence-time regressions because TimeTree places it outside the Millettioids; kept elsewhere.", "before analysis"),
    ("D2", "deviation", "Step 3", "Registered paralog statistic (Spearman rho of D - N with Ks) shown by simulation to be biased; replaced by the Pearson correlation of the mean product m1*m2 with Ks; registered version still reported.", "before Step 3 data"),
    ("D3", "deviation", "Step 3", "OrthoFinder 3 writes root-node HOGs as its orthogroups; analysis C used the legume crown node (N1).", "before analysis C produced results"),
    ("E1-E4", "exploratory", "Step 2", "Post hoc exploratory analyses of the Step 2 results (dilution, attenuation reference, depth, joint model).", "after Step 2 results"),
    ("Step 3 class-by-Ks", "exploratory (pre-declared)", "Step 3", "Collinear versus other paralogs by Ks bin, pooled over species (Fig. S7).", "declared in the Step 3 specification without decision rule"),
    ("2B robustness", "exploratory", "Step 5", "Four robustness checks of the ohnolog-asymmetry trend (Fig. S8).", "after the 2B result"),
    ("R1", "exploratory", "Step 5", "Relative dosage by quartile of reference expression, against the value full dosage gives under log2(x + 1) (Table S11).", "after all registered results; prompted by an internal pre-submission review"),
    ("R2-R3", "exploratory", "Steps 2 and 4", "Absolute-scale copy-number gradients within studies, between species and between studies, and their change with divergence time (Tables S7, S8).", "after all registered results; prompted by an internal pre-submission review"),
    ("R4", "exploratory", "Step 2", "Secondary test with the within-study gradient recomputed on each pair's orthogroup set (Table S7).", "after all registered results; prompted by an internal pre-submission review"),
    ("R5", "exploratory", "Step 5", "Variance-based study share F2 (Table S8).", "after all registered results; prompted by an internal pre-submission review"),
    ("R6", "exploratory", "Step 5", "Permutation (Mantel-type) test of the ohnolog-asymmetry trend, species labels of the divergence-time matrix permuted 10,000 times.", "after all registered results; prompted by an internal pre-submission review"),
], columns=["id", "type", "step", "description", "timing"])
T["S3 Deviations"] = (s12, "Deviations from the specifications and analyses labelled exploratory.")

T = dict(sorted(T.items(), key=lambda kv: (int(re.match(r"S(\d+)", kv[0]).group(1)), kv[0])))
with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    readme = pd.DataFrame([(k, v[1], len(v[0])) for k, v in T.items()], columns=["sheet", "description", "rows"])
    readme.to_excel(xw, sheet_name="README", index=False)
    for k, (df, _) in T.items():
        df.to_excel(xw, sheet_name=k[:31], index=False)
print(f"wrote {OUT}")
print(readme.to_string(index=False, max_colwidth=60))
