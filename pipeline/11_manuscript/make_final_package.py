#!/usr/bin/env python3
"""Build FINAL_SUBMISSION/: (1) journal files, (2) Zenodo upload files with metadata and checksums, (3) a clean GitHub
repository folder with only this paper's code, specifications, metadata and results. Copies committed files only.
Usage (repository root): python pipeline/11_manuscript/make_final_package.py
"""
import hashlib, json, os, shutil, subprocess, zipfile

OUT = "FINAL_SUBMISSION"
J, Z, G = f"{OUT}/1_Journal_MPE", f"{OUT}/2_Zenodo_upload", f"{OUT}/3_GitHub_repository"
os.makedirs(OUT, exist_ok=True)
for x in os.listdir(OUT):                      # empty the folder but keep it (it may be open in File Explorer)
    p = os.path.join(OUT, x)
    shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
for d in (J, Z, G):
    os.makedirs(d)

TITLE = open("manuscript_mpe/manuscript.md", encoding="utf-8").readline().strip().strip("*")
AUTHORS = [  # name, ORCID, affiliation (as in the manuscript)
    ("Hossen, Md. Isrfil", "0009-0001-2189-1039", "State Key Laboratory of Soil and Sustainable Agriculture, Institute of Soil Science, Chinese Academy of Sciences, Nanjing, China"),
    ("Albaghdady, Dosouky Mohammad Dosouky", "0000-0001-9963-4465", "Department of Agricultural Zoology and Nematology, Faculty of Agriculture, Al-Azhar University, Cairo, Egypt"),
    ("Kabir, Md. Golam", "0000-0003-3387-9958", "State Key Laboratory of Soil and Sustainable Agriculture, Institute of Soil Science, Chinese Academy of Sciences, Nanjing, China"),
    ("Rahman, Md. Mosaddekur", "0009-0002-1274-1837", "Department of Agriculture, Faculty of Environmental Science, King Abdulaziz University, Jeddah, Saudi Arabia"),
    ("Karim, Md. Parvez", "0009-0007-3668-6373", "Department of Genetics and Plant Breeding, Bangladesh Agricultural University, Mymensingh, Bangladesh"),
    ("Helal, Md. Mostofa Uddin", "0000-0002-7033-3578", "Plant Molecular Biotechnology Laboratory, Department of Biochemistry and Molecular Biology, University of Rajshahi, Rajshahi, Bangladesh"),
]
DESCRIPTION = ("Code, specifications, metadata, quantification outputs and results for the study “" + TITLE + "”. "
    "142 public leaf RNA-seq libraries from 12 papilionoid legumes (15 studies) were processed uniformly (fastp, salmon), orthology was inferred with "
    "OrthoFinder and a species tree with IQ-TREE. The archive contains the analysis code for every step, the dated analysis specifications (Steps 2–5), "
    "the run tables and study metadata, salmon quantifications (quant.sf) for all 142 libraries, OrthoFinder results, the species-tree alignments and support "
    "analyses, paralog Ks tables, per-orthogroup divergence tables, all test results, figure source data and the numbers ledger from which every value in the "
    "article is generated. The between-species analyses were pre-registered at https://doi.org/10.17605/OSF.IO/S5QXK.")
KEYWORDS = ["gene duplication", "comparative transcriptomics", "expression evolution", "whole-genome duplication", "ohnologs", "Fabaceae",
            "batch effects", "pre-registration", "RNA-seq"]

# ---------------------------------------------------------------- 1. journal files
for f in os.listdir("manuscript_mpe/submission_MPE"):
    src = os.path.join("manuscript_mpe/submission_MPE", f)
    (shutil.copytree if os.path.isdir(src) else shutil.copy)(src, os.path.join(J, f))

# ---------------------------------------------------------------- 3. clean GitHub repository (only this paper)
CODE = ["pipeline/server", "pipeline/08_step2", "pipeline/09_step3", "pipeline/10_step4", "pipeline/11_manuscript", "pipeline/12_step5"]
DOCS = ["docs/step2_spec.md", "docs/step3_spec.md", "docs/step4_spec.md", "docs/step5_spec.md", "docs/preregistration.md"]
META = ["metadata/step2_runs.tsv", "metadata/step4_runs.tsv", "metadata/step2_divergence_times.tsv", "metadata/software_versions_server.tsv",
        "metadata/software_versions_local.tsv", "metadata/study_metadata"]
RESULTS = ["results/results_step2", "results/results_step3", "results/results_step3C", "results/results_step4", "results/results_step5",
           "results/exploratory_step2", "results/exploratory_step5", "results/exploratory_review", "results/tree_support",
           "results/manuscript_numbers", "results/server_logs"]
MS = ["manuscript_mpe/figures/source_data"]   # manuscript text and rendered figures stay private until publication
ign = shutil.ignore_patterns("__pycache__", "*.pyc", "*.tar.gz", "final_package_docs")
for p in CODE + DOCS + META + RESULTS + MS:
    dst = os.path.join(G, p)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    (shutil.copytree(p, dst, ignore=ign) if os.path.isdir(p) else shutil.copy(p, dst))
big = [os.path.join(r, f) for r, _, fs in os.walk(G) for f in fs if os.path.getsize(os.path.join(r, f)) > 95e6]
assert not big, f"files too large for GitHub: {big}"

prov = subprocess.run(["git", "log", "--date=iso", "--format=%h  %ad  %s", "--", *DOCS, "pipeline/08_step2", "pipeline/09_step3",
                       "pipeline/10_step4", "pipeline/12_step5/step5_analysis.py"], capture_output=True, text=True, encoding="utf-8").stdout
open(f"{G}/PROVENANCE_git_log.txt", "w", encoding="utf-8").write(
    "Dated commits of the analysis specifications and analysis code (from the authors' working repository; newest first).\n"
    "Steps 2 and 3 were also registered on OSF (https://doi.org/10.17605/OSF.IO/S5QXK); Steps 4 and 5 are documented by these commits.\n\n" + prov)

readme = f"""# {TITLE}

Code, specifications, metadata and results accompanying the manuscript above.

- Pre-registration (Steps 2 and 3): https://doi.org/10.17605/OSF.IO/S5QXK
- Archive with the large data files (salmon quantifications, OrthoFinder results, alignments): https://doi.org/10.5281/zenodo.23290174
- Every number in the article is generated by `pipeline/11_manuscript/make_ledger.R` into `results/manuscript_numbers/ledger.tsv`

## Contents

| Path | What it is |
|---|---|
| `docs/step2_spec.md` … `step5_spec.md` | Dated analysis specifications with decision rules (see `PROVENANCE_git_log.txt` for commit dates) |
| `pipeline/server/` | Cluster scripts: genomes, salmon indices, download, quantification, OrthoFinder, species tree |
| `pipeline/08_step2/` | Step 2: copy-number gradient within and between species, divergence-time test (R) |
| `pipeline/09_step3/` | Step 3: paralog Ks and noise floor, Brownian-motion model, hierarchical orthogroups; validation simulations |
| `pipeline/10_step4/` | Step 4: between-study comparison (R) |
| `pipeline/12_step5/` | Step 5: study share, expression level and read ambiguity, ohnolog asymmetry, relative dosage (Python); exploratory checks |
| `pipeline/11_manuscript/` | Numbers ledger, figures, supplementary tables, graphical abstract, manuscript rendering, tree-support analyses |
| `metadata/` | Run tables (ENA accessions, MD5), study characteristics, divergence times, software versions |
| `results/` | All test results, per-orthogroup divergence tables, exploratory analyses, tree support, server logs |
| `manuscript_mpe/figures/source_data/` | Source data for every figure panel |

## Reproducing the analyses

1. Cluster (Linux, PBS): `pipeline/server/step01_genomes.sh` → `step02_submit_index.sh` → `step03_download.sh` → `step04_submit_quant.sh`; OrthoFinder (`orthofinder.pbs`); species tree (`pipeline/08_step2/step05b_tree.pbs`).
2. Step 2: `Rscript pipeline/08_step2/step05_analysis.R` (needs the quantifications; on the cluster via `pipeline/server/step05.pbs`).
3. Step 3: `pipeline/server/step3AB.pbs`, `step3C.pbs`. Step 4: `pipeline/10_step4/step08_step4.pbs`.
4. Step 5 (local): unpack the Zenodo quantification archive to a folder `DATA`, then `python pipeline/12_step5/step5_analysis.py DATA`.
5. Numbers, figures and supplementary tables: `Rscript pipeline/11_manuscript/make_ledger.R`; `python pipeline/11_manuscript/make_figures.py`;
   `python pipeline/11_manuscript/make_supp_tables.py`.

Software versions are listed in `metadata/software_versions_*.tsv`.

## Licence

Code: MIT (see `LICENSE`). Data, specifications and documentation: CC BY 4.0.

## Citation

Hossen MI, Albaghdady DMD, Kabir MG, Rahman MM, Karim MP, Helal MMU. {TITLE}. Manuscript; the published version will be linked here.
"""
open(f"{G}/README.md", "w", encoding="utf-8").write(readme)
open(f"{G}/LICENSE", "w", encoding="utf-8").write(
    "MIT License\n\nCopyright (c) 2026 Md. Isrfil Hossen, Dosouky Mohammad Dosouky Albaghdady, Md. Golam Kabir, Md. Mosaddekur Rahman, "
    "Md. Parvez Karim, Md. Mostofa Uddin Helal\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of this software and "
    "associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, "
    "modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, "
    "subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all copies or substantial portions "
    "of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF "
    "MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, "
    "DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR "
    "OTHER DEALINGS IN THE SOFTWARE.\n")
cff = ["cff-version: 1.2.0", "message: \"If you use this material, please cite the article.\"", f"title: \"{TITLE}\"", "authors:"]
for name, orcid, aff in AUTHORS:
    fam, giv = name.split(", ")
    cff += [f"  - family-names: \"{fam}\"", f"    given-names: \"{giv}\"", f"    orcid: \"https://orcid.org/{orcid}\"", f"    affiliation: \"{aff}\""]
cff += ["license: MIT", "keywords:"] + [f"  - {k}" for k in KEYWORDS]
open(f"{G}/CITATION.cff", "w", encoding="utf-8").write("\n".join(cff) + "\n")
zen = {"title": TITLE + ": code, data and results", "upload_type": "dataset", "description": DESCRIPTION, "access_right": "open", "license": "cc-by-4.0",
       "creators": [{"name": n, "orcid": o, "affiliation": a} for n, o, a in AUTHORS], "keywords": KEYWORDS,
       "related_identifiers": [{"identifier": "10.17605/OSF.IO/S5QXK", "relation": "isDocumentedBy", "scheme": "doi"}]}
json.dump(zen, open(f"{G}/.zenodo.json", "w", encoding="utf-8"), indent=2, ensure_ascii=False)
open(f"{G}/.gitignore", "w").write("__pycache__/\n*.pyc\n*.tar.gz\n.Rhistory\n")

# ---------------------------------------------------------------- 2. Zenodo upload
with zipfile.ZipFile(f"{Z}/01_code_specifications_metadata_results.zip", "w", zipfile.ZIP_DEFLATED) as zf:
    for r, _, fs in os.walk(G):
        for f in fs:
            full = os.path.join(r, f)
            zf.write(full, os.path.join("legume12_expression_evolution", os.path.relpath(full, G)))
shutil.copy("results/legume2_deposit.tar.gz", f"{Z}/02_quantification_orthology_tree_paralogs.tar.gz")
shutil.copy("results/ambig_dup.tar.gz", f"{Z}/03_salmon_read_ambiguity.tar.gz")
json.dump(zen, open(f"{Z}/zenodo_metadata.json", "w", encoding="utf-8"), indent=2, ensure_ascii=False)
files = sorted(f for f in os.listdir(Z) if f.endswith((".zip", ".tar.gz")))
sums = []
for f in files:
    h = hashlib.sha256()
    with open(os.path.join(Z, f), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    sums.append(f"{h.hexdigest()}  {f}")
open(f"{Z}/SHA256SUMS.txt", "w").write("\n".join(sums) + "\n")
sizes = {f: os.path.getsize(os.path.join(Z, f)) for f in files}
print("journal files:", len(os.listdir(J)), "| GitHub files:", sum(len(fs) for _, _, fs in os.walk(G)), "| Zenodo:", {k: f"{v/1e6:.0f} MB" for k, v in sizes.items()})

# ---------------------------------------------------------------- instructions (kept in the repository so a rebuild keeps them)
shutil.copy("pipeline/11_manuscript/final_package_docs/00_README_FIRST.md", f"{OUT}/00_README_FIRST.md")
shutil.copy("pipeline/11_manuscript/final_package_docs/ZENODO_FORM_FIELDS.md", f"{Z}/ZENODO_FORM_FIELDS.md")
