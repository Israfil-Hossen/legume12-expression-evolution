#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Step 5 - analysis for docs/step2_spec.md (12 legumes, control leaf expression)
##
## Implements spec sections 5-8 as written, plus deviation D1 (Cajanus out of the
## primary regression). Details the spec leaves open are marked "IMPLEMENTATION".
##
## Usage:  Rscript step05_analysis.R [project_dir]
## Output: <project_dir>/results_step2/
## ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
P <- if (length(args)) args[1] else "/home/liuchuanfa/israfil/legume2"
## Step 3 analysis C (docs/step3_spec.md): STEP5_HOG=1 repeats the primary test with OrthoFinder
## hierarchical orthogroups at the legume crown node (deviation D3) in place of orthogroups;
## output to results_step3C/.
HOG <- Sys.getenv("STEP5_HOG") == "1"
OUT <- file.path(P, if (HOG) "results_step3C" else "results_step2"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20260923)

NDRAW <- 50; NBOOT <- 2000; MATERIAL <- 0.050; NEGLIGIBLE <- 0.025
SPECIES <- c("Gmax","Pvul","Vrad","Vung","Ccaj","Mtru","Cari","Psat","Tpra","Ljap","Lang","Adur")
EXCL_PRIMARY <- "Ccaj"                                    # deviation D1
THIN_TO <- 5e6                                            # sensitivity 4

con <- file(file.path(OUT, "report.txt"), open = "wt")
say <- function(fmt, ...) { s <- if (length(list(...))) sprintf(fmt, ...) else fmt
  cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }
say("Step 5 analysis   %s   project: %s", format(Sys.time()), P)

## ============================================================== 1. samples + QC
runs <- fread(file.path(P, "metadata/step2_runs.tsv"))
S <- unique(runs[, .(species, group, biorep, layout)])
num_after <- function(f, section, key) {                 # first `key` after `section` (fastp JSON)
  x <- readLines(f, warn = FALSE); i <- grep(sprintf('"%s"', section), x)[1]
  j <- grep(sprintf('"%s"', key), x); j <- j[j > i][1]
  as.numeric(gsub("[^0-9.eE+-]", "", sub(".*:", "", x[j])))
}
meta_num <- function(f, key) {
  x <- paste(readLines(f, warn = FALSE), collapse = " ")
  m <- regmatches(x, regexpr(sprintf('"%s": *[0-9.eE+-]+', key), x))
  if (length(m)) as.numeric(sub(".*: *", "", m)) else NA_real_
}
rrna_ids <- function(sp) {                                # nuclear rRNA transcripts: NR_/XR_ and "ribosomal RNA|rRNA"
  f <- file.path(P, "tables", sprintf("%s_rRNA_tx.txt", sp))
  if (!file.exists(f)) {
    h <- system(sprintf("grep '^>' %s", shQuote(file.path(P, "refs", sp, "transcripts.fa"))), intern = TRUE)
    h <- h[grepl("^>(NR|XR)_", h) & grepl("ribosomal RNA|rRNA", h)]
    writeLines(sub("^>(\\S+).*", "\\1", h), f)
  }
  readLines(f)
}
Q <- list()                                               # transcript-level quant, kept for thinning
S[, done := file.exists(file.path(P, "quant", biorep, ".done"))]
S[, `:=`(pct_mapped = NA_real_, reads_after = NA_real_, rrna_frac = NA_real_, mapped_reads = NA_real_)]
for (i in which(S$done)) {
  b <- S$biorep[i]; sp <- S$species[i]
  q <- fread(file.path(P, "quant", b, "quant.sf"), select = c("Name","EffectiveLength","TPM","NumReads"))
  Q[[b]] <- q
  rr <- rrna_ids(sp)
  set(S, i, "pct_mapped", meta_num(file.path(P, "quant", b, "aux_info/meta_info.json"), "percent_mapped"))
  set(S, i, "reads_after", num_after(file.path(P, "qc", paste0(b, ".fastp.json")), "after_filtering", "total_reads") /
                           ifelse(S$layout[i] == "PAIRED", 2, 1))
  set(S, i, "rrna_frac", sum(q$NumReads[q$Name %in% rr]) / sum(q$NumReads))
  set(S, i, "mapped_reads", sum(q$NumReads))
}
S[, qc_reason := fifelse(!done, "not quantified",
                 fifelse(pct_mapped < 70, "mapping < 70%",
                 fifelse(reads_after < 5e6, "< 5M reads after fastp",
                 fifelse(rrna_frac > 0.5, "rRNA > 50%", "pass"))))]
S[, qc_pass := qc_reason == "pass"]
fwrite(S, file.path(OUT, "qc_samples.tsv"), sep = "\t")
say("\n== QC (spec section 5)")
say("samples quantified: %d of %d;  passing QC: %d", sum(S$done), nrow(S), sum(S$qc_pass))
for (r in S[, .(n = .N, pass = sum(qc_pass)), by = species][order(match(species, SPECIES))][
       , sprintf("  %-5s %2d of %2d pass", species, pass, n)]) say(r)
bad <- S[qc_reason != "pass"]
if (nrow(bad)) { say("excluded:"); for (r in bad[, sprintf("  %s  %s", biorep, qc_reason)]) say(r) }

## ============================================================== 2. orthogroups
ogf <- tail(sort(Sys.glob(file.path(P, "proteomes/OrthoFinder/Results_*", "Orthogroups/Orthogroups.tsv"))), 1)
if (!length(ogf)) stop("orthogroup file not found")
if (HOG) {
  ## deviation D3 (docs/step3_spec.md): OrthoFinder >= 3.1 writes the root HOGs (N0) as
  ## Orthogroups/Orthogroups.tsv, so N0 is the Step 2 input itself. Use the legume crown node:
  ## the node of OrthoFinder's rooted species tree whose descendants are exactly the 12 legumes.
  ## (Every node file keeps all species columns, so the node cannot be told from the header.)
  res <- dirname(dirname(ogf))
  tr <- ape::read.tree(file.path(res, "Species_Tree", "SpeciesTree_rooted_node_labels.txt"))
  pp <- ape::prop.part(tr); lab <- attr(pp, "labels")    # pp[[i]]: tips below node Ntip + i
  k <- which(vapply(pp, function(x) setequal(lab[x], SPECIES), logical(1)))
  if (length(k) != 1) stop(sprintf("expected one species-tree node with exactly the 12 legumes, found %d", length(k)))
  ogf <- file.path(res, "Phylogenetic_Hierarchical_Orthogroups", paste0(tr$node.label[k], ".tsv"))
  if (!file.exists(ogf)) stop("node file not found: ", ogf)
}
og <- fread(ogf, sep = "\t", header = TRUE, colClasses = "character")
n_out <- if (HOG) sum(rowSums(as.matrix(og[, setdiff(names(og), c("HOG", "OG", "Gene Tree Parent Clade", SPECIES)), with = FALSE]) != "", na.rm = TRUE) > 0) else 0
if (HOG) { og <- og[, c("HOG", intersect(SPECIES, names(og))), with = FALSE]; setnames(og, "HOG", "Orthogroup") }
long <- rbindlist(lapply(SPECIES, function(sp) {
  v <- og[[sp]]; v[is.na(v)] <- ""
  g <- strsplit(v, ", ", fixed = TRUE)
  data.table(Orthogroup = rep(og$Orthogroup, lengths(g)), species = sp,
             gene = sub(paste0("^", sp, "_"), "", unlist(g)))
}))
cn <- long[, .(copies = .N), by = .(Orthogroup, species)]
cnw <- dcast(cn, Orthogroup ~ species, value.var = "copies", fill = 0)
for (sp in setdiff(SPECIES, names(cnw))) cnw[, (sp) := 0L]
SC <- cnw[rowSums(as.matrix(cnw[, ..SPECIES]) == 1) == length(SPECIES), Orthogroup]
fwrite(cnw, file.path(OUT, "orthogroup_copy_numbers.tsv"), sep = "\t")
say("\n== Orthogroups: %s (%d rows in the file)", basename(dirname(dirname(ogf))), nrow(og))
if (HOG) say("hierarchical orthogroups (deviation D3): %s;  rows with genes outside the 12 legumes: %d", basename(ogf), n_out)
say("orthogroups with >= 1 copy in any legume: %d;  single-copy in all 12: %d", nrow(cnw), length(SC))

## ============================================================== 3. expression
tx2g <- lapply(setNames(SPECIES, SPECIES), function(sp)
  fread(file.path(P, "tables", sprintf("%s_tx2gene.tsv", sp)), header = FALSE, col.names = c("Name", "gene")))
og_tpm <- function(b, sp, thin) {                         # orthogroup TPM for one sample
  q <- Q[[b]]
  if (thin) {                                             # sensitivity 4: binomial thinning, TPM recomputed
    p <- min(1, THIN_TO / sum(q$NumReads))
    n <- rbinom(nrow(q), round(q$NumReads), p); r <- n / pmax(q$EffectiveLength, 1)
    q <- data.table(Name = q$Name, TPM = r / sum(r) * 1e6)
  }
  g <- merge(q[, .(Name, TPM)], tx2g[[sp]], by = "Name")[, .(TPM = sum(TPM)), by = gene]
  merge(long[species == sp], g, by = "gene", all.x = TRUE)[, .(tpm = sum(TPM, na.rm = TRUE)), by = Orthogroup]
}
## per species: log2 matrix, orthogroups with >= 1 copy in that species x QC-passed samples
build_expr <- function(thin = FALSE, per_copy = FALSE) {
  raw <- list()
  for (sp in SPECIES) {
    bs <- S[species == sp & qc_pass == TRUE, biorep]; if (!length(bs)) next
    ogs <- cn[species == sp, Orthogroup]
    M <- matrix(0, length(ogs), length(bs), dimnames = list(ogs, bs))
    for (b in bs) { x <- og_tpm(b, sp, thin); M[x$Orthogroup, b] <- x$tpm }
    raw[[sp]] <- M
  }
  ## IMPLEMENTATION: each sample scaled to the grand median of the per-sample single-copy medians
  sf <- unlist(unname(lapply(raw, function(M) apply(M[intersect(SC, rownames(M)), , drop = FALSE], 2, median))))
  C <- median(sf)
  lapply(setNames(names(raw), names(raw)), function(sp) {
    M <- raw[[sp]]; M <- sweep(M, 2, sf[colnames(M)], "/") * C
    if (per_copy) M <- M / cn[species == sp][match(rownames(M), Orthogroup), copies]
    log2(M + 1)
  })
}

## ============================================================== 4. draws
## IMPLEMENTATION: stratified split. Slots go one at a time to a random group that can still
## supply two more samples, so both halves always have the same group composition.
split_draw <- function(groups, k) {
  gs <- split(names(groups), groups); cap <- floor(lengths(gs) / 2)
  if (sum(cap) < k) return(NULL)
  al <- setNames(integer(length(gs)), names(gs))
  for (i in seq_len(k)) { av <- names(cap)[cap - al > 0]; pk <- av[sample.int(length(av), 1)]; al[pk] <- al[pk] + 1L }
  h1 <- h2 <- character(0)
  for (g in names(al)[al > 0]) { s <- gs[[g]][sample.int(length(gs[[g]]), 2 * al[g])]
    h1 <- c(h1, s[seq_len(al[g])]); h2 <- c(h2, s[al[g] + seq_len(al[g])]) }
  list(h1, h2)
}
within_div <- function(X, groups, k) {
  acc <- 0
  for (d in seq_len(NDRAW)) { h <- split_draw(groups, k); if (is.null(h)) return(NULL)
    acc <- acc + abs(rowMeans(X[, h[[1]], drop = FALSE]) - rowMeans(X[, h[[2]], drop = FALSE])) / sqrt(2) }
  acc / NDRAW
}
between_div <- function(XA, XB, k) {
  acc <- 0
  for (d in seq_len(NDRAW)) {
    a <- colnames(XA)[sample.int(ncol(XA), k)]; b <- colnames(XB)[sample.int(ncol(XB), k)]
    acc <- acc + abs(rowMeans(XA[, a, drop = FALSE]) - rowMeans(XB[, b, drop = FALSE])) / sqrt(2)
  }
  acc / NDRAW
}
srho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))

## ============================================================== 5. one full analysis
run_variant <- function(name, E, k = 2, species_ok = names(E), per_copy = FALSE, keep_div = FALSE) {
  sp_use <- intersect(SPECIES, species_ok)
  W <- rbindlist(lapply(sp_use, function(sp) {
    X <- E[[sp]]; grp <- setNames(S[match(colnames(X), biorep), group], colnames(X))
    if (ncol(X) < 2 * k) return(NULL)
    d <- within_div(X, grp, k); if (is.null(d)) return(NULL)
    c1 <- cn[species == sp][match(names(d), Orthogroup), copies]
    if (keep_div) fwrite(data.table(Orthogroup = names(d), copies = c1, div = d),
                         file.path(OUT, sprintf("div_within_%s_%s.tsv.gz", name, sp)), sep = "\t")
    data.table(species = sp, n_samples = ncol(X), n_og = length(d), rho = srho(c1, d))
  }))
  sp_b <- sp_use[sapply(sp_use, function(s) ncol(E[[s]]) >= k)]
  pairs <- t(combn(sp_b, 2))
  B <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    a <- pairs[i, 1]; b <- pairs[i, 2]
    m <- merge(cn[species == a, .(Orthogroup, ca = copies)], cn[species == b, .(Orthogroup, cb = copies)], by = "Orthogroup")
    if (per_copy) { pred <- sqrt(m$ca * m$cb)             # IMPLEMENTATION: geometric mean of the two copy numbers
    } else { m <- m[ca == cb]; pred <- m$ca }             # primary: equal copy number; shared number is the predictor
    d <- between_div(E[[a]][m$Orthogroup, , drop = FALSE], E[[b]][m$Orthogroup, , drop = FALSE], k)
    if (keep_div) fwrite(data.table(Orthogroup = m$Orthogroup, copies = pred, div = d),
                         file.path(OUT, sprintf("div_between_%s_%s_%s.tsv.gz", name, a, b)), sep = "\t")
    data.table(sp_a = a, sp_b = b, n_og = nrow(m), rho = srho(pred, d))
  }))
  fwrite(W, file.path(OUT, sprintf("rho_within_%s.tsv", name)), sep = "\t")
  fwrite(B, file.path(OUT, sprintf("rho_between_%s.tsv", name)), sep = "\t")
  list(W = W, B = B)
}

## ============================================================== 6. tests
TT <- fread(file.path(P, "metadata/step2_divergence_times.tsv"))
TT <- rbind(TT[, .(sp_a, sp_b, time = median_mya)], TT[, .(sp_a = sp_b, sp_b = sp_a, time = median_mya)])
pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b))
boot_species <- function(sp) sample(sp, length(sp), replace = TRUE)
pairs_of <- function(s) { ij <- t(combn(length(s), 2)); ij <- ij[s[ij[, 1]] != s[ij[, 2]], , drop = FALSE]
  pair_key(s[ij[, 1]], s[ij[, 2]]) }

primary_test <- function(B, excl, times = TT) {
  d <- merge(B, times, by = c("sp_a", "sp_b"))[!(sp_a %in% excl | sp_b %in% excl) & is.finite(rho)]
  d[, key := pair_key(sp_a, sp_b)]; setkey(d, key)
  rng <- diff(range(d$time))                              # fixed: observed range of the full data
  fit <- function(dd) if (nrow(dd) < 3 || var(dd$time) == 0) NA_real_ else coef(lm(rho ~ time, dd))[[2]] * rng
  est <- fit(d); sp <- unique(c(d$sp_a, d$sp_b))
  bt <- replicate(NBOOT, fit(d[J(pairs_of(boot_species(sp))), nomatch = 0L]))
  ci <- quantile(bt, c(.025, .975), na.rm = TRUE)
  v <- if (est >= MATERIAL && ci[1] > 0) "supports evolutionary account" else
       if (est <= NEGLIGIBLE) "supports technical account" else "partial"
  flag <- if (est <= -MATERIAL) "NOTE: material DECREASE with time - matches neither prediction" else ""
  data.table(test = "primary: rho_between ~ divergence time", n = nrow(d), estimate = est,
             lo = ci[[1]], hi = ci[[2]], verdict = v, flag = flag)
}
secondary_test <- function(W, B) {
  B <- copy(B)[is.finite(rho)][, key := pair_key(sp_a, sp_b)]; setkey(B, key)
  W <- W[is.finite(rho)]; sp <- unique(c(B$sp_a, B$sp_b))
  stat <- function(s) mean(B[J(pairs_of(s)), rho, nomatch = 0L]) - mean(W[match(s, species), rho], na.rm = TRUE)
  est <- mean(B$rho) - mean(W$rho)
  bt <- replicate(NBOOT, stat(boot_species(sp)))
  ci <- quantile(bt, c(.025, .975), na.rm = TRUE)
  v <- if (est >= MATERIAL && ci[1] > 0) "supports evolutionary account" else
       if (abs(est) <= NEGLIGIBLE) "supports technical account" else "partial"
  data.table(test = "secondary: mean rho_between - mean rho_within", n = nrow(B), estimate = est,
             lo = ci[[1]], hi = ci[[2]], verdict = v, flag = "")
}

## ============================================================== 7. run
variants <- list(
  list(name = "primary",           thin = FALSE, per_copy = FALSE, k = 2, drop = character(0)),
  list(name = "S1_per_copy",       thin = FALSE, per_copy = TRUE,  k = 2, drop = character(0)),
  list(name = "S2_three_per_side", thin = FALSE, per_copy = FALSE, k = 3, drop = character(0)),
  list(name = "S3_no_Pisum",       thin = FALSE, per_copy = FALSE, k = 2, drop = "Psat"),
  list(name = "S4_thinned_5M",     thin = TRUE,  per_copy = FALSE, k = 2, drop = character(0)))
if (HOG) variants <- variants[1]                            # analysis C: primary only
E_cache <- list(); tests <- list()
for (v in variants) {
  key <- paste(v$thin, v$per_copy)
  if (is.null(E_cache[[key]])) E_cache[[key]] <- build_expr(thin = v$thin, per_copy = v$per_copy)
  E <- E_cache[[key]]
  ok <- setdiff(names(E), v$drop)
  if (v$k == 3) ok <- ok[sapply(ok, function(s) ncol(E[[s]]) >= 6)]   # S2: species with >= 6 replicates
  say("\n== variant %s: %d species, %d samples per side", v$name, length(ok), v$k)
  R <- run_variant(v$name, E, k = v$k, species_ok = ok, per_copy = v$per_copy, keep_div = v$name == "primary")
  for (r in R$W[, sprintf("  within  %-5s n=%2d  rho = %+.4f  (%d orthogroups)", species, n_samples, rho, n_og)]) say(r)
  say("  between: %d pairs, rho %+.4f to %+.4f (median %+.4f)", nrow(R$B), min(R$B$rho, na.rm = TRUE),
      max(R$B$rho, na.rm = TRUE), median(R$B$rho, na.rm = TRUE))
  tt <- rbind(primary_test(R$B, EXCL_PRIMARY), secondary_test(R$W, R$B))[, variant := v$name]
  for (r in tt[, sprintf("  %-48s est %+.4f  [%+.4f, %+.4f]  -> %s %s", test, estimate, lo, hi, verdict, flag)]) say(r)
  tests[[v$name]] <- tt
  if (v$name == "primary") B_primary <- R$B
}
treef <- file.path(P, "tree", "concat.fa.treefile")
if (!HOG && file.exists(treef) && requireNamespace("ape", quietly = TRUE)) {
  D <- ape::cophenetic.phylo(ape::read.tree(treef))
  PD <- as.data.table(as.table(D))[, .(sp_a = as.character(V1), sp_b = as.character(V2), time = N)][sp_a != sp_b]
  t5 <- primary_test(B_primary, character(0), PD)[, `:=`(test = "S5: rho_between ~ patristic distance", variant = "S5_patristic")]
  say("\n== variant S5_patristic: primary between-species values vs patristic distance (Cajanus included, D1)")
  for (r in t5[, sprintf("  %-48s est %+.4f  [%+.4f, %+.4f]  -> %s %s", test, estimate, lo, hi, verdict, flag)]) say(r)
  tests[["S5_patristic"]] <- t5
} else if (HOG) say("%s", "Sensitivity 5 is not part of analysis C.") else say("Sensitivity 5 skipped: %s not found (run the tree job first).", treef)
fwrite(rbindlist(tests), file.path(OUT, "tests.tsv"), sep = "\t")
say("\nwrote %s", OUT)
close(con)
