#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Step 4 - two studies of one species, the "between-study" rung (docs/step4_spec.md).
## Details the spec leaves open are marked IMPLEMENTATION.
##
## Inputs: Step 2 outputs in <project>/results_step2 (QC table, per-orthogroup between- and
##         within-species divergences), Step 2 and Step 4 quantification, metadata/step4_runs.tsv.
## Usage:  Rscript step08_step4.R [project_dir]
## Output: <project_dir>/results_step4/
## ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
P <- if (length(args)) args[1] else "/home/liuchuanfa/israfil/legume2"
R2 <- file.path(P, "results_step2")
OUT <- file.path(P, "results_step4"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20260924)

NDRAW <- 50; NBOOT <- 2000; MATERIAL <- 0.050; NEGLIGIBLE <- 0.025; K <- 2
SPECIES <- c("Gmax","Pvul","Vrad","Vung","Ccaj","Mtru","Cari","Psat","Tpra","Ljap","Lang","Adur")
FOCAL <- c("Gmax", "Pvul", "Mtru")

con <- file(file.path(OUT, "report_step4.txt"), open = "wt")
say <- function(fmt, ...) { s <- if (length(list(...))) sprintf(fmt, ...) else fmt
  cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }
srho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
say("Step 4 analysis   %s   project: %s", format(Sys.time()), P)

## ============================================================== 1. samples + QC (as Step 2)
S1 <- fread(file.path(R2, "qc_samples.tsv"))[qc_pass == TRUE]        # study 1 = the Step 2 samples
runs <- fread(file.path(P, "metadata/step4_runs.tsv"))
S2 <- unique(runs[, .(species, group, biorep, layout)])
num_after <- function(f, section, key) {
  x <- readLines(f, warn = FALSE); i <- grep(sprintf('"%s"', section), x)[1]
  j <- grep(sprintf('"%s"', key), x); j <- j[j > i][1]
  as.numeric(gsub("[^0-9.eE+-]", "", sub(".*:", "", x[j])))
}
meta_num <- function(f, key) {
  x <- paste(readLines(f, warn = FALSE), collapse = " ")
  m <- regmatches(x, regexpr(sprintf('"%s": *[0-9.eE+-]+', key), x))
  if (length(m)) as.numeric(sub(".*: *", "", m)) else NA_real_
}
rrna_ids <- function(sp) readLines(file.path(P, "tables", sprintf("%s_rRNA_tx.txt", sp)))   # written by Step 2
Q <- list()
S2[, done := file.exists(file.path(P, "quant", biorep, ".done"))]
S2[, `:=`(pct_mapped = NA_real_, reads_after = NA_real_, rrna_frac = NA_real_, mapped_reads = NA_real_)]
for (i in which(S2$done)) {
  b <- S2$biorep[i]; sp <- S2$species[i]
  q <- fread(file.path(P, "quant", b, "quant.sf"), select = c("Name", "TPM", "NumReads"))
  Q[[b]] <- q[, .(Name, TPM)]
  set(S2, i, "pct_mapped", meta_num(file.path(P, "quant", b, "aux_info/meta_info.json"), "percent_mapped"))
  set(S2, i, "reads_after", num_after(file.path(P, "qc", paste0(b, ".fastp.json")), "after_filtering", "total_reads") /
                            ifelse(S2$layout[i] == "PAIRED", 2, 1))
  set(S2, i, "rrna_frac", sum(q$NumReads[q$Name %in% rrna_ids(sp)]) / sum(q$NumReads))
  set(S2, i, "mapped_reads", sum(q$NumReads))
}
S2[, qc_reason := fifelse(!done, "not quantified",
                  fifelse(pct_mapped < 70, "mapping < 70%",
                  fifelse(reads_after < 5e6, "< 5M reads after fastp",
                  fifelse(rrna_frac > 0.5, "rRNA > 50%", "pass"))))]
S2[, qc_pass := qc_reason == "pass"]
fwrite(S2, file.path(OUT, "qc_samples_step4.tsv"), sep = "\t")
say("\n== QC of study 2 (Step 2 thresholds)")
say("samples quantified: %d of %d;  passing QC: %d", sum(S2$done), nrow(S2), sum(S2$qc_pass))
for (r in S2[, .(n = .N, pass = sum(qc_pass)), by = species][, sprintf("  %-5s %2d of %2d pass", species, pass, n)]) say(r)
bad <- S2[qc_reason != "pass"]
if (nrow(bad)) { say("excluded:"); for (r in bad[, sprintf("  %s  %s", biorep, qc_reason)]) say(r) }
if (!all(FOCAL %in% S2$species)) stop("step4_runs.tsv does not cover the three focal species")

## ============================================================== 2. orthogroups (the Step 2 run)
ogf <- tail(sort(Sys.glob(file.path(P, "proteomes/OrthoFinder/Results_*", "Orthogroups/Orthogroups.tsv"))), 1)
if (!length(ogf)) stop("orthogroup file not found")
og <- fread(ogf, sep = "\t", header = TRUE, colClasses = "character")
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
say("\n== Orthogroups: %s;  single-copy in all 12: %d", basename(dirname(dirname(ogf))), length(SC))

## ============================================================== 3. expression (Step 2 scaling)
tx2g <- lapply(setNames(SPECIES, SPECIES), function(sp)
  fread(file.path(P, "tables", sprintf("%s_tx2gene.tsv", sp)), header = FALSE, col.names = c("Name", "gene")))
og_tpm <- function(b, sp) {                               # orthogroup TPM for one sample, all orthogroups of sp
  q <- if (is.null(Q[[b]])) fread(file.path(P, "quant", b, "quant.sf"), select = c("Name", "TPM")) else Q[[b]]
  g <- merge(q, tx2g[[sp]], by = "Name")[, .(TPM = sum(TPM)), by = gene]
  x <- merge(long[species == sp], g, by = "gene", all.x = TRUE)[, .(tpm = sum(TPM, na.rm = TRUE)), by = Orthogroup]
  setNames(x$tpm, x$Orthogroup)
}
## Step 2 scaled every sample to C, the grand median of the per-sample single-copy medians of all
## QC-passed Step 2 samples. IMPLEMENTATION: C is recomputed here from the same samples, so study-1
## values equal the Step 2 values; study-2 samples are divided by their own single-copy median, times C.
raw1 <- list(); sf1 <- c()
for (i in seq_len(nrow(S1))) {
  b <- S1$biorep[i]; sp <- S1$species[i]; x <- og_tpm(b, sp)
  sf1[b] <- median(x[SC])
  if (sp %in% FOCAL) raw1[[b]] <- x
}
C <- median(sf1)
say("scaling constant C (median of %d Step 2 single-copy medians): %.4f", length(sf1), C)
expr <- function(sp, bioreps, cache) {                    # log2 matrix: orthogroups of sp x samples
  ogs <- cn[species == sp, Orthogroup]
  M <- sapply(bioreps, function(b) { x <- if (is.null(cache[[b]])) og_tpm(b, sp) else cache[[b]]; x[ogs] / median(x[SC]) * C })
  dimnames(M) <- list(ogs, bioreps); log2(M + 1)
}

## ============================================================== 4. draws (as Step 2)
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

## ============================================================== 5. per focal species
## IMPLEMENTATION: one set of 50 between-study draws per focal species gives the per-orthogroup
## between-study divergence; it is then restricted to the orthogroup set of each pair (A, B).
pair_file <- function(a, b) { s <- SPECIES[sort(match(c(a, b), SPECIES))]
  file.path(R2, sprintf("div_between_primary_%s_%s.tsv.gz", s[1], s[2])) }
within_file <- function(sp) file.path(R2, sprintf("div_within_primary_%s.tsv.gz", sp))
atten <- function(copies, noise, div) {                    # the manuscript's attenuation reference
  rn <- srho(copies, noise); list(rho_noise = rn, rho_exp = rn * sd(noise) / sd(div)) }
rho2 <- fread(file.path(R2, "rho_between_primary.tsv"))
PAIRS <- list(); DESC <- list(); BOOT <- list()
for (A in FOCAL) {
  b1 <- S1[species == A, biorep]; b2 <- S2[species == A & qc_pass == TRUE, biorep]
  if (length(b2) < 2 * K) stop(sprintf("%s: fewer than %d QC-passed study-2 samples", A, 2 * K))
  X1 <- expr(A, b1, raw1); X2 <- expr(A, b2, list())
  cA <- cn[species == A][match(rownames(X1), Orthogroup), copies]
  dS <- between_div(X1, X2, K)                             # between-study divergence, all orthogroups of A
  w2 <- within_div(X2, setNames(S2[match(b2, biorep), group], b2), K)
  w1f <- fread(within_file(A)); w1 <- w1f$div[match(rownames(X1), w1f$Orthogroup)]
  fwrite(data.table(Orthogroup = rownames(X1), copies = cA, div_study = dS, div_within_study1 = w1,
                    div_within_study2 = if (is.null(w2)) NA_real_ else w2),
         file.path(OUT, sprintf("div_study_%s.tsv.gz", A)), sep = "\t")
  at <- if (is.null(w2)) list(rho_noise = NA_real_, rho_exp = NA_real_) else atten(cA, (w1 + w2) / 2, dS)
  DESC[[A]] <- data.table(species = A, n_study1 = length(b1), n_study2 = length(b2), n_og = length(dS),
    rho_study_all = srho(cA, dS), rho_within_study1 = srho(cA, w1),
    rho_within_study2 = if (is.null(w2)) NA_real_ else srho(cA, w2),
    study_rho_expected = at$rho_exp, study_obs_over_exp = srho(cA, dS) / at$rho_exp)
  sp_oe <- c()
  for (B in setdiff(SPECIES, A)) {
    f <- pair_file(A, B); if (!file.exists(f)) { say("  %s-%s: no Step 2 between-species file - skipped", A, B); next }
    d <- fread(f)                                          # Orthogroup, copies, div: the Step 2 set and values
    ds <- dS[match(d$Orthogroup, rownames(X1))]
    if (anyNA(ds)) stop(sprintf("%s-%s: orthogroups of the Step 2 set missing from the %s matrix", A, B, A))
    r_sp <- srho(d$copies, d$div); r_st <- srho(d$copies, ds)
    reg <- rho2[(sp_a == A & sp_b == B) | (sp_a == B & sp_b == A), rho]
    if (length(reg) != 1 || abs(reg - r_sp) > 1e-6) stop(sprintf("%s-%s: rho from the divergence file differs from the Step 2 table", A, B))
    n <- nrow(d)
    bt <- vapply(seq_len(NBOOT), function(i) { j <- sample.int(n, n, replace = TRUE)
      c(srho(d$copies[j], d$div[j]), srho(d$copies[j], ds[j])) }, numeric(2))
    dl <- bt[1, ] - bt[2, ]; BOOT[[paste(A, B)]] <- dl
    PAIRS[[paste(A, B)]] <- data.table(focal = A, other = B, n_og = n, rho_species = r_sp, rho_study = r_st,
      delta = r_sp - r_st, lo = quantile(dl, .025)[[1]], hi = quantile(dl, .975)[[1]])
    if (file.exists(within_file(B))) {                     # species-rung attenuation reference (as exploratory E2)
      wB <- fread(within_file(B)); m <- data.table(Orthogroup = d$Orthogroup, copies = d$copies, div = d$div,
        wa = w1f$div[match(d$Orthogroup, w1f$Orthogroup)], wb = wB$div[match(d$Orthogroup, wB$Orthogroup)])
      m <- m[is.finite(wa) & is.finite(wb)]
      sp_oe[B] <- r_sp / atten(m$copies, (m$wa + m$wb) / 2, m$div)$rho_exp
    }
  }
  DESC[[A]][, `:=`(species_obs_over_exp_median = median(sp_oe), species_obs_over_exp_n = length(sp_oe))]
}
PAIRS <- rbindlist(PAIRS); DESC <- rbindlist(DESC)
fwrite(PAIRS, file.path(OUT, "pairs_step4.tsv"), sep = "\t")
fwrite(DESC, file.path(OUT, "descriptive_step4.tsv"), sep = "\t")

## ============================================================== 6. Delta and the decision rule
BM <- do.call(cbind, BOOT)                                 # NBOOT x pairs; pairs resampled independently
per_sp <- rbindlist(lapply(FOCAL, function(A) { k <- which(PAIRS$focal == A); b <- rowMeans(BM[, k, drop = FALSE])
  data.table(level = A, n_pairs = length(k), delta = mean(PAIRS$delta[k]), lo = quantile(b, .025)[[1]], hi = quantile(b, .975)[[1]]) }))
ball <- rowMeans(sapply(FOCAL, function(A) rowMeans(BM[, PAIRS$focal == A, drop = FALSE])))
overall <- data.table(level = "all", n_pairs = nrow(PAIRS), delta = mean(per_sp$delta),
                      lo = quantile(ball, .025)[[1]], hi = quantile(ball, .975)[[1]])
verdict <- function(est, lo) if (est >= MATERIAL && lo > 0) "between-species gradient exceeds between-study" else
  if (abs(est) <= NEGLIGIBLE) "between-species gradient no larger than between-study (study account)" else "partial"
TESTS <- rbind(per_sp, overall)
TESTS[, verdict := mapply(verdict, delta, lo)]
TESTS[, flag := fifelse(delta < -NEGLIGIBLE, "NOTE: between-study gradient is the larger one", "")]
fwrite(TESTS, file.path(OUT, "tests_step4.tsv"), sep = "\t")

say("\n== Pairs (rho on the Step 2 orthogroup set of each pair; interval: %d orthogroup bootstraps)", NBOOT)
for (r in PAIRS[, sprintf("  %-4s vs %-4s  n %5d  rho_species %+.4f  rho_study %+.4f  delta %+.4f [%+.4f, %+.4f]",
                          focal, other, n_og, rho_species, rho_study, delta, lo, hi)]) say(r)
say("\n== Delta (section 5) and decision (section 6)")
for (r in TESTS[, sprintf("  %-4s  pairs %2d  delta %+.4f [%+.4f, %+.4f]%s  %s", level, n_pairs, delta, lo, hi,
                          fifelse(level == "all", paste0("  -> ", verdict, "   (PRIMARY)"), paste0("  (", verdict, ")")), flag)]) say(r)
say("\n== Descriptive (section 7; no decision rule)")
for (r in DESC[, sprintf("  %-4s  samples %d + %d  orthogroups %d  rho_study (all) %+.4f | half-split study 1 %+.4f | study 2 %+.4f",
                         species, n_study1, n_study2, n_og, rho_study_all, rho_within_study1, rho_within_study2)]) say(r)
for (r in DESC[, sprintf("  %-4s  attenuation reference, observed/expected:  study rung %.2f  |  species rung %.2f (median of %d pairs)",
                         species, study_obs_over_exp, species_obs_over_exp_median, species_obs_over_exp_n)]) say(r)
say("\nwrote %s", OUT)
close(con)
