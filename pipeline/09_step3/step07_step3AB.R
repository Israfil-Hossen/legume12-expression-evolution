#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Step 3, analyses A and B (docs/step3_spec.md). Details the spec leaves open are
## marked IMPLEMENTATION.
##  A  paralog-pair expression divergence D, its half-split noise floor N, and D - N, against Ks
##  B  Brownian-motion rate per orthogroup with and without within-species variance
## Inputs: Step 2 QC table and quantification, OrthoFinder orthogroups, step3/<sp>/pairs_ks.tsv,
##         tree/concat.fa.treefile.  Output: <project>/results_step3/
## Usage: Rscript step07_step3AB.R [project_dir]
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(ape) })
args <- commandArgs(TRUE)
P <- if (length(args)) args[1] else "/home/liuchuanfa/israfil/legume2"
OUT <- file.path(P, "results_step3"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20260924)
NDRAW <- 50; NBOOT <- 2000; MATERIAL <- 0.050; NEGLIGIBLE <- 0.025; KS_MAX <- 3
SPECIES <- c("Gmax","Pvul","Vrad","Vung","Ccaj","Mtru","Cari","Psat","Tpra","Ljap","Lang","Adur")
con <- file(file.path(OUT, "report_step3AB.txt"), open = "wt")
say <- function(fmt, ...) { s <- if (length(list(...))) sprintf(fmt, ...) else fmt
  cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }
srho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
say("Step 3 analyses A and B   %s   project: %s", format(Sys.time()), P)

## ---------------------------------------------------------------- shared inputs (as Step 2)
S <- fread(file.path(P, "results_step2", "qc_samples.tsv"))[qc_pass == TRUE]
ogf <- tail(sort(Sys.glob(file.path(P, "proteomes/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv"))), 1)
og <- fread(ogf, sep = "\t", header = TRUE, colClasses = "character")
long <- rbindlist(lapply(SPECIES, function(sp) {
  v <- og[[sp]]; v[is.na(v)] <- ""; g <- strsplit(v, ", ", fixed = TRUE)
  data.table(Orthogroup = rep(og$Orthogroup, lengths(g)), species = sp, gene = sub(paste0("^", sp, "_"), "", unlist(g)))
}))
cn <- long[, .(copies = .N), by = .(Orthogroup, species)]
cnw <- dcast(cn, Orthogroup ~ species, value.var = "copies", fill = 0)
SC <- cnw[rowSums(as.matrix(cnw[, ..SPECIES]) == 1) == length(SPECIES), Orthogroup]

## gene TPM matrices (genes in orthogroups x QC-passed samples), per species
G <- list()
for (sp in SPECIES) {
  bs <- S[species == sp, biorep]; if (!length(bs)) next
  t2g <- fread(file.path(P, "tables", sprintf("%s_tx2gene.tsv", sp)), header = FALSE, col.names = c("Name", "gene"))
  genes <- long[species == sp, unique(gene)]
  M <- matrix(0, length(genes), length(bs), dimnames = list(genes, bs))
  for (b in bs) {
    q <- fread(file.path(P, "quant", b, "quant.sf"), select = c("Name", "TPM"))
    x <- merge(q, t2g, by = "Name")[, .(TPM = sum(TPM)), by = gene][gene %in% genes]
    M[x$gene, b] <- x$TPM
  }
  G[[sp]] <- M
}
## sample scaling exactly as Step 2: orthogroup sums, median over single-copy orthogroups,
## each sample scaled to the grand median of those medians
og_sum <- function(sp) { M <- G[[sp]]; L <- long[species == sp]
  rowsum(M[L$gene, , drop = FALSE], L$Orthogroup) }
OGS <- lapply(setNames(names(G), names(G)), og_sum)
sf <- unlist(unname(lapply(OGS, function(M) apply(M[intersect(SC, rownames(M)), , drop = FALSE], 2, median))))
C <- median(sf)
say("samples: %d QC-passed; single-copy orthogroups used for scaling: %d", nrow(S), length(SC))

split_draw <- function(groups, k) {                       # identical to Step 2
  gs <- split(names(groups), groups); cap <- floor(lengths(gs) / 2)
  if (sum(cap) < k) return(NULL)
  al <- setNames(integer(length(gs)), names(gs))
  for (i in seq_len(k)) { av <- names(cap)[cap - al > 0]; pk <- av[sample.int(length(av), 1)]; al[pk] <- al[pk] + 1L }
  h1 <- h2 <- character(0)
  for (g in names(al)[al > 0]) { s <- gs[[g]][sample.int(length(gs[[g]]), 2 * al[g])]
    h1 <- c(h1, s[seq_len(al[g])]); h2 <- c(h2, s[al[g] + seq_len(al[g])]) }
  list(h1, h2)
}
boot_median <- function(v) { v <- v[is.finite(v)]
  bt <- replicate(NBOOT, median(sample(v, length(v), replace = TRUE)))
  c(est = median(v), lo = unname(quantile(bt, .025)), hi = unname(quantile(bt, .975))) }

## ================================================================ A. paralog pairs
say("\n== A. Paralog expression divergence against Ks (Ks <= %g; pairs with both genes mean TPM < 1 excluded)", KS_MAX)
A <- list(); PAIRS <- list()
for (sp in names(G)) {
  f <- file.path(P, "step3", sp, "pairs_ks.tsv"); if (!file.exists(f)) { say("  %-5s no pairs_ks.tsv - skipped", sp); next }
  pr <- fread(f)[, `:=`(g1 = sub(paste0("^", sp, "_"), "", gene1), g2 = sub(paste0("^", sp, "_"), "", gene2))]
  M <- G[[sp]]
  pr <- pr[is.finite(ks) & ks <= KS_MAX & g1 %in% rownames(M) & g2 %in% rownames(M)]
  ## IMPLEMENTATION: the expression filter uses unscaled TPM averaged over QC-passed samples
  pr <- pr[!(rowMeans(M[g1, , drop = FALSE]) < 1 & rowMeans(M[g2, , drop = FALSE]) < 1)]
  X <- log2(sweep(M, 2, sf[colnames(M)], "/") * C + 1)
  Z <- X[pr$g1, , drop = FALSE] - X[pr$g2, , drop = FALSE]
  pr[, D := abs(rowMeans(Z))]
  grp <- setNames(S[match(colnames(M), biorep), group], colnames(M))
  if (ncol(M) >= 4) {
    acc <- 0; accP <- 0
    for (d in seq_len(NDRAW)) { h <- split_draw(grp, 2)
      m1 <- rowMeans(Z[, h[[1]], drop = FALSE]); m2 <- rowMeans(Z[, h[[2]], drop = FALSE])
      acc <- acc + abs(m1 - m2) / sqrt(2)
      accP <- accP + m1 * m2 }                            # deviation D2: E[m1*m2] = true squared divergence
    pr[, `:=`(N = acc / NDRAW, P = accP / NDRAW)]
  } else pr[, `:=`(N = NA_real_, P = NA_real_)]
  pr[, excess := D - N]
  PAIRS[[sp]] <- pr[, .(species = sp, gene1, gene2, ks, collinear, D, N, excess, P)]
  A[[sp]] <- data.table(species = sp, n_pairs = nrow(pr), rho_D_ks = srho(pr$D, pr$ks),
                        rho_N_ks = srho(pr$N, pr$ks), rho_excess_ks = srho(pr$excess, pr$ks),
                        r_P_ks = if (all(is.na(pr$P))) NA_real_ else cor(pr$P, pr$ks))
}
A <- rbindlist(A); PAIRS <- rbindlist(PAIRS)
fwrite(A, file.path(OUT, "A_per_species.tsv"), sep = "\t")
fwrite(PAIRS, file.path(OUT, "A_pairs.tsv.gz"), sep = "\t")
for (r in A[, sprintf("  %-5s pairs %5d   rho(D,Ks) %+.3f   rho(N,Ks) %+.3f   rho(D-N,Ks) %+.3f   r(m1*m2,Ks) %+.3f",
                     species, n_pairs, rho_D_ks, rho_N_ks, rho_excess_ks, r_P_ks)]) say(r)
AN <- A[is.finite(rho_N_ks)]                               # species with a half-split (>= 4 samples)
a1 <- boot_median(AN$rho_D_ks); a2 <- boot_median(AN$rho_N_ks)
a3reg <- boot_median(AN$rho_excess_ks)                    # as registered; biased (deviation D2)
a3 <- boot_median(AN$r_P_ks)                              # deviation D2: Pearson r(mean m1*m2, Ks)
v2 <- if (a2["est"] <= -MATERIAL && a2["hi"] < 0) "ambiguity noise confirmed" else "not confirmed"
v3 <- if (a3["est"] >= MATERIAL && a3["lo"] > 0) "evolutionary divergence detectable at gene level" else
      if (a3["est"] <= NEGLIGIBLE) "not detectable" else "partial"
say("  median over %d species (95%% interval, species bootstrap):", nrow(AN))
say("  A1 rho(D, Ks)     %+.4f [%+.4f, %+.4f]   (reported)", a1[1], a1[2], a1[3])
say("  A2 rho(N, Ks)     %+.4f [%+.4f, %+.4f]   -> %s", a2[1], a2[2], a2[3], v2)
say("  A3 as registered, rho(D-N, Ks)  %+.4f [%+.4f, %+.4f]   BIASED - reported only (deviation D2)", a3reg[1], a3reg[2], a3reg[3])
say("  A3 r(mean m1*m2, Ks)             %+.4f [%+.4f, %+.4f]   -> %s   (PRIMARY, deviation D2)", a3[1], a3[2], a3[3], v3)
## exploratory: collinear (WGD/segmental) versus other pairs within Ks bins
EX <- PAIRS[is.finite(excess)][, `:=`(class = fifelse(collinear == 1, "WGD/segmental", "other"),
        ks_bin = cut(ks, c(0, 0.1, 0.3, 0.6, 1, 3), include.lowest = TRUE))][
        , .(n = .N, median_D_minus_N = median(excess)), by = .(ks_bin, class)][order(ks_bin, class)]
fwrite(EX, file.path(OUT, "A_exploratory_class_by_ks.tsv"), sep = "\t")
say("  EXPLORATORY (no decision rule): median D - N by duplicate class within Ks bins")
for (r in EX[, sprintf("    Ks %-10s %-14s n %6d   D-N %+.3f", ks_bin, class, n, median_D_minus_N)]) say(r)

## ================================================================ B. Brownian-motion rate
tr <- read.tree(file.path(P, "tree", "concat.fa.treefile"))
tr <- drop.tip(tr, setdiff(tr$tip.label, names(G)))
sp_b <- intersect(tr$tip.label, names(G))
Cm <- vcv(tr)[sp_b, sp_b]
allog <- cnw[rowSums(as.matrix(cnw[, ..sp_b]) >= 1) == length(sp_b), Orthogroup]
say("\n== B. Brownian-motion rate per orthogroup (%d species on the IQ-TREE tree; %d orthogroups present in all)",
    length(sp_b), length(allog))
mean_se <- lapply(setNames(sp_b, sp_b), function(sp) {
  M <- OGS[[sp]][allog, , drop = FALSE]
  cp <- cn[species == sp][match(allog, Orthogroup), copies]
  X <- log2(sweep(M, 2, sf[colnames(M)], "/") * C / cp + 1)     # per-copy expression, log2
  list(m = rowMeans(X), se2 = apply(X, 1, var) / ncol(X))
})
Ym <- sapply(mean_se, `[[`, "m"); Ys <- sapply(mean_se, `[[`, "se2")
bm_fit <- function(y, se2) {                              # ML sigma^2 with GLS root; se2 on the diagonal
  nll <- function(ls2) { V <- exp(ls2) * Cm + diag(se2, length(y))
    ch <- tryCatch(chol(V), error = function(e) NULL); if (is.null(ch)) return(1e10)
    Vi <- chol2inv(ch); a <- sum(Vi %*% y) / sum(Vi); r <- y - a
    sum(log(diag(ch))) + 0.5 * drop(t(r) %*% Vi %*% r) }
  exp(optimize(nll, c(-20, 5))$minimum)
}
s2_naive <- sapply(seq_along(allog), function(i) bm_fit(Ym[i, ], rep(0, length(sp_b))))
s2_me    <- sapply(seq_along(allog), function(i) bm_fit(Ym[i, ], Ys[i, ]))
pred <- rowMeans(log2(as.matrix(cnw[match(allog, Orthogroup), ..sp_b])))
BT <- data.table(Orthogroup = allog, mean_log2_copies = pred, s2_naive = s2_naive, s2_me = s2_me)
fwrite(BT, file.path(OUT, "B_rates.tsv.gz"), sep = "\t")
b1 <- srho(s2_naive, pred); b_me <- srho(s2_me, pred); dr <- b_me - b1
bt <- replicate(NBOOT, { i <- sample.int(length(pred), replace = TRUE); srho(s2_me[i], pred[i]) - srho(s2_naive[i], pred[i]) })
ci <- quantile(bt, c(.025, .975))
vb <- if (dr <= -MATERIAL) "within-species variance accounts for a material part of the gradient (technical)" else
      if (abs(dr) <= NEGLIGIBLE) "no change (evolutionary)" else "partial"
say("  B1 rho(sigma2_naive, copies) %+.4f   (reported)", b1)
say("     rho(sigma2_me,    copies) %+.4f", b_me)
say("  B2 delta rho %+.4f [%+.4f, %+.4f]   -> %s   (PRIMARY)", dr, ci[1], ci[2], vb)
say("  orthogroups whose fitted rate reached the lower bound: naive %d, with variance %d",
    sum(s2_naive < 1e-8), sum(s2_me < 1e-8))

fwrite(data.table(analysis = c("A1", "A2", "A3_registered_biased", "A3_D2", "B1", "B2"),
                  estimate = c(a1[1], a2[1], a3reg[1], a3[1], b1, dr), lo = c(a1[2], a2[2], a3reg[2], a3[2], NA, ci[1]),
                  hi = c(a1[3], a2[3], a3reg[3], a3[3], NA, ci[2]),
                  verdict = c("reported", v2, "not used (biased, D2)", v3, "reported", vb)),
       file.path(OUT, "tests_step3AB.tsv"), sep = "\t")
say("\nwrote %s", OUT)
close(con)
