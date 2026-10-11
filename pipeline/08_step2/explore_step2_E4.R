#!/usr/bin/env Rscript
## EXPLORATORY (post hoc) E4: does rho_between track the two species' own within-species
## gradients (a species-level measurement property) better than it tracks divergence time?
suppressPackageStartupMessages(library(data.table))
R <- "F:/legume_work/results/results_step2"; OUT <- "F:/legume_work/results/exploratory_step2"
X <- fread(file.path(OUT, "pairs_dilution_attenuation.tsv"))
W <- fread(file.path(R, "rho_within_primary.tsv"))[, .(species, rho_w = rho)]
X <- merge(X, W[, .(sp_a = species, wa = rho_w)], by = "sp_a")
X <- merge(X, W[, .(sp_b = species, wb = rho_w)], by = "sp_b")
X[, mean_within := (wa + wb) / 2]
M <- X[sp_a != "Ccaj" & sp_b != "Ccaj"]
srho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
con <- file(file.path(OUT, "report_E4.txt"), "wt")
say <- function(fmt, ...) { s <- sprintf(fmt, ...); cat(s, "\n"); cat(s, "\n", file = con) }
say("E4 (EXPLORATORY): %d pairs with a half-split in both species, Cajanus excluded", nrow(M))
say("  Spearman rho(rho_between, mean within rho of the two species) = %+.3f", srho(M$rho_obs, M$mean_within))
say("  Spearman rho(rho_between, divergence time)                   = %+.3f", srho(M$rho_obs, M$time))
f1 <- lm(rho_obs ~ time, M); f2 <- lm(rho_obs ~ mean_within, M); f3 <- lm(rho_obs ~ time + mean_within, M)
say("  R2: time alone %.3f | mean within alone %.3f | both %.3f", summary(f1)$r.squared, summary(f2)$r.squared, summary(f3)$r.squared)
cf <- coef(summary(f3))
say("  joint model: time %+.5f per Myr (P = %.3g) -> change over %.1f Myr = %+.3f;  mean within %+.3f (P = %.3g)",
    cf["time", 1], cf["time", 4], diff(range(M$time)), cf["time", 1] * diff(range(M$time)), cf["mean_within", 1], cf["mean_within", 4])
close(con)
