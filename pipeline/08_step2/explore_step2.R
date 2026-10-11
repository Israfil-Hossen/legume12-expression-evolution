#!/usr/bin/env Rscript
## EXPLORATORY (post hoc, not pre-registered) follow-up of the Step 2 results.
## Reads the registered per-orthogroup outputs in results/results_step2; changes nothing there.
##  E1 dilution: does between-species divergence grow with divergence time?
##  E2 attenuation reference (manuscript's method): expected rho_between if the part of
##     between-species divergence beyond within-species noise carries no copy-number gradient
##  E3 within-species rho vs read depth and replicate number
suppressPackageStartupMessages(library(data.table))
R  <- "F:/legume_work/results/results_step2"
OUT <- "F:/legume_work/results/exploratory_step2"; dir.create(OUT, showWarnings = FALSE)
con <- file(file.path(OUT, "report_exploratory.txt"), "wt")
say <- function(fmt, ...) { s <- if (length(list(...))) sprintf(fmt, ...) else fmt; cat(s, "\n"); cat(s, "\n", file = con) }
say("EXPLORATORY analyses (post hoc; not part of the pre-registration)  %s", format(Sys.time()))
srho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))

TT <- fread("F:/legume_work/metadata/step2_divergence_times.tsv")[, .(sp_a, sp_b, time = median_mya)]
TT <- rbind(TT, TT[, .(sp_a = sp_b, sp_b = sp_a, time)])
B  <- merge(fread(file.path(R, "rho_between_primary.tsv")), TT, by = c("sp_a", "sp_b"))
W  <- fread(file.path(R, "rho_within_primary.tsv"))
within <- lapply(setNames(W$species, W$species), function(s) fread(file.path(R, sprintf("div_within_primary_%s.tsv.gz", s))))

rows <- list()
for (i in seq_len(nrow(B))) {
  a <- B$sp_a[i]; b <- B$sp_b[i]
  d <- fread(file.path(R, sprintf("div_between_primary_%s_%s.tsv.gz", a, b)))
  r <- data.table(sp_a = a, sp_b = b, time = B$time[i], rho_obs = B$rho[i], median_div = median(d$div))
  if (!is.null(within[[a]]) && !is.null(within[[b]])) {
    m <- merge(d, within[[a]][, .(Orthogroup, wa = div)], by = "Orthogroup")
    m <- merge(m, within[[b]][, .(Orthogroup, wb = div)], by = "Orthogroup")
    m[, noise := (wa + wb) / 2]
    rho_noise <- srho(m$copies, m$noise)
    r[, `:=`(rho_noise = rho_noise, sd_ratio = sd(m$noise) / sd(m$div),
             rho_exp = rho_noise * sd(m$noise) / sd(m$div), n_og = nrow(m))]
    r[, obs_over_exp := rho_obs / rho_exp]
  }
  rows[[i]] <- r
}
X <- rbindlist(rows, fill = TRUE)
fwrite(X, file.path(OUT, "pairs_dilution_attenuation.tsv"), sep = "\t")

main <- X[sp_a != "Ccaj" & sp_b != "Ccaj"]          # same pair set as the registered primary test (D1)
say("\n== E1 dilution: median between-species divergence vs divergence time (%d pairs, Cajanus excluded)", nrow(main))
say("  Spearman rho(median divergence, time) = %+.3f", srho(main$median_div, main$time))
say("  Spearman rho(rho_between, median divergence) = %+.3f", srho(main$rho_obs, main$median_div))
for (r in main[, .(n = .N, median_div = median(median_div), rho_between = median(rho_obs)), by = .(time = round(time, 1))][order(time)][
       , sprintf("  time %5.1f Myr  pairs %2d  median divergence %.3f  median rho_between %+.4f", time, n, median_div, rho_between)]) say(r)

E <- main[is.finite(obs_over_exp)]
say("\n== E2 attenuation reference (%d pairs with a half-split in both species; Lotus has none)", nrow(E))
say("  median rho_between observed %+.4f | expected if extra component copy-neutral %+.4f", median(E$rho_obs), median(E$rho_exp))
say("  observed / expected: median %.2f  (IQR %.2f - %.2f)", median(E$obs_over_exp), quantile(E$obs_over_exp, .25), quantile(E$obs_over_exp, .75))
fit <- lm(obs_over_exp ~ time, E)
say("  obs/exp ~ time: slope %+.4f per Myr  ->  change over %.1f-%.1f Myr = %+.3f",
    coef(fit)[[2]], min(E$time), max(E$time), coef(fit)[[2]] * diff(range(E$time)))
say("  Spearman rho(obs/exp, time) = %+.3f", srho(E$obs_over_exp, E$time))

Q <- fread(file.path(R, "qc_samples.tsv"))[qc_pass == TRUE, .(median_reads_M = median(reads_after) / 1e6, n = .N), by = species]
Z <- merge(W[, .(species, rho)], Q, by = "species")
fwrite(Z, file.path(OUT, "within_rho_vs_depth.tsv"), sep = "\t")
say("\n== E3 within-species rho vs read depth and replicates (11 species; descriptive only)")
for (r in Z[order(rho)][, sprintf("  %-5s rho %+.4f   median reads %5.1f M   n %2d", species, rho, median_reads_M, n)]) say(r)
say("  Spearman rho(within rho, depth) = %+.3f;  rho(within rho, n) = %+.3f", srho(Z$rho, Z$median_reads_M), srho(Z$rho, Z$n))
close(con)
