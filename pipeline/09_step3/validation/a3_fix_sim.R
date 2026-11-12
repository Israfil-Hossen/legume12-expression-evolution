## Compare estimators of "divergence in excess of noise" when noise falls with Ks (A2's premise).
sim <- function(delta_fun, label, n = 10, npair = 4000, ndraw = 50, reps = 20) {
  out <- replicate(reps, {
    ks <- runif(npair, 0.01, 3); sd_z <- 0.6 * exp(-ks) + 0.15
    Z <- matrix(rnorm(npair * n, rep(delta_fun(ks), n), rep(sd_z, n)), npair)
    D_all <- abs(rowMeans(Z)); aN <- aD2 <- aP <- 0
    for (d in 1:ndraw) { s <- sample(n, 4); m1 <- rowMeans(Z[, s[1:2]]); m2 <- rowMeans(Z[, s[3:4]])
      aN <- aN + abs(m1 - m2) / sqrt(2); aD2 <- aD2 + (abs(m1) + abs(m2)) / 2; aP <- aP + m1 * m2 }
    N <- aN / ndraw; r <- function(x) cor(x, ks, method = "spearman")
    c(registered = r(D_all - N), product_spearman = r(aP / ndraw), product_pearson = cor(aP / ndraw, ks))
  })
  m <- rowMeans(out)
  cat(sprintf("%-38s registered %+.3f | product Spearman %+.3f | product PEARSON %+.3f\n", label, m[1], m[2], m[3]))
}
set.seed(2)
sim(function(k) 0,             "no true divergence (answer 0)")
sim(function(k) 0.5,           "equal divergence, all ages (answer 0)")
sim(function(k) 0.25 * k,      "divergence grows with Ks (answer > 0)")
