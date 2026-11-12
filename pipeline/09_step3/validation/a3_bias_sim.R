## Pure noise: two paralogs with NO true expression difference; noise sd falls with Ks (A2's premise).
set.seed(1); n <- 10; npair <- 3000; ndraw <- 50
ks <- runif(npair, 0.01, 3); sd_z <- 0.6 * exp(-ks) + 0.15            # young pairs noisier
Z <- matrix(rnorm(npair * n, 0, rep(sd_z, n)), npair)                    # gene1 - gene2 per sample, true diff = 0
D_reg <- abs(rowMeans(Z))                                                 # as registered: all samples
acc_N <- 0; acc_D2 <- 0
for (d in 1:ndraw) { s <- sample(n, 4); h1 <- s[1:2]; h2 <- s[3:4]
  m1 <- rowMeans(Z[, h1]); m2 <- rowMeans(Z[, h2])
  acc_N <- acc_N + abs(m1 - m2) / sqrt(2); acc_D2 <- acc_D2 + (abs(m1) + abs(m2)) / 2 }
N <- acc_N / ndraw; D2 <- acc_D2 / ndraw
r <- function(x) cor(x, ks, method = "spearman")
cat(sprintf("true divergence = 0 for every pair\n"))
cat(sprintf("registered A3: rho(D_all - N, Ks)       = %+.3f   <- spurious 'evolutionary' signal\n", r(D_reg - N)))
cat(sprintf("corrected A3:  rho(D_2samples - N, Ks)  = %+.3f   <- ~0, as it should be\n", r(D2 - N)))
cat(sprintf("A2: rho(N, Ks)                           = %+.3f\n", r(N)))
