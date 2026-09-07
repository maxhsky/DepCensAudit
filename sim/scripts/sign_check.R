# Sign verification for SAP revision C1 (Yan Wuxia).
# Derivation: T = -log(U1)/lambda is DECREASING in U1, C = qweibull(U2) is
# increasing in U2, and the Gaussian copula couples (U1, U2) with
# tau(U1, U2) = kappa. Hence tau(T, C) = -kappa exactly.
# This script verifies empirically at the generator level.

pkg_root <- sub("sim.*$", "", getwd())
source(file.path(pkg_root, "sim", "R", "generate.R"))

set.seed(20260908)
n <- 2e4  # Kendall is O(n^2): keep n small (2e4 -> ~2e8 pairs, seconds)
for (k in c(0.3, 0.6)) {
  U1 <- runif(n)
  rho <- sin(pi * k / 2)
  z1 <- qnorm(U1)
  U2 <- pnorm(rho * z1 + sqrt(1 - rho^2) * rnorm(n))
  Tfull <- -log(U1) / .lambda0
  Cfull <- qweibull(U2, shape = .cens_shape, scale = .cens_scale)
  tau_tc <- cor(Tfull, Cfull, method = "kendall")
  cat(sprintf("kappa = %.1f : empirical tau(T, C) = %+.3f (predicted -%.1f)\n",
              k, tau_tc, k))
}
