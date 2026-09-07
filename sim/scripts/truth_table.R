# ------------------------------------------------------------------------------
# Analytic truth table for the M3a grid (SAP item P5).
#
# Truth for dimension C is graded ANALYTICALLY from the marginal true
# survival, not from grid-cell membership:
#   theta_C = (RMST(horizon) - RMST(t_cutoff)) / RMST(horizon)
# where t_cutoff solves S_marg(t_cutoff) = 1 - maturity (the administrative
# cutoff used by generate_scenario). theta_C is the share of the decision
# window's restricted mean survival that lies beyond data support.
#
# Grading (aligned with psa_config defaults rho_maturity_yellow/red):
#   theta_C <  0.20 -> none   (expected: green)
#   0.20 <= theta_C < 0.50 -> mild  (expected: yellow / PARTIAL_MATURITY)
#   theta_C >= 0.50 -> severe      (expected: red / IMMATURE_DATA)
# ------------------------------------------------------------------------------

pkg_root <- sub("sim.*$", "", getwd())
source(file.path(pkg_root, "sim", "R", "generate.R"))

lambda0 <- log(2) / 12
horizon <- 36

s_marg <- function(t, shape) {
  0.5 * (exp(-lambda0 * t) + exp(-lambda0 * cumHR(t, shape)))
}
rmst <- function(a, shape) integrate(s_marg, 0, a, shape = shape)$value

grid <- expand.grid(shape = shape_levels, maturity = c(0.4, 0.6, 0.8))
grid$t_cutoff <- mapply(function(sh, m) {
  uniroot(function(t) s_marg(t, sh) - (1 - m), c(0, 240))$root
}, grid$shape, grid$maturity)
grid$rmst_horizon <- mapply(rmst, horizon, grid$shape)
grid$rmst_cutoff  <- mapply(rmst, grid$t_cutoff, grid$shape)
grid$theta_C <- with(grid, (rmst_horizon - rmst_cutoff) / rmst_horizon)
grid$truth_C_grade <- cut(grid$theta_C, c(-Inf, 0.20, 0.50, Inf),
                          labels = c("none", "mild", "severe"))
# Expected S_marg at the median-event tau (tau ~ 8-9 months; report at the
# analytic median of the pooled event distribution as a reference point)
grid$t_median <- mapply(function(sh) {
  uniroot(function(t) s_marg(t, sh) - 0.5, c(0, 240))$root
}, grid$shape)

print(grid, digits = 3)
outdir <- file.path(pkg_root, "sim", "outputs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(grid, file.path(outdir, "truth_table.csv"), row.names = FALSE)
cat("\nWritten:", file.path(outdir, "truth_table.csv"), "\n")
cat("Kappa <-> Kendall tau (Gaussian copula): tau = kappa exactly.\n")
