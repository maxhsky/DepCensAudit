pkg_root <- sub("sim.*$", "", getwd())
for (f in c("sim/R/generate.R", "sim/R/runners.R", "sim/R/analysis.R"))
  source(file.path(pkg_root, f), local = TRUE)

# pools build correctly
stopifnot(length(.pool_spec_e2) == 4, length(.pool_severe) == 11)
stopifnot(all(grepl("^G_.*m0\\.8$", .pool_spec_e2)),
          all(.pool_severe[1:9] != ""),
          sum(grepl("^OS_", .pool_severe)) == 2)
cat("pools OK: spec_e2", length(.pool_spec_e2), "| severe", length(.pool_severe), "\n")
print(.pool_severe)

# E2 three-dimension non-green + errored-in-denominator, synthetic rows
fake <- data.frame(
  scenario_id = rep("G_ph_k0.0_m0.8", 6),
  eval_flag_codes = c("", "", "IPCW_WEIGHT_INSTABILITY", "", "", ""),
  pos_flag_codes  = c("", "", "", "SPARSE_TAIL_RISKSET", "", ""),
  pos_worst_level = c("green", "green", "green", "yellow", "green", NA),
  ext_perf_level  = c("green", "yellow", "green", "green", "green", NA),
  ext_root_codes  = c("", "", "", "", "EXTRAP_INSTABILITY", NA),
  error = c("none", "none", "none", "none", "none", "flexsurv failed"),
  truth_C = FALSE, kappa = 0, stringsAsFactors = FALSE)
res <- e2_clean_scene_rate(fake)
stopifnot(res$n_clean == 6, res$n_flagged == 4)  # rows 1 and 6 unflagged
                                                  # (row 6 = errored, in denom)
stopifnot(abs(res$rate_valid_only - 4/5) < 1e-12) # error-free rows only
# 1/6 errors = 16.7% > 2% bound -> downgrade clause fires (by design)
stopifnot(res$judgement_status == "downgraded_descriptive_error_rate_gt_2pct")
cat("E2 OK: n =", res$n_clean, "| flagged =", res$n_flagged,
    "| rate =", round(res$rate, 3), "| rate_valid =",
    round(res$rate_valid_only, 3), "| status:", res$judgement_status, "\n")

# detect() + sens_spec C_mat dimension
fake2 <- data.frame(
  scenario_id = "G_ph_k0.0_m0.4", rep = 1:4,
  error = c("none", "none", "boom", "none"),
  truth_B = TRUE, truth_C = TRUE, truth_A = FALSE,
  ext_perf_level = c("green", "yellow", NA, "green"),
  ext_perf_code = c("ALL_CLEAR", "PARTIAL_MATURITY", "", "ALL_CLEAR"),
  ext_root_codes = c("", "PARTIAL_MATURITY", "", ""),
  pos_flag_codes = c("", "", "", ""),
  eval_weight_unstable = c(0, 0, NA, 0))
d <- detect(fake2)
stopifnot(sum(d$det_C_mat) == 1, sum(d$det_C) == 1,
          all(!d$det_C[3]), d$errored[3])
ss <- sens_spec(d, "C_mat")
stopifnot(ss$sensitivity == 0.25, ss$error_rate == 0.25)  # 1 of 4 (error = miss)
cat("detect/sens_spec C_mat OK: sens =", ss$sensitivity,
    "| error_rate =", ss$error_rate, "\n")
cat("ALL V1.2 CHECKS PASSED\n")
