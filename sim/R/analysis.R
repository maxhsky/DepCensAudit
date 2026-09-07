# ------------------------------------------------------------------------------
# M3a simulation framework :: analysis helpers (SAP v1.1 estimands)
#
# FROZEN-IN-CANDIDATE - detection rules and judgement criteria are frozen in
# the SAP before any performance readout (project doc v1.0 r2, S5; SAP v1.1
# revision list R1-R9). Errors count as NON-DETECTED (denominator retained;
# per-cell error rates reported alongside) - no silent na.rm dropping.
# ------------------------------------------------------------------------------

library(survival)

is_err <- function(rows) !is.na(rows$error) & rows$error != "none"

# Detection rules (SAP section 4; consume TOOL OUTPUT ONLY):
#   det_C      : extrapolation performance flag non-green OR root code
#                matches IMMATURE|INSTAB   (SPARSE branch removed: R5 -
#                SPARSE_TAIL_RISKSET is a positivity-dimension code and
#                never appears in ext_root_codes)
#   det_C_mat  : maturity-coded subset (IMMATURE|MATURITY codes only)
#   det_C_inst : instability-coded subset (INSTAB codes only), reported
#                separately to keep the E1 sensitivity semantics clean
#   det_B      : positivity root flag non-empty (SPARSE_TAIL_RISKSET)
#   det_A      : evaluation weight-instability flag raised
# Reps with error != "none" count as non-detected (R4/B3).
detect <- function(rows) {
  err <- is_err(rows)
  perf_ng <- !is.na(rows$ext_perf_level) & rows$ext_perf_level != "green"
  data.frame(
    scenario_id = rows$scenario_id, rep = rows$rep, errored = err,
    truth_B = rows$truth_B, truth_C = rows$truth_C, truth_A = rows$truth_A,
    det_C = !err & (perf_ng | grepl("IMMATURE|INSTAB", rows$ext_root_codes)),
    det_C_mat = !err & (grepl("IMMATURE|MATURITY", rows$ext_root_codes) |
                          rows$ext_perf_code %in% c("PARTIAL_MATURITY",
                                                    "IMMATURE_DATA")),
    det_C_inst = !err & (grepl("INSTAB", rows$ext_root_codes) |
                           rows$ext_perf_code == "EXTRAP_INSTABILITY"),
    det_B = !err & (!is.na(rows$pos_flag_codes) & nzchar(rows$pos_flag_codes)),
    det_A = !err & !is.na(rows$eval_weight_unstable) &
      rows$eval_weight_unstable == 1)
}

# Per-scenario sensitivity / specificity for one dimension. Errors stay in
# the denominator as non-detections (SAP B3). Dimension "C_mat" consumes
# det_C_mat - the E1 PRIMARY judgement statistic (SAP v1.2 blocker 5).
# NOTE: the truth_* columns are coarse grid booleans (e.g. truth_C =
# maturity == 0.4), NOT the analytic grade; they are only valid inside
# whitelist pools whose analytic grades are homogeneous (SAP v1.2 r2).
sens_spec <- function(d, dimension = c("B", "C", "C_mat", "A")) {
  dimension <- match.arg(dimension)
  det <- switch(dimension, B = d$det_B, C = d$det_C,
                C_mat = d$det_C_mat, A = d$det_A)
  tru <- switch(dimension, B = d$truth_B, C = d$truth_C,
                C_mat = d$truth_C, A = d$truth_A)
  tp <- sum(det & tru, na.rm = TRUE);  fn <- sum(!det & tru, na.rm = TRUE)
  fp <- sum(det & !tru, na.rm = TRUE); tn <- sum(!det & !tru, na.rm = TRUE)
  data.frame(dimension = dimension,
             sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA,
             specificity = if (fp + tn > 0) tn / (fp + tn) else NA,
             error_rate = mean(d$errored),
             tp = tp, fp = fp, fn = fn, tn = tn)
}

# ---- Judgement pools (SAP v1.2, blocker 4/5): explicit scenario_id
# whitelists so the readout has ZERO discretionary pooling. The n = 200
# sensitivity layer (S_*) never enters any judgement pool. ----

# E1 sensitivity / E2 clean scene: m0.8 x kappa 0 x 4 shapes (analytic
# truth_C = none only at m0.8; the grid boolean truth_C does NOT encode the
# analytic grade and must not be used for pooling).
.pool_spec_e2 <- c("G_ph_k0.0_m0.8", "G_early_k0.0_m0.8",
                   "G_late_k0.0_m0.8", "G_crossing_k0.0_m0.8")

# E1 severe pool: main grid m0.4 x {ph, early, late} x ALL kappa
# {0, 0.3, 0.6} (9 cells) + severe heterologous sublayer ph cells
# (logit_dropout, comp_miscode; staggered excluded - effective maturity
# between m0.6 and m0.8). n = 11 cells x 500 = 5500.
.pool_severe <- c(
  outer(c("ph", "early", "late"), c("0.0", "0.3", "0.6"),
        function(sh, k) sprintf("G_%s_k%s_m0.4", sh, k)),
  "OS_logit_dropout_ph", "OS_comp_miscode_ph")

subset_pool <- function(rows, ids) rows[rows$scenario_id %in% ids, ]

# E2 (false-positive control, SAP v1.2 blockers 1/2/3):
#   scene    = .pool_spec_e2 (truth_C none AND kappa 0, whitelist);
#   non-green = ANY dimension, ANY layer, evaluated over FIVE physical row
#             columns (SAP v1.2 r2, verbatim enumeration):
#             (1) eval_flag_codes non-empty
#             (2) pos_flag_codes non-empty
#             (3) pos_worst_level != "green"
#             (4) ext_perf_level != "green"
#             (5) ext_root_codes non-empty
#             Equivalence note: dimension-A/B root flags are currently
#             yellow-only, so flag existence <=> non-green. This equivalence
#             is guaranteed by the recorder construction in runners.R
#             (.flag_codes / .worst_level are derived from the same flag
#             list); if the recorder changes, re-verify before readout.
#             overall_level alone only aggregates extrapolation performance
#             flags and must NOT be used alone for E2.
#   errors   = retained in the denominator as NOT flagged (frozen rule B3;
#             NOTE the asymmetry vs E1: for E2 this direction is
#             anti-conservative - it eases passage - hence rate_valid and
#             the 2% error-rate downgrade clause below).
#   Execution standard: the PRESET CONSERVATIVE POOLED threshold is the
#   ONLY execution standard (pooled n = 2000 -> rate <= 0.0866). Per-cell
#   non-green rates are reported DESCRIPTIVELY with Wilson 95% CIs;
#   0.072 is a per-cell reference bound only and is NOT a judgement
#   criterion (SAP v1.2 r2 blocker B2). Exact binomial p-value: reference
#   only. Downgrade clause: if the pooled error rate in the spec pool
#   exceeds 2% (preset bound), E2 judgement downgrades to descriptive and
#   any conclusion proceeds via a SAP amendment.
e2_any_nongreen <- function(rows) {
  ng <- (!is.na(rows$eval_flag_codes) & nzchar(rows$eval_flag_codes)) |
    (!is.na(rows$pos_flag_codes) & nzchar(rows$pos_flag_codes)) |
    (!is.na(rows$pos_worst_level) & rows$pos_worst_level != "green") |
    (!is.na(rows$ext_perf_level) & rows$ext_perf_level != "green") |
    (!is.na(rows$ext_root_codes) & nzchar(rows$ext_root_codes))
  ng[is.na(ng)] <- FALSE   # errored / missing rows count as not flagged
  ng
}

e2_clean_scene_rate <- function(rows) {
  d <- subset_pool(rows, .pool_spec_e2)
  err <- is_err(d)
  ng <- e2_any_nongreen(d)                     # errors stay in denominator
  n <- nrow(d); x <- sum(ng)
  p_hat <- x / n
  # Anti-conservative-direction safeguard (SAP v1.2 r2, B3 strengthening):
  # errors ease E2 passage, so report the error-free rate alongside and
  # downgrade to descriptive if pooled error rate exceeds the 2% bound.
  valid <- !err
  rate_valid <- if (any(valid)) sum(ng[valid]) / sum(valid) else NA_real_
  error_rate <- mean(err)
  pval <- tryCatch(
    stats::binom.test(x, n, p = 0.10, alternative = "less")$p.value,
    error = function(e) NA_real_)
  data.frame(n_clean = n, n_flagged = x, rate = p_hat,
             rate_valid_only = rate_valid, error_rate = error_rate,
             exact_p_ref_only = pval,
             pass_threshold_0866 = p_hat <= 0.0866 && error_rate <= 0.02,
             judgement_status = if (error_rate > 0.02)
               "downgraded_descriptive_error_rate_gt_2pct" else "active")
}

# E3(iii): dataset-level OLS of KM bias on kappa with shape fixed effects,
# cluster-robust SE by cell (SAP open-point 2 verdict; Cochran-Armitage
# dropped: KM bias is continuous, CA tests proportions).
e3_bias_trend <- function(rows) {
  ok <- !is_err(rows) & !is.na(rows$km_bias) & !is.na(rows$kappa)
  d <- rows[ok, ]
  d$shape <- factor(d$shape)
  fit <- lm(km_bias ~ kappa + shape, data = d)
  # cluster-robust SE by scenario cell (sandwich, no extra deps)
  cells <- factor(d$scenario_id)
  X <- model.matrix(fit); e <- residuals(fit)
  meat <- t(X) %*% (e^2 * X)
  k <- max(2, nlevels(cells))
  G <- lapply(levels(cells), function(g) X[cells == g, , drop = FALSE] %*%
                e[cells == g])
  M <- Reduce(`+`, lapply(G, function(g) g %*% t(g)))
  bread <- solve(t(X) %*% X)
  V <- bread %*% M %*% bread * (k / (k - 1)) *
    ((nrow(d) - 1) / (nrow(d) - ncol(X)))
  se <- sqrt(diag(V))
  est <- coef(fit)
  data.frame(term = names(est), estimate = est, se = se,
             t = est / se,
             p_one_sided_negative = stats::pt(est / se, df = nrow(d) - ncol(X)))
}

# Paired McNemar test (reserved for Phase B comparator arms; NOT used in any
# Phase A readout - SAP B1 removed the CMH cross-cell test: truth is a
# cell-level constant, so within-stratum 2x2 tables degenerate).
mcnemar_p <- function(det_int, det_cmp) {
  ok <- !is.na(det_int) & !is.na(det_cmp)
  if (!any(ok)) return(NA_real_)
  stats::mcnemar.test(table(det_int[ok], det_cmp[ok]), correct = FALSE)$p.value
}
