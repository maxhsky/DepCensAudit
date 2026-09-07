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
    det_B = !err & nzchar(rows$pos_flag_codes),
    det_A = !err & rows$eval_weight_unstable == 1)
}

# Per-scenario sensitivity / specificity for one dimension. Errors stay in
# the denominator as non-detections (SAP B3).
sens_spec <- function(d, dimension = c("B", "C", "A")) {
  dimension <- match.arg(dimension)
  det <- switch(dimension, B = d$det_B, C = d$det_C, A = d$det_A)
  tru <- switch(dimension, B = d$truth_B, C = d$truth_C, A = d$truth_A)
  tp <- sum(det & tru, na.rm = TRUE);  fn <- sum(!det & tru, na.rm = TRUE)
  fp <- sum(det & !tru, na.rm = TRUE); tn <- sum(!det & !tru, na.rm = TRUE)
  data.frame(dimension = dimension,
             sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA,
             specificity = if (fp + tn > 0) tn / (fp + tn) else NA,
             error_rate = mean(d$errored),
             tp = tp, fp = fp, fn = fn, tn = tn)
}

# E2 (false-positive control): clean scene = truth_C none AND kappa == 0.
# Criterion: pooled non-green rate over ANY dimension (root or performance
# layer) <= 0.10, judged by an exact one-sided binomial test; execution
# thresholds (SAP A1): pooled n = 2000 -> observed rate <= 0.0866;
# per-cell n = 500 -> <= 0.072.
e2_clean_scene_rate <- function(rows) {
  clean <- !is_err(rows) & !rows$truth_C & rows$kappa == 0
  ng <- !is.na(rows$overall_level) & rows$overall_level != "green"
  n <- sum(clean); x <- sum(clean & ng)
  p_hat <- x / n
  pval <- tryCatch(
    stats::binom.test(x, n, p = 0.10, alternative = "less")$p.value,
    error = function(e) NA_real_)
  data.frame(n_clean = n, n_flagged = x, rate = p_hat,
             exact_p_less_0.10 = pval,
             pass_pooled_0866 = p_hat <= 0.0866)
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
