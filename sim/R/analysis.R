# ------------------------------------------------------------------------------
# M3a simulation framework :: analysis helpers (SAP estimands)
#
# DRAFT - the concrete detection criteria and tests are FROZEN in the OSF
# pre-registration before any performance readout (project doc v1.0 r2, S5).
# These helpers implement the estimands; nothing here is final until frozen.
# ------------------------------------------------------------------------------

library(survival)

# Dataset-level detection rule (DRAFT - to be calibrated on pilot runs and
# frozen in the SAP before the confirmatory batch):
#   B-defect detected  : positivity sensitivity-interval width >= w_B
#   C-defect detected  : extrapolation performance flag in {yellow, red}
#                        OR any C-type root flag code
#   A-defect detected  : evaluation weight-instability flag raised
detect <- function(rows, w_B = NULL) {
  if (is.null(w_B)) {
    # uncalibrated placeholder: 90th percentile of kappa = 0 cells in the
    # SAME batch (cell-matched specificity anchor; SAP freezes the rule)
    w_B <- stats::quantile(rows$pos_interval_width[rows$kappa == 0 &
                                                     !is.na(rows$pos_interval_width)],
                           probs = 0.90, na.rm = TRUE)
  }
  data.frame(
    scenario_id = rows$scenario_id, rep = rows$rep,
    truth_B = rows$truth_B, truth_C = rows$truth_C, truth_A = rows$truth_A,
    det_B = rows$pos_interval_width >= w_B,
    det_C = rows$ext_perf_level %in% c("yellow", "red") |
      grepl("SPARSE|IMMATURE|INSTAB", rows$ext_root_codes),
    det_A = rows$eval_weight_unstable == 1)
}

# Per-scenario sensitivity / specificity for one dimension.
sens_spec <- function(d, dimension = c("B", "C", "A")) {
  dimension <- match.arg(dimension)
  det <- switch(dimension, B = d$det_B, C = d$det_C, A = d$det_A)
  tru <- switch(dimension, B = d$truth_B, C = d$truth_C, A = d$truth_A)
  tp <- sum(det & tru, na.rm = TRUE);  fn <- sum(!det & tru, na.rm = TRUE)
  fp <- sum(det & !tru, na.rm = TRUE); tn <- sum(!det & !tru, na.rm = TRUE)
  data.frame(dimension = dimension,
             sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA,
             specificity = if (fp + tn > 0) tn / (fp + tn) else NA,
             tp = tp, fp = fp, fn = fn, tn = tn)
}

# Paired McNemar test within one scenario (integrated vs single-dimension
# comparator arm; comparator wiring lands with the M3a comparator layer).
mcnemar_p <- function(det_int, det_cmp) {
  ok <- !is.na(det_int) & !is.na(det_cmp)
  if (!any(ok)) return(NA_real_)
  stats::mcnemar.test(table(det_int[ok], det_cmp[ok]), correct = FALSE)$p.value
}

# Stratified (across-scenario) CMH test of detection superiority.
cmh_p <- function(det, tru, strata) {
  ok <- !is.na(det) & !is.na(tru)
  tbl <- array(0, dim = c(2, 2, nlevels(factor(strata))))
  for (s in levels(factor(strata))) {
    idx <- ok & strata == s
    tbl[, , s] <- table(det[idx], tru[idx])
  }
  tryCatch(
    stats::mantelhaen.test(tbl, correct = FALSE)$p.value,
    error = function(e) NA_real_)
}
