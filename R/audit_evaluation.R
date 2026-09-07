#' Evaluation audit (dimension A)
#'
#' Model-agnostic IPCW-based Brier / integrated Brier score wrapper, reusing
#' the IBS-Dep idea (DependentEVAL, arXiv:2502.19460; UAI 2026 per the
#' arXiv comment field) without re-implementing its internals. If the
#' DependentEVAL package is installed it can be used by the caller upstream;
#' this function only requires a prediction function.
#'
#' \strong{Status: experimental (work item A-1).} The operating
#' characteristics of IBS-style metrics used \emph{as audit signals} are under
#' validation in the M3 simulation study. IPCW assumes a correctly specified
#' censoring model with bounded weights and is fragile precisely in the
#' dependent-censoring scenarios this package targets; competing risks are out
#' of scope for v0.0.x. Use the output as a screening signal, not as proof of
#' model adequacy.
#'
#' @param data A data frame.
#' @param time,event Bare column names (1 = event, 0 = censored).
#' @param pred_fun Function \code{function(newdata, t0)} returning, for each
#'   row of \code{newdata}, the predicted survival probability at time
#'   \code{t0}. Any model class is acceptable (Cox, ML, etc.).
#' @param eval_times Numeric vector of time points at which the Brier score is
#'   evaluated; IBS is the trapezoidal integral over this grid. Defaults to 25
#'   quantiles of observed event times.
#' @param .interaction Frozen Phase B interface. Must be \code{NULL}.
#'
#' @return An object of class \code{depcens_evaluation_audit} with
#'   \code{brier} (per-time-point IPCW Brier scores), \code{ibs} (integrated
#'   Brier score), and \code{flags} (currently limited to weight-instability
#'   diagnostics).
#' @export
audit_evaluation <- function(data, time, event, pred_fun,
                             eval_times = NULL,
                             .interaction = NULL) {
  .check_interaction_reserved(.interaction)
  stopifnot(is.function(pred_fun))
  time_v <- data[[deparse(substitute(time))]]
  event_v <- data[[deparse(substitute(event))]]
  stopifnot(is.numeric(time_v), all(event_v %in% c(0, 1)))

  if (is.null(eval_times)) {
    ev_times <- time_v[event_v == 1]
    eval_times <- as.numeric(stats::quantile(
      ev_times, probs = seq(0.05, 0.95, length.out = 25), names = FALSE))
  }

  # Censoring KM for IPCW weights (Graf et al. 1999 scheme).
  cens_fit <- survival::survfit(survival::Surv(time_v, 1 - event_v) ~ 1)
  G <- function(t) {
    s <- summary(cens_fit, times = pmax(t, 1e-12), extend = TRUE)$surv
    pmax(as.numeric(s), 1e-6)
  }

  brier <- numeric(length(eval_times))
  max_weight <- 0
  for (i in seq_along(eval_times)) {
    t0 <- eval_times[i]
    s_hat <- as.numeric(pred_fun(data, t0))
    stopifnot(length(s_hat) == length(time_v))
    y <- as.numeric(time_v > t0)
    w <- ifelse(time_v <= t0 & event_v == 1, 1 / G(time_v),
                ifelse(time_v > t0, 1 / G(t0), 0))
    max_weight <- max(max_weight, w)
    brier[i] <- mean(w * (y - s_hat)^2)
  }
  ibs <- sum(diff(eval_times) *
               (head(brier, -1) + tail(brier, -1)) / 2) /
    (max(eval_times) - min(eval_times))

  flags <- list()
  if (max_weight > 10)
    flags <- c(flags, list(new_flag(
      "yellow", "evaluation", "root",
      code = "IPCW_WEIGHT_INSTABILITY",
      message = sprintf(
        paste0("Maximum IPCW weight is %.1f (>10); the censoring model ",
               "yields near-positivity-violating weights and IBS is unstable."),
        max_weight),
      evidence = list(max_weight = max_weight))))

  out <- list(
    brier = data.frame(time = eval_times, brier = brier),
    ibs = ibs,
    flags = list(root = flags, performance = list()),
    meta = list(
      status = paste0(
        "experimental (A-1): operating characteristics as an audit signal ",
        "under validation; IPCW assumes correct censoring model, bounded ",
        "weights, no competing risks."),
      session_info = utils::sessionInfo(),
      call = match.call()
    )
  )
  class(out) <- "depcens_evaluation_audit"
  out
}

#' @export
print.depcens_evaluation_audit <- function(x, ...) {
  cat("DepCens-Audit :: evaluation audit (dimension A)\n")
  cat("  IBS (IPCW):", signif(x$ibs, 4), "\n")
  if (length(x$flags$root))
    cat(paste0(format_flags(x$flags$root), collapse = "\n"), "\n")
  cat("  status:", x$meta$status, "\n")
  invisible(x)
}
