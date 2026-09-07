#' Positivity audit (dimension B)
#'
#' By design this dimension does \emph{not} emit a green/yellow/red verdict on
#' the censoring mechanism: censoring independence is not directly testable
#' (preprint evidence E4). Instead it reports:
#'
#' \enumerate{
#'   \item a \strong{continuous sensitivity interval}: how a target survival
#'     estimate moves when censoring weights are tilted by a sensitivity
#'     parameter \eqn{\kappa} over \code{kappa_grid}. Tilting is a
#'     sensitivity-analysis device only; it is not an identification strategy,
#'     and copula-family misspecification is out of scope;
#'   \item \strong{detectable positivity sub-problems} (the only place flags
#'     are allowed): sparse tail risk sets.
#' }
#'
#' @param data A data frame.
#' @param time,event Bare column names (1 = event, 0 = censored).
#' @param kappa_grid Numeric vector of tilting strengths to scan. Default
#'   \code{seq(0, 0.6, by = 0.2)}.
#' @param tilting Character. \code{"risk"} tilts censoring weights by a Cox
#'   risk score when covariates are supplied via \code{risk}; otherwise falls
#'   back to \code{"time"} (tilt by follow-up time rank).
#' @param risk Optional numeric vector (same length as data rows) used as the
#'   tilting score, e.g. a prognostic risk score.
#' @param tau Time point at which the sensitivity interval of S(t) is
#'   reported. Defaults to the median observed event time.
#' @param tail_window Proportion of follow-up defining the tail window for the
#'   sparse-risk-set check. Default 0.1.
#' @param min_tail_atrisk Minimum proportion of subjects expected at risk in
#'   the tail window before a flag is raised. Default 0.05.
#' @param .interaction Frozen Phase B interface. Must be \code{NULL}.
#'
#' @return An object of class \code{depcens_positivity_audit} with
#'   \code{sensitivity} (per-kappa weighted S(tau) and the resulting interval)
#'   and \code{flags} (restricted to detectable positivity sub-problems).
#' @export
audit_positivity <- function(data, time, event,
                             kappa_grid = seq(0, 0.6, by = 0.2),
                             tilting = c("risk", "time"),
                             risk = NULL,
                             tau = NULL,
                             tail_window = 0.1,
                             min_tail_atrisk = 0.05,
                             .interaction = NULL) {
  .check_interaction_reserved(.interaction)
  tilting <- match.arg(tilting)
  time_v <- data[[deparse(substitute(time))]]
  event_v <- data[[deparse(substitute(event))]]
  stopifnot(is.numeric(time_v), all(event_v %in% c(0, 1)))
  kappa_grid <- sort(unique(c(0, kappa_grid)))

  if (is.null(tau))
    tau <- stats::median(time_v[event_v == 1], na.rm = TRUE)

  score <- if (tilting == "risk" && !is.null(risk)) {
    stopifnot(length(risk) == length(time_v))
    as.numeric(scale(risk))
  } else {
    as.numeric(scale(rank(time_v, ties.method = "average")))
  }

  # Weighted Kaplan-Meier S(tau) under tilting: censored subjects receive
  # weight exp(kappa * score); events keep weight 1. Sensitivity device only.
  wkm_surv <- function(kappa) {
    w <- ifelse(event_v == 0, exp(kappa * score), 1)
    ord <- order(time_v)
    tt <- time_v[ord]; ee <- event_v[ord]; ww <- w[ord]
    d_times <- sort(unique(tt[ee == 1 & tt <= tau]))
    s <- 1
    for (tj in d_times) {
      at_risk <- sum(ww[tt >= tj])
      d_j <- sum(ww[tt == tj & ee == 1])
      if (at_risk > 0) s <- s * (1 - d_j / at_risk)
    }
    s
  }

  sens <- data.frame(
    kappa = kappa_grid,
    surv_tau = vapply(kappa_grid, wkm_surv, numeric(1))
  )
  interval <- range(sens$surv_tau, na.rm = TRUE)

  flags <- list()
  t_max <- max(time_v, na.rm = TRUE)
  tail_start <- t_max * (1 - tail_window)
  at_risk_tail <- mean(time_v >= tail_start)
  if (at_risk_tail < min_tail_atrisk)
    flags <- c(flags, list(new_flag(
      "yellow", "positivity", "root",
      code = "SPARSE_TAIL_RISKSET",
      message = sprintf(
        paste0("Only %.1f%% of subjects remain at risk in the final %.0f%% ",
               "of follow-up; tail estimates are data-sparse."),
        100 * at_risk_tail, 100 * tail_window),
      evidence = list(at_risk_tail = at_risk_tail, tail_start = tail_start))))

  out <- list(
    sensitivity = sens,
    interval = interval,
    tau = tau,
    tilting = tilting,
    flags = list(root = flags, performance = list()),
    meta = list(
      kappa_grid = kappa_grid,
      note = paste0(
        "Continuous sensitivity interval; NOT an identification result. ",
        "Censoring independence is untestable; no traffic-light verdict is ",
        "issued on the mechanism itself."),
      session_info = utils::sessionInfo(),
      call = match.call()
    )
  )
  class(out) <- "depcens_positivity_audit"
  out
}

#' @export
print.depcens_positivity_audit <- function(x, ...) {
  cat("DepCens-Audit :: positivity audit (dimension B)\n")
  cat("  sensitivity interval for S(tau = ", signif(x$tau, 3), "): [",
      signif(x$interval[1], 3), ", ", signif(x$interval[2], 3), "]\n", sep = "")
  cat("  kappa grid:", paste(x$meta$kappa_grid, collapse = ", "), "\n")
  if (length(x$flags$root))
    cat(paste0(format_flags(x$flags$root), collapse = "\n"), "\n")
  cat("  note:", x$meta$note, "\n")
  invisible(x)
}
