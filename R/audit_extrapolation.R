# ---- internal helpers ---------------------------------------------------

# Survival function of a fitted candidate model at times t.
.surv_of <- function(fit, t) {
  if (inherits(fit, "depcens_constail")) {
    base <- .surv_of(fit$fit, pmin(t, fit$t_max))
    ext <- exp(-fit$h_tail * pmax(0, t - fit$t_max))
    return(base * ext)
  }
  s <- summary(fit, t = t, type = "survival")[[1]]$est
  as.numeric(s)
}

# Hazard of a flexsurv fit at a single time point.
.hazard_of <- function(fit, t) {
  as.numeric(summary(fit, t = t, type = "hazard")[[1]]$est)
}

# RMST to `horizon` for any candidate model (numeric integration).
.rmst_of <- function(fit, horizon) {
  f <- function(t) .surv_of(fit, t)
  out <- tryCatch(
    stats::integrate(Vectorize(f), lower = 0, upper = horizon,
                     subdivisions = 200L)$value,
    error = function(e) NA_real_
  )
  out
}

# Fit one candidate model. Returns a fit object (flexsurvreg or depcens_constail).
.fit_candidate <- function(model, time, event) {
  df <- data.frame(time = time, event = event)
  fml <- survival::Surv(time, event) ~ 1
  switch(
    model,
    rp2 = flexsurv::flexsurvspline(fml, data = df, k = 2, scale = "hazard"),
    psm = tryCatch(
      flexsurv::flexsurvreg(fml, data = df, dist = "gengamma"),
      error = function(e) flexsurv::flexsurvreg(fml, data = df, dist = "weibull")
    ),
    constail = {
      # Simplified constant-hazard tail: RP(2) up to the last observed event
      # time, constant hazard beyond it. Documented simplification for v0.0.1;
      # the NMA-grade constant-tail-HR variant is out of scope for single-
      # cohort IPD (evidence E1 is NMA-context).
      base <- flexsurv::flexsurvspline(fml, data = df, k = 2, scale = "hazard")
      t_max <- max(time[event == 1], na.rm = TRUE)
      structure(list(fit = base, t_max = t_max,
                     h_tail = .hazard_of(base, t_max)),
                class = "depcens_constail")
    },
    stop("unknown model tag: ", model)
  )
}

.aic_of <- function(fit) {
  if (inherits(fit, "depcens_constail")) return(stats::AIC(fit$fit))
  stats::AIC(fit)
}

# Administratively truncate the data at quantile `q` of observed times.
.truncate_data <- function(time, event, q) {
  t_star <- as.numeric(stats::quantile(time, probs = q, type = 1))
  list(time = pmin(time, t_star),
       event = ifelse(time <= t_star, event, 0L),
       t_star = t_star)
}

# Lightweight reproducibility fingerprint (no external digest dependency).
.data_fingerprint <- function(time, event) {
  s <- paste0(format(time, digits = 12), collapse = ",")
  s2 <- paste0(event, collapse = "")
  h <- (sum(utf8ToInt(s)) * 1315423911 + sum(utf8ToInt(s2)) * 2654435761) %% 2^31
  sprintf("fp%08x", as.integer(h))
}

.check_interaction_reserved <- function(.interaction) {
  if (!is.null(.interaction))
    stop(
      "`.interaction` is a frozen interface reserved for the Phase B ",
      "methodology release and is not implemented in this version. ",
      "Pass `.interaction = NULL`.",
      call. = FALSE
    )
  invisible(TRUE)
}

# ---- main function --------------------------------------------------------

#' Extrapolation audit (dimension C)
#'
#' Fits three candidate extrapolation models (Royston-Parmar 2-knot spline,
#' a standard parametric model, and a constant-tail variant), re-fits them on
#' administratively truncated data, and computes:
#'
#' \itemize{
#'   \item \strong{rho instability}: max over models of the relative RMST
#'     difference (at \code{horizon}) between full-data and truncated-refit
#'     fits. Yellow at \code{>= config$rho_yellow} (default 0.10).
#'   \item \strong{rho maturity}: \eqn{(E[S] - RMST_{cutoff}) / E[S]} from the
#'     best-AIC model, the share of expected survival contributed by the
#'     unobserved tail. Yellow at 0.20, red at 0.50 (defaults).
#'   \item \strong{model disagreement}: relative spread of decision-window
#'     RMST across the three candidates; exceeding
#'     \code{config$disagreement_escalate} escalates the flag one level.
#' }
#'
#' The resulting flag maps to a PSA strategy via \code{\link{psa_mapping}}.
#' Output is diagnostic only; it does not judge regulatory compliance.
#'
#' @param data A data frame.
#' @param time,event Bare column names: follow-up time and event indicator
#'   (1 = event, 0 = censored).
#' @param horizon Numeric. Decision-analytic time horizon for RMST.
#' @param models Character vector subset of \code{c("rp2","psm","constail")}.
#' @param trunc_quantile Quantile of observed times at which the refit data
#'   are administratively truncated. Default 0.8.
#' @param config A \code{\link{psa_config}} object.
#' @param .interaction Frozen Phase B interface. Must be \code{NULL}.
#' @param seed Optional seed stored in the report metadata.
#'
#' @return An object of class \code{depcens_extrapolation_audit} with elements
#'   \code{fits}, \code{fits_trunc}, \code{metrics} (per-model and summary),
#'   \code{flags} (root-cause and performance layers), \code{psa} (strategy
#'   from \code{\link{psa_mapping}}), and \code{meta} (seed, session info,
#'   data fingerprint, config version, call).
#' @export
audit_extrapolation <- function(data, time, event, horizon,
                                models = c("rp2", "psm", "constail"),
                                trunc_quantile = 0.8,
                                config = psa_config(),
                                .interaction = NULL,
                                seed = NULL) {
  .check_interaction_reserved(.interaction)
  if (!requireNamespace("flexsurv", quietly = TRUE))
    stop("Package 'flexsurv' is required for audit_extrapolation().")
  if (!is.null(seed)) set.seed(seed)

  time_v <- data[[deparse(substitute(time))]]
  event_v <- data[[deparse(substitute(event))]]
  stopifnot(is.numeric(time_v), all(event_v %in% c(0, 1)), horizon > 0)
  models <- match.arg(models, several.ok = TRUE)

  fits <- stats::setNames(lapply(models, .fit_candidate,
                                 time = time_v, event = event_v), models)

  tr <- .truncate_data(time_v, event_v, trunc_quantile)
  fits_trunc <- stats::setNames(
    lapply(models, .fit_candidate, time = tr$time, event = tr$event), models)

  rmst_full <- vapply(fits, .rmst_of, numeric(1), horizon = horizon)
  rmst_trunc <- vapply(fits_trunc, .rmst_of, numeric(1), horizon = horizon)
  rho_inst_per_model <- abs(rmst_full - rmst_trunc) / rmst_full
  rho_instability <- max(rho_inst_per_model, na.rm = TRUE)

  # Maturity: share of expected survival beyond the data cutoff (best model).
  best <- names(which.min(vapply(fits, .aic_of, numeric(1))))
  t_cut <- max(time_v, na.rm = TRUE)
  rmst_cut <- .rmst_of(fits[[best]], t_cut)
  e_surv <- .rmst_of(fits[[best]], max(horizon, 5 * t_cut))
  rho_maturity <- if (is.finite(e_surv) && e_surv > 0)
    max(0, (e_surv - rmst_cut) / e_surv) else NA_real_

  disagreement <- if (length(models) > 1) {
    rng <- range(rmst_full, na.rm = TRUE)
    (rng[2] - rng[1]) / mean(rmst_full, na.rm = TRUE)
  } else 0

  # ---- root-cause layer flags ----
  flags_root <- list()
  if (is.finite(rho_instability) && rho_instability >= config$rho_yellow)
    flags_root <- c(flags_root, list(new_flag(
      "yellow", "extrapolation", "root",
      code = "EXTRAP_INSTABILITY",
      message = sprintf(
        paste0("Truncated-refit RMST differs by %.1f%% (threshold %.0f%%); ",
               "tail behaviour is sensitive to the last portion of follow-up."),
        100 * rho_instability, 100 * config$rho_yellow),
      evidence = list(rho_instability = rho_instability,
                      per_model = rho_inst_per_model, t_star = tr$t_star))))
  if (is.finite(rho_maturity)) {
    if (rho_maturity >= config$rho_maturity_red)
      flags_root <- c(flags_root, list(new_flag(
        "red", "extrapolation", "root",
        code = "IMMATURE_DATA",
        message = sprintf(
          paste0("%.0f%% of expected survival lies beyond the data cutoff; ",
                 "data support is insufficient to constrain tail behaviour."),
          100 * rho_maturity),
        evidence = list(rho_maturity = rho_maturity, best_model = best))))
    else if (rho_maturity >= config$rho_maturity_yellow)
      flags_root <- c(flags_root, list(new_flag(
        "yellow", "extrapolation", "root",
        code = "PARTIAL_MATURITY",
        message = sprintf(
          "%.0f%% of expected survival lies beyond the data cutoff.",
          100 * rho_maturity),
        evidence = list(rho_maturity = rho_maturity, best_model = best))))
  }

  # ---- performance layer (conditional) ----
  perf <- combine_flags(flags_root)
  if (is.finite(disagreement) && disagreement > config$disagreement_escalate &&
      perf$level != "green") {
    perf <- escalate_flag(perf, 1L, reason = sprintf(
      "model disagreement %.1f%% > %.0f%%",
      100 * disagreement, 100 * config$disagreement_escalate))
    perf$triggered_by <- c(perf$triggered_by, "MODEL_DISAGREEMENT")
  }
  if (perf$level == "green" &&
      is.finite(disagreement) && disagreement > config$disagreement_escalate)
    perf <- new_flag(
      "yellow", "extrapolation", "performance",
      code = "MODEL_DISAGREEMENT",
      message = sprintf(
        paste0("Candidate models' decision-window RMST spread is %.1f%% ",
               "(> %.0f%%); structural uncertainty is non-trivial."),
        100 * disagreement, 100 * config$disagreement_escalate),
      evidence = list(disagreement = disagreement, rmst_full = rmst_full))

  psa <- psa_mapping(perf$level, config)

  out <- list(
    fits = fits, fits_trunc = fits_trunc,
    metrics = list(
      rmst_full = rmst_full, rmst_trunc = rmst_trunc,
      rho_instability_per_model = rho_inst_per_model,
      rho_instability = rho_instability,
      rho_maturity = rho_maturity,
      model_disagreement = disagreement,
      best_model_by_aic = best,
      t_star = tr$t_star, horizon = horizon
    ),
    flags = list(root = flags_root, performance = list(perf)),
    psa = psa,
    meta = list(
      seed = seed,
      config_version = config$version,
      config = config,
      data_fingerprint = .data_fingerprint(time_v, event_v),
      session_info = utils::sessionInfo(),
      call = match.call(),
      audit_time = Sys.time()
    )
  )
  class(out) <- "depcens_extrapolation_audit"
  out
}

#' @export
print.depcens_extrapolation_audit <- function(x, ...) {
  cat("DepCens-Audit :: extrapolation audit (dimension C)\n")
  cat("  horizon:", x$metrics$horizon,
      " | truncation refit at t* =", round(x$metrics$t_star, 2), "\n")
  cat("  rho instability:", signif(x$metrics$rho_instability, 3),
      " | rho maturity:", signif(x$metrics$rho_maturity, 3),
      " | disagreement:", signif(x$metrics$model_disagreement, 3), "\n")
  cat("  performance flag: [", toupper(x$flags$performance[[1]]$level), "] ",
      x$flags$performance[[1]]$code, "\n", sep = "")
  cat("  PSA strategy:", x$psa$strategy, "\n")
  cat("  diagnostic:", x$psa$diagnostic_wording, "\n")
  invisible(x)
}
