#' PSA configuration object (versioned thresholds)
#'
#' Creates the versioned configuration governing how extrapolation-audit
#' results map to probabilistic sensitivity analysis (PSA) strategies.
#' All thresholds are initial defaults and must be calibrated against the
#' M3 simulation study before regulatory-facing use; the \code{version}
#' field is embedded in every audit report for traceability.
#'
#' @param rho_yellow Numeric. Instability criterion: relative RMST difference
#'   between full-data and truncated-refit fits at which a yellow flag is
#'   raised. Default 0.10.
#' @param rho_maturity_yellow Numeric. Maturity criterion
#'   \eqn{\rho = (E[S] - RMST_{cutoff}) / E[S]} yellow threshold. Default 0.20.
#' @param rho_maturity_red Numeric. Maturity red threshold. Default 0.50.
#' @param disagreement_escalate Numeric. If the relative spread of the three
#'   candidate models' decision-window RMST exceeds this value, the flag is
#'   escalated by at least one level. Default 0.15.
#' @param mixture_weights Character. Weights for the structural-mixture PSA
#'   used under a yellow flag. \code{"akaike"} (default) uses Akaike weights.
#'   \code{"audit"} (SSE-O audit-weighted) is available as an option but is
#'   \emph{not} the default because its validity is unverified by simulation.
#' @param version Character. Threshold-set version tag.
#'
#' @return An object of class \code{psa_config}.
#' @export
psa_config <- function(rho_yellow = 0.10,
                       rho_maturity_yellow = 0.20,
                       rho_maturity_red = 0.50,
                       disagreement_escalate = 0.15,
                       mixture_weights = c("akaike", "audit"),
                       version = "0.1.0") {
  mixture_weights <- match.arg(mixture_weights)
  stopifnot(
    rho_yellow > 0, rho_maturity_yellow > 0,
    rho_maturity_red > rho_maturity_yellow,
    disagreement_escalate > 0
  )
  structure(
    list(
      rho_yellow = rho_yellow,
      rho_maturity_yellow = rho_maturity_yellow,
      rho_maturity_red = rho_maturity_red,
      disagreement_escalate = disagreement_escalate,
      mixture_weights = mixture_weights,
      version = version
    ),
    class = "psa_config"
  )
}

#' @export
print.psa_config <- function(x, ...) {
  cat("PSA configuration (version", x$version, ")\n")
  cat("  rho instability (yellow):     ", x$rho_yellow, "\n")
  cat("  rho maturity (yellow / red):  ", x$rho_maturity_yellow, "/",
      x$rho_maturity_red, "\n")
  cat("  model disagreement escalate:  ", x$disagreement_escalate, "\n")
  cat("  mixture weights:              ", x$mixture_weights, "\n")
  invisible(x)
}

#' Map an extrapolation flag to a PSA strategy
#'
#' Implements the calibrated mapping rules (HTA advisor input, 2026-09-07):
#'
#' \itemize{
#'   \item \strong{green}: standard TSD14 workflow; single best curve with
#'     asymptotic-covariance parameter sampling; no interval widening.
#'   \item \strong{yellow}: structural-mixture PSA (sample a curve by weight,
#'     then sample its parameters), \emph{or} doubling of parameter
#'     uncertainty; base case remains the single best curve, mixture enters
#'     scenario analysis; TSD26 expert-elicitation priors may be mounted in
#'     the sensitivity layer only (never re-fit the base case).
#'   \item \strong{red}: do \emph{not} widen intervals. Output each candidate
#'     curve's full cost-effectiveness results, stratified CEACs, and an
#'     ICER tipping-point analysis; explicitly state that data maturity is
#'     insufficient to support a stable base case; consider evidence
#'     generation / data-cut updates / managed entry agreements.
#' }
#'
#' @param flag Character, one of \code{"green"}, \code{"yellow"}, \code{"red"}.
#' @param config A \code{\link{psa_config}} object.
#'
#' @return A list describing the PSA strategy: \code{strategy} (tag),
#'   \code{base_case}, \code{sampling}, \code{sensitivity_layer},
#'   \code{reporting} (character vector of required outputs), and
#'   \code{diagnostic_wording} (mandated phrasing).
#' @export
psa_mapping <- function(flag = c("green", "yellow", "red"),
                        config = psa_config()) {
  flag <- match.arg(flag)
  switch(
    flag,
    green = list(
      strategy = "standard_tsd14",
      base_case = "single best curve (pre-specified; data-driven selection prohibited)",
      sampling = "asymptotic covariance parameter sampling; no interval widening",
      sensitivity_layer = "none required by audit",
      reporting = c("single curve PSA", "CEAC"),
      diagnostic_wording = "Extrapolation behaviour is adequately constrained by the observed data."
    ),
    yellow = list(
      strategy = if (config$mixture_weights == "akaike")
        "structural_mixture_psa" else "structural_mixture_psa_audit_weighted",
      base_case = "single best curve retained as base case",
      sampling = paste0(
        "two-stage: sample curve by ", config$mixture_weights,
        " weight, then sample curve parameters; alternative: inflate parameter uncertainty x2"
      ),
      sensitivity_layer = "TSD26 expert-elicitation priors mounted in sensitivity layer only",
      reporting = c("mixture PSA (scenario analysis)", "CEAC per scenario",
                    "flag criteria and thresholds (with config version)"),
      diagnostic_wording = "Tail behaviour is partially constrained; structural uncertainty should be propagated via mixture PSA in scenario analyses."
    ),
    red = list(
      strategy = "dual_scenario_no_widening",
      base_case = "no stable base case; report each candidate scenario (e.g., cure-like vs standard) in full",
      sampling = "per-curve parameter sampling within each scenario; interval widening is prohibited",
      sensitivity_layer = "expert priors optional in sensitivity layer",
      reporting = c("full cost-effectiveness results per candidate curve",
                    "stratified CEACs", "ICER tipping-point analysis",
                    "explicit statement of insufficient data maturity",
                    "options: evidence generation / data-cut update / managed entry"),
      diagnostic_wording = "Data support is insufficient to constrain tail behaviour; decision uncertainty should be represented as discrete scenarios rather than widened intervals."
    )
  )
}
