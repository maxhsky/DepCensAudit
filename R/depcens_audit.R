#' Run the integrated DepCens-Audit pipeline
#'
#' Chains the three dimensions (evaluation, positivity, extrapolation) and
#' assembles the two-layer flag report (root-cause + performance layers).
#'
#' \strong{Phase A scope}: the dimensions run independently and their outputs
#' are concatenated. The cross-dimensional interaction pathway
#' (\code{.interaction}) is a \emph{frozen interface} reserved for the Phase B
#' methodology release: passing a non-\code{NULL} value is an error in this
#' version, and no cross-dimensional performance claim is made anywhere in
#' this package.
#'
#' @param data A data frame.
#' @param time,event Bare column names (1 = event, 0 = censored).
#' @param horizon Decision-analytic horizon for the extrapolation audit.
#' @param dimensions Character vector subset of
#'   \code{c("evaluation","positivity","extrapolation")}.
#' @param pred_fun Prediction function for dimension A (see
#'   \code{\link{audit_evaluation}}); required only if \code{"evaluation"} is
#'   requested.
#' @param config A \code{\link{psa_config}} object.
#' @param .interaction Frozen Phase B interface. Must be \code{NULL}.
#' @param seed Optional seed stored in report metadata.
#' @param ... Passed to the dimension auditors.
#'
#' @return An object of class \code{depcens_audit}: a list with one element
#'   per requested dimension, plus \code{flags} (combined two-layer report),
#'   \code{psa} (strategy mapped from the worst performance flag), and
#'   \code{meta} (seed, config version, session info, call).
#' @export
depcens_audit <- function(data, time, event, horizon,
                          dimensions = c("evaluation", "positivity",
                                         "extrapolation"),
                          pred_fun = NULL,
                          config = psa_config(),
                          .interaction = NULL,
                          seed = NULL,
                          ...) {
  .check_interaction_reserved(.interaction)
  dimensions <- match.arg(dimensions, several.ok = TRUE)
  time_sym <- as.name(deparse(substitute(time)))
  event_sym <- as.name(deparse(substitute(event)))

  res <- list()
  if ("evaluation" %in% dimensions) {
    if (is.null(pred_fun))
      stop("`pred_fun` is required when dimension 'evaluation' is requested.")
    res$evaluation <- do.call(
      audit_evaluation,
      c(list(data = data, time = time_sym, event = event_sym,
             pred_fun = pred_fun, .interaction = NULL), list(...)))
  }
  if ("positivity" %in% dimensions)
    res$positivity <- do.call(
      audit_positivity,
      c(list(data = data, time = time_sym, event = event_sym,
             .interaction = NULL), list(...)))
  if ("extrapolation" %in% dimensions)
    res$extrapolation <- do.call(
      audit_extrapolation,
      c(list(data = data, time = time_sym, event = event_sym,
             horizon = horizon, config = config, .interaction = NULL,
             seed = seed), list(...)))

  all_flags <- unlist(lapply(res, function(r) r$flags$root),
                      recursive = FALSE)
  perf_flags <- unlist(lapply(res, function(r) r$flags$performance),
                       recursive = FALSE)
  worst <- combine_flags(perf_flags)
  psa <- if ("extrapolation" %in% dimensions)
    res$extrapolation$psa else psa_mapping(worst$level, config)

  out <- list(
    dimensions = res,
    flags = list(root = all_flags, performance = perf_flags,
                 overall = worst),
    psa = psa,
    meta = list(
      phase = "A (engineering integration; no cross-dimensional claims)",
      seed = seed,
      config_version = config$version,
      session_info = utils::sessionInfo(),
      call = match.call(),
      audit_time = Sys.time()
    )
  )
  class(out) <- "depcens_audit"
  out
}

#' Test for DepCens-Audit objects
#'
#' @param x Object to test.
#' @return TRUE if \code{x} inherits from class \code{depcens_audit}.
#' @export
is.depcens_audit <- function(x) inherits(x, "depcens_audit")

#' @export
print.depcens_audit <- function(x, ...) {
  cat("DepCens-Audit :: integrated audit (Phase A)\n")
  cat("  dimensions run:", paste(names(x$dimensions), collapse = ", "), "\n")
  cat("  overall performance flag: [", toupper(x$flags$overall$level), "] ",
      x$flags$overall$code, "\n", sep = "")
  if (length(x$flags$root))
    cat(paste0(format_flags(x$flags$root), collapse = "\n"), "\n")
  cat("  PSA strategy:", x$psa$strategy, "\n")
  cat("  phase note:", x$meta$phase, "\n")
  invisible(x)
}

#' @export
summary.depcens_audit <- function(object, ...) {
  data.frame(
    dimension = names(object$dimensions),
    n_root_flags = vapply(object$dimensions,
                          function(d) length(d$flags$root), integer(1)),
    worst_flag = vapply(object$dimensions, function(d) {
      # Performance layer if present, else fall back to root-cause layer so
      # diagnostic root flags (e.g. sparse tail risk sets) stay visible.
      pf <- d$flags$performance
      src <- if (length(pf)) pf else d$flags$root
      if (!length(src)) return("green")
      lvls <- vapply(src, function(f) match(f$level, .flag_levels),
                     integer(1))
      .flag_levels[max(lvls)]
    }, character(1)),
    row.names = NULL
  )
}
