# Two-layer flag system (methodology review, 2026-09-07):
#   root-cause layer: censoring-dependence / positivity-support / maturity
#   performance layer: conditional flags annotated with their root-cause triggers
# Flags are diagnostic, never adjudicative.

.flag_levels <- c("green", "yellow", "red")

new_flag <- function(level = c("green", "yellow", "red"),
                     dimension = c("evaluation", "positivity", "extrapolation"),
                     layer = c("root", "performance"),
                     code, message, evidence = list(),
                     triggered_by = character(0)) {
  level <- match.arg(level)
  dimension <- match.arg(dimension)
  layer <- match.arg(layer)
  structure(
    list(level = level, dimension = dimension, layer = layer,
         code = code, message = message, evidence = evidence,
         triggered_by = triggered_by),
    class = "depcens_flag"
  )
}

# Escalate a flag by n levels (capped at red).
escalate_flag <- function(flag, n = 1L, reason = NULL) {
  i <- match(flag$level, .flag_levels)
  j <- min(i + n, length(.flag_levels))
  if (j > i) {
    flag$level <- .flag_levels[j]
    if (!is.null(reason))
      flag$message <- paste0(flag$message, " [escalated: ", reason, "]")
  }
  flag
}

# Worst flag wins when combining across dimensions.
combine_flags <- function(flags) {
  if (length(flags) == 0L)
    return(new_flag("green", "extrapolation", "performance",
                    code = "ALL_CLEAR",
                    message = "No diagnostic flags raised."))
  lvls <- vapply(flags, function(f) match(f$level, .flag_levels), integer(1))
  flags[[which.max(lvls)]]
}

format_flags <- function(flags) {
  if (length(flags) == 0L) return("  (none)")
  vapply(flags, function(f) {
    trig <- if (length(f$triggered_by))
      paste0(" | triggered by: ", paste(f$triggered_by, collapse = ", ")) else ""
    sprintf("  [%s] %s (%s/%s, %s)%s\n    %s",
            toupper(f$level), f$code, f$dimension, f$layer,
            paste(f$triggered_by, collapse = ""), trig, f$message)
  }, character(1))
}
