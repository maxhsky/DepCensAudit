# ------------------------------------------------------------------------------
# M3a simulation framework :: runners
#
# One rep = generate data -> run DepCensAudit (3 dimensions) -> record one row
# of flags + continuous statistics. Batches are checkpointed per scenario as
# CSV, so any interrupted run resumes where it stopped.
# ------------------------------------------------------------------------------

# Worst flag level within a flag list (green if empty).
.worst_level <- function(flags, levels = c("green", "yellow", "red")) {
  if (!length(flags)) return("green")
  idx <- vapply(flags, function(f) match(f$level, levels), integer(1))
  levels[max(idx)]
}

.flag_codes <- function(flags) {
  if (!length(flags)) return("")
  paste(vapply(flags, function(f) f$code, character(1)), collapse = ";")
}

# Position-weighted rolling hash (31^i weights) of "scenario_id#rep":
# anagram-safe, dependency-free, deterministic across sessions.
.scenario_seed <- function(scenario_id, rep) {
  s <- as.integer(charToRaw(paste0(scenario_id, "#", rep)))
  h <- 0
  for (b in s) h <- (h * 31 + b) %% 2147483647
  as.integer((h + rep * 2654435761) %% .Machine$integer.max) + 1L
}

# Cox-based prediction function fitted on the OBSERVED data (a realistic user
# workflow); falls back to marginal KM if the Cox fit fails.
.make_pred_fun <- function(obs) {
  km_fit <- survival::survfit(survival::Surv(time, event) ~ 1, data = obs)
  marginal <- function(newdata, t0)
    rep(as.numeric(summary(km_fit, times = t0, extend = TRUE)$surv),
        nrow(newdata))
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ arm, data = obs),
    error = function(e) NULL)
  if (is.null(fit)) return(marginal)
  bh <- survival::basehaz(fit, centered = FALSE)
  beta <- stats::coef(fit)
  function(newdata, t0) {
    s0 <- exp(-stats::approx(bh$time, bh$hazard, xout = t0,
                             method = "constant", rule = 2)$y)
    s0^exp(newdata$arm * beta)
  }
}

# Per-rep engine factory. `generate_scenario` must be in scope (source
# sim/R/generate.R first); DepCensAudit must be installed.
make_rep_runner <- function(horizon = 36) {
  function(sc, n = 1000, rep = 1, seed = NULL) {
    if (is.null(seed)) seed <- .scenario_seed(sc$scenario_id, rep)
    base_row <- data.frame(
      scenario_id = sc$scenario_id, rep = rep, seed = seed, n = n,
      pkg_version = as.character(utils::packageVersion("DepCensAudit")),
      shape = sc$shape, kappa = sc$kappa, maturity = sc$maturity,
      injection = sc$injection,
      truth_A = sc$truth_A, truth_B = sc$truth_B, truth_C = sc$truth_C,
      stringsAsFactors = FALSE)
    res <- tryCatch({
      obs <- generate_scenario(sc, n = n, seed = seed)
      info <- attr(obs, "info")
      pred_fun <- .make_pred_fun(obs)
      audit <- depcens_audit(
        obs, time, event,
        dimensions = c("evaluation", "positivity", "extrapolation"),
        pred_fun = pred_fun, horizon = horizon, seed = seed)
      ev <- audit$dimensions$evaluation
      po <- audit$dimensions$positivity
      ex <- audit$dimensions$extrapolation
      km_at_tau <- po$sensitivity$surv_tau[po$sensitivity$kappa == 0]
      s_true_tau <- true_surv_marginal(po$tau, sc$shape)
      data.frame(
        event_rate = info$event_rate, cens_rate = info$cens_rate,
        ibs = ev$ibs,
        eval_weight_unstable = as.integer(length(ev$flags$root) > 0),
        tau_used = po$tau,
        km_at_tau = km_at_tau,
        s_true_at_tau = s_true_tau,
        km_bias = km_at_tau - s_true_tau,
        pos_covers_truth = as.integer(s_true_tau >= po$interval[1] &&
                                        s_true_tau <= po$interval[2]),
        pos_interval_lo = po$interval[1], pos_interval_hi = po$interval[2],
        pos_interval_width = diff(po$interval),
        pos_flag_codes = .flag_codes(po$flags$root),
        pos_worst_level = .worst_level(po$flags$root),
        rho_instability = ex$metrics$rho_instability,
        rho_maturity = ex$metrics$rho_maturity,
        model_disagreement = ex$metrics$model_disagreement,
        ext_perf_level = ex$flags$performance[[1]]$level,
        ext_perf_code = ex$flags$performance[[1]]$code,
        ext_root_codes = .flag_codes(ex$flags$root),
        eval_flag_codes = .flag_codes(ev$flags$root),
        overall_level = audit$flags$overall$level,
        overall_code = audit$flags$overall$code,
        psa_strategy = audit$psa$strategy,
        config_version = audit$meta$config_version,
        error = "none")
    }, error = function(e) {
      data.frame(event_rate = NA, cens_rate = NA, ibs = NA,
                 eval_weight_unstable = NA,
                 tau_used = NA, km_at_tau = NA, s_true_at_tau = NA,
                 km_bias = NA, pos_covers_truth = NA,
                 pos_interval_lo = NA, pos_interval_hi = NA,
                 pos_interval_width = NA, pos_flag_codes = "",
                 pos_worst_level = NA, rho_instability = NA,
                 rho_maturity = NA, model_disagreement = NA,
                 ext_perf_level = NA, ext_perf_code = "",
                 ext_root_codes = "", eval_flag_codes = "",
                 overall_level = NA, overall_code = "",
                 psa_strategy = "", config_version = NA_character_,
                 error = conditionMessage(e))
    })
    cbind(base_row, res)
  }
}

# Resumable batch: one CSV per scenario under outdir; appends missing reps in
# chunks (default 25) so an interrupted run loses at most one chunk of work.
run_batch <- function(cells, n_reps, n_subj, outdir, runner,
                      verbose = TRUE, chunk = 25) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  for (i in seq_len(nrow(cells))) {
    sc <- cells[i, ]
    f <- file.path(outdir, paste0(sc$scenario_id, ".csv"))
    done <- if (file.exists(f)) sum(count.fields(f, sep = ",")) - 1 else 0
    if (done >= n_reps) next
    if (verbose && done == 0)
      message(sprintf("[%d/%d] %s : starting", i, nrow(cells), sc$scenario_id))
    while (done < n_reps) {
      m <- min(chunk, n_reps - done)
      rows <- lapply((done + 1):(done + m), function(r)
        runner(sc, n = n_subj, rep = r))
      df <- do.call(rbind, rows)
      utils::write.table(df, f, sep = ",", row.names = FALSE,
                         col.names = done == 0, append = done > 0,
                         qmethod = "double")
      done <- done + m
    }
  }
  invisible(outdir)
}
