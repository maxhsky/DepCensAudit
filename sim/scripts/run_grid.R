# ------------------------------------------------------------------------------
# M3a full batch launcher (parallel, resumable).
#
#   Rscript sim/scripts/run_grid.R [n_reps] [n_subj] [outdir_suffix] [workers]
#
# Defaults: 500 reps, n = 1000 (main + orthogonal), n = 200 (sensitivity
# layer), workers = 4. Checkpointed per scenario under sim/outputs/<suffix>/
# - kill and re-run any time; it picks up where it stopped.
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
n_reps  <- as.integer(if (length(args) >= 1) args[1] else 500)
n_subj  <- as.integer(if (length(args) >= 2) args[2] else 1000)
suffix  <- if (length(args) >= 3) args[3] else "grid"
workers <- as.integer(if (length(args) >= 4) args[4] else 4)

pkg_root <- sub("sim.*$", "", normalizePath(getwd()))
outdir <- file.path(pkg_root, "sim", "outputs", suffix)

for (f in c("sim/R/generate.R", "sim/R/runners.R"))
  source(file.path(pkg_root, f), local = TRUE)
library(DepCensAudit); library(survival)

cells <- rbind(make_main_grid(), make_orthogonal_grid())
sens_cells <- make_sensitivity_grid()
runner <- make_rep_runner(horizon = 36)

t0 <- Sys.time()
if (workers > 1) {
  chunks <- split(seq_len(nrow(cells)),
                  cut(seq_len(nrow(cells)), workers, labels = FALSE))
  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("pkg_root", "outdir", "n_reps", "n_subj",
                                "cells", "chunks", "runner"))
  invisible(parallel::clusterEvalQ(cl, {
    source(file.path(pkg_root, "sim", "R", "generate.R"), local = TRUE)
    source(file.path(pkg_root, "sim", "R", "runners.R"),  local = TRUE)
    library(DepCensAudit); library(survival)
  }))
  invisible(parallel::parLapply(cl, chunks, function(idx) {
    run_batch(cells[idx, , drop = FALSE], n_reps, n_subj,
              outdir, runner, verbose = FALSE)
    NULL
  }))
} else {
  run_batch(cells, n_reps, n_subj, outdir, runner)
}

# Sensitivity layer: n = 200, sequential (cheap).
run_batch(sens_cells, n_reps, 200, outdir, runner, verbose = FALSE)

cat("\nDone in", round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 1),
    "hours. Results in", outdir, "\n")
