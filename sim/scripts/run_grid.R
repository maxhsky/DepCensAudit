# ------------------------------------------------------------------------------
# M3a full batch launcher (resumable; run in chunks if needed).
#
#   Rscript sim/scripts/run_grid.R [n_reps] [n_subj] [outdir_suffix]
#
# Defaults: 500 reps, n = 1000 (main + orthogonal), n = 200 (sensitivity
# layer). Checkpointed per scenario under sim/outputs/<suffix>/ - kill and
# re-run any time; it picks up where it stopped.
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
n_reps <- as.integer(if (length(args) >= 1) args[1] else 500)
n_subj <- as.integer(if (length(args) >= 2) args[2] else 1000)
suffix <- if (length(args) >= 3) args[3] else "grid"

pkg_root <- sub("sim.*$", "", getwd())
for (f in c("sim/R/generate.R", "sim/R/runners.R"))
  source(file.path(pkg_root, f))

library(DepCensAudit)
library(survival)

runner <- make_rep_runner(horizon = 36)
outdir <- file.path(pkg_root, "sim", "outputs", suffix)

run_batch(make_main_grid(), n_reps, n_subj, outdir, runner)
run_batch(make_orthogonal_grid(), n_reps, n_subj, outdir, runner)
run_batch(make_sensitivity_grid(), n_reps, 200, outdir, runner)

cat("\nDone. Results in", outdir, "\n")
