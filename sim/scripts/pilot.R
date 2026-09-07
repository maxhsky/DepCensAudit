# ------------------------------------------------------------------------------
# M3a pilot: smoke test of the full simulation pipeline.
# 8 scenario cells x 5 reps at n = 200. Purpose:
#   1. prove the pipeline runs end to end;
#   2. check that truth labels and audit signals co-move in the right
#      direction (kappa > 0 widens the positivity interval; maturity 0.4
#      fires extrapolation flags; orthogonal mechanisms all execute).
# Results land in sim/outputs/pilot/.
# ------------------------------------------------------------------------------

pkg_root <- sub("sim.*$", "", getwd())
for (f in c("sim/R/generate.R", "sim/R/runners.R"))
  source(file.path(pkg_root, f))

library(DepCensAudit)
library(survival)

mg <- make_main_grid()
pick <- c("G_ph_k0.0_m0.8", "G_ph_k0.6_m0.8", "G_late_k0.3_m0.4",
          "G_crossing_k0.6_m0.6")
cells <- rbind(mg[mg$scenario_id %in% pick, ],
               make_orthogonal_grid()[c(1, 5, 9), ])  # ph shape x 3 mechanisms

runner <- make_rep_runner(horizon = 36)
outdir <- file.path(pkg_root, "sim", "outputs", "pilot")
run_batch(cells, n_reps = 5, n_subj = 200, outdir = outdir, runner = runner)

# ---- quick directional summary -----------------------------------------------
files <- list.files(outdir, full.names = TRUE)
all <- do.call(rbind, lapply(files, read.csv))
cat("\n==== pilot directional check (5 reps x 8 cells, n=200) ====\n")
agg <- aggregate(
  cbind(event_rate, km_bias, pos_covers_truth, pos_interval_width,
        rho_maturity, as.integer(ext_perf_level != "green")) ~
    scenario_id + injection + truth_B + truth_C,
  data = all, FUN = mean, na.rm = TRUE)
print(agg, digits = 3)
cat("\nerrors:", sum(all$error != "none", na.rm = TRUE), "of", nrow(all), "reps\n")
