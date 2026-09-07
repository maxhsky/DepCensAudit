sim_cohort <- function(n = 400, rate = 0.08, cens_rate = 0.03, seed = 1) {
  set.seed(seed)
  t_event <- rexp(n, rate)
  t_cens <- rexp(n, cens_rate)
  data.frame(time = pmin(t_event, t_cens),
             event = as.integer(t_event <= t_cens))
}

test_that("extrapolation audit runs end-to-end on a simulated cohort", {
  skip_if_not_installed("flexsurv")
  d <- sim_cohort()
  a <- audit_extrapolation(d, time, event, horizon = 30, seed = 42)
  expect_s3_class(a, "depcens_extrapolation_audit")
  expect_true(is.finite(a$metrics$rho_instability))
  expect_true(a$metrics$rho_maturity >= 0)
  expect_equal(length(a$fits), 3L)
  expect_true(a$flags$performance[[1]]$level %in% c("green", "yellow", "red"))
  expect_true(nchar(a$meta$data_fingerprint) > 0)
})

test_that("frozen .interaction interface errors in Phase A", {
  d <- sim_cohort()
  expect_error(
    audit_positivity(d, time, event, .interaction = list(feedback = TRUE)),
    "frozen interface")
  expect_error(
    audit_extrapolation(d, time, event, horizon = 30,
                        .interaction = list(feedback = TRUE)),
    "frozen interface")
  skip_if_not_installed("flexsurv")
  expect_error(
    depcens_audit(d, time, event, horizon = 30,
                  dimensions = "extrapolation",
                  .interaction = list(feedback = TRUE)),
    "frozen interface")
})

test_that("positivity audit returns a sensitivity interval, not a verdict", {
  d <- sim_cohort()
  b <- audit_positivity(d, time, event)
  expect_s3_class(b, "depcens_positivity_audit")
  expect_equal(nrow(b$sensitivity), length(b$meta$kappa_grid))
  expect_true(all(diff(b$interval) >= 0))
  # no performance-layer verdict on the censoring mechanism
  expect_equal(length(b$flags$performance), 0L)
})

test_that("evaluation audit computes IPCW IBS from a prediction function", {
  d <- sim_cohort()
  km <- survival::survfit(survival::Surv(time, event) ~ 1, data = d)
  pred_fun <- function(newdata, t0) {
    rep(as.numeric(summary(km, times = t0, extend = TRUE)$surv),
        nrow(newdata))
  }
  a <- audit_evaluation(d, time, event, pred_fun = pred_fun)
  expect_s3_class(a, "depcens_evaluation_audit")
  expect_true(a$ibs >= 0 && a$ibs <= 1)
  expect_equal(nrow(a$brier), 25L)
})
