test_that("psa_config validates thresholds and prints", {
  cfg <- psa_config()
  expect_s3_class(cfg, "psa_config")
  expect_equal(cfg$rho_yellow, 0.10)
  expect_equal(cfg$version, "0.1.0")
  expect_error(psa_config(rho_maturity_red = 0.1))  # red must exceed yellow
  expect_output(print(cfg), "PSA configuration")
})

test_that("psa_mapping implements the calibrated three-level rules", {
  g <- psa_mapping("green")
  expect_equal(g$strategy, "standard_tsd14")
  expect_match(g$diagnostic_wording, "adequately constrained")

  y <- psa_mapping("yellow")
  expect_equal(y$strategy, "structural_mixture_psa")
  expect_match(y$sensitivity_layer, "sensitivity layer only")
  expect_match(y$base_case, "single best curve")

  r <- psa_mapping("red")
  expect_equal(r$strategy, "dual_scenario_no_widening")
  expect_true(any(grepl("tipping-point", r$reporting)))
  # red must never widen intervals
  expect_match(r$sampling, "widening is prohibited")
  # diagnostic, not adjudicative wording
  expect_false(any(grepl("not reimbursable|non-compliant", r$diagnostic_wording)))
})

test_that("flag escalation is capped at red", {
  f <- new_flag("yellow", "extrapolation", "root", code = "X", message = "m")
  f2 <- escalate_flag(f, 5L, reason = "test")
  expect_equal(f2$level, "red")
  expect_match(f2$message, "escalated")
})
