# DepCensAudit

<!-- badges: start -->
[![R-CMD-check](https://github.com/maxhsky/DepCensAudit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/maxhsky/DepCensAudit/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Integrated **diagnostic audit** framework for censored survival data:
evaluation × positivity × extrapolation.

> **Phase A release (v0.0.1, engineering integration).** This package chains
> three single-dimension audits into one reproducible workflow. It makes
> **no claim** of statistical performance superiority over manually chaining
> single-dimension tools; the cross-dimensional interaction pathway is a
> frozen interface (`.interaction`) reserved for a later methodology release.
> All outputs are **diagnostic** ("data support is insufficient to constrain
> tail behaviour"), never adjudicative.

## Quick example (NCCTG lung cohort)

A fully reproducible three-dimension audit on `survival::lung` (n = 227,
72.2% events, 3-year horizon, Cox model for the evaluation dimension):

```r
library(DepCensAudit)
library(survival)

d <- lung[complete.cases(lung[, c("time","status","age","sex","ph.ecog")]), ]
d$event <- as.integer(d$status == 2)

fit  <- coxph(Surv(time, event) ~ age + sex + ph.ecog, data = d)
bh   <- basehaz(fit, centered = FALSE)
beta <- coef(fit)
pred_fun <- function(newdata, t0) {
  s0 <- exp(-approx(bh$time, bh$hazard, xout = t0,
                    method = "constant", rule = 2)$y)
  mm <- model.matrix(~ age + sex + ph.ecog, data = newdata)[, -1, drop = FALSE]
  s0^exp(as.numeric(mm %*% beta))
}

res <- depcens_audit(d, time, event,
                     pred_fun = pred_fun, horizon = 365 * 3,
                     seed = 20260908)
render_audit_report(res, "html")
```

Real-data result (2026-09-08 run): evaluation GREEN, positivity
**YELLOW** (`SPARSE_TAIL_RISKSET` — only 1.3% of subjects at risk in the
final 10% of follow-up), extrapolation GREEN → PSA strategy
`standard_tsd14`. Full script and rendered reports:
[`inst/examples/`](inst/examples/).

## What it audits

| Dimension | Function | Output | Basis |
|-----------|----------|--------|-------|
| A. Evaluation | `audit_evaluation()` | IPCW Brier / IBS (model-agnostic `pred_fun` interface) | reuses the IBS-Dep idea (DependentEVAL, arXiv:2502.19460). **Experimental (work item A-1)**: operating characteristics as an audit signal under validation |
| B. Positivity | `audit_positivity()` | **continuous sensitivity interval** for S(τ) under tilting κ ∈ [0, 0.6]; flags only detectable positivity sub-problems (sparse tail risk sets) | censoring independence is untestable → no traffic-light verdict on the mechanism |
| C. Extrapolation | `audit_extrapolation()` | three-model comparison (RP-2knot / parametric / constant-tail), ρ instability + ρ maturity + model disagreement, PSA strategy mapping | calibrated thresholds (`psa_config()`, versioned); NMA-context evidence is labelled as such |

## Install (development)

```r
# remotes::install_local("DepCensAudit")   # from this source tree
# Requires: survival, flexsurv; Suggests: rmarkdown (reports), testthat
```

## Quick start

```r
library(DepCensAudit)

d <- survival::lung   # single-cohort IPD only in v0.0.x
d$event <- as.integer(d$status == 2)

# Dimension C: extrapolation audit
aud <- audit_extrapolation(d, time, event, horizon = 1000, seed = 1)
aud$psa            # mapped PSA strategy (green/yellow/red rules)

# Full pipeline (dimension A needs a prediction function)
fit <- survival::coxph(survival::Surv(time, event) ~ age + sex, data = d)
pred_fun <- function(newdata, t0) {
  sf <- survival::survfit(fit, newdata = newdata)
  vapply(seq_len(nrow(newdata)),
         function(i) as.numeric(summary(sf[i], times = t0,
                                        extend = TRUE)$surv),
         numeric(1))
}
full <- depcens_audit(d, time, event, horizon = 1000,
                      pred_fun = pred_fun, seed = 1)

# Three-format report: pdf (dossier appendix) / docx (chapter merge) / html
render_audit_report(full, "pdf",  "audit_report.pdf")
render_audit_report(full, "docx", "audit_report.docx")
render_audit_report(full, "html", "audit_report.html")
```

## HTA-facing design decisions

- **Red flag ≠ widened PSA interval.** Red triggers per-scenario full
  cost-effectiveness results, stratified CEACs and an ICER tipping-point
  analysis instead.
- **Expert-elicitation priors (TSD26) mount in the sensitivity layer only**;
  they never re-fit the base case.
- Thresholds live in a versioned `psa_config()` object and must be
  calibrated against the project simulation study before regulatory-facing
  use. The config version is embedded in every report.
- Reports embed seed, session info, data fingerprint and diagnostic wording.

## Scope limits (v0.0.x)

- Single-cohort IPD only; NMA-level extrapolation is out of scope.
- Competing risks are out of scope (planned v1.x; empirical demonstration
  datasets with high competing mortality are used as illustrations only).
- `constail` is a simplified constant-hazard tail for single cohorts, not
  the NMA-grade constant-tail-HR variant.

## Reproducibility

```r
# devtools::test()   # testthat suite: PSA rules, frozen interface, audits
```

License: MIT.
