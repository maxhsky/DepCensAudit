# M3a Pilot Findings (2026-09-08, 8 cells x 5 reps, n = 200)

Pilot purpose: pipeline smoke test + directional checks. **Not confirmatory.**
All findings below are inputs to the SAP freeze, not performance claims.

## P1. Framework runs end to end — 0 errors / 35 reps (after fixes)

- Generation -> 3-dimension audit -> CSV checkpointing works.
- Error-column bug fixed: empty strings read as logical NA inflated error
  counts (reported 35/35 false errors); runner now writes `"none"` sentinel.
- Orthogonal column-mismatch bug fixed (mechanism column dropped to align
  with main grid).

## P2. Copula dependence is genuinely generated (validated)

Empirical Kendall's tau at generation kappa: 0.029 / -0.278 / -0.587.
KM bias at median-event tau scales monotonically:
-0.004 (k=0) / -0.032 (k=0.3) / -0.066 (k=0.6). Injection mechanism valid.

## P3. One-sided tilting interval was a real defect — FIXED in v0.0.2

Original default `kappa_grid = seq(0, 0.6)` only tilted upward; whenever KM
overestimated the truth (clean-cell noise, logit_dropout mechanisms), the
truth fell outside the interval (coverage 0.2–0.4). Default is now the
symmetric `seq(-0.6, 0.6, by = 0.2)` — the direction of untestable
dependence is unknown in practice, so the sensitivity interval is two-sided
by design. Post-fix coverage: clean cell 1.0, copula cells 0.8–1.0,
logit_dropout cells 0.6 (heterologous mechanism; characterize at n = 1000).

## P4. OPEN — interval width does not separate kappa levels

Positivity interval width is nearly flat in kappa (0.085–0.097 across cells).
The draft detection rule `det_B: width >= w_B` will therefore have
sensitivity near zero. SAP implication: dimension B should be framed as
*sensitivity quantification* (interval coverage / KM-bias response under
specified tilt), consistent with reviewer item 2.2 — NOT as a dependent-
censoring *detector*. Estimands for H1/H3 in dimension B need redefinition
before freezing.

## P5. OPEN — truth labels for dimension C are mis-graded

The clean m=0.8 cell fires PARTIAL_MATURITY in 3/5 reps. Analysis: at
m=0.8 the analytic administrative cutoff lands at ~33/36 months of the
horizon with ~23% of expected survival mass beyond data support — the audit
is arguably RIGHT and the draft truth label (`truth_C = maturity == 0.4`)
is WRONG. SAP fix: grade truth_C analytically per cell from
`true_surv_marginal()` (expected-mass-beyond-cutoff ratio), not from grid
cell membership.

## Next actions

1. SAP draft: redefine dimension-B estimands (P4) and analytic truth-C
   grading (P5); freeze before any confirmatory batch.
2. Confirmatory batch: n = 1000, >= 500 reps, all 36 cells + orthogonal
   sublayer (`sim/scripts/run_grid.R`).
3. Characterize logit_dropout coverage at n = 1000 (P3 residual).
