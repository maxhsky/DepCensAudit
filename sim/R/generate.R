# ------------------------------------------------------------------------------
# M3a simulation framework :: data generation
#
# Design (project doc v1.0 r2, section 5; reviewer round 2 amendments):
#   * Main grid: non-PH shape {ph, early, late, crossing}
#                x censoring dependence kappa {0, 0.3, 0.6}   (Kendall's tau)
#                x follow-up maturity {0.40, 0.60, 0.80}
#     -> 36 cells; kappa is realised through a GAUSSIAN copula coupling the
#        latent event and censoring uniforms (a DIFFERENT family from the
#        audit's exponential-tilting detector -> honest "heterologous"
#        performance, reviewer item 2.3).
#   * Orthogonal injection sublayer (non-copula mechanisms, >=50% of injected
#     scenarios overall, per pre-registration commitment):
#       1. logit_dropout  - dropout driven by an unobserved prognostic frailty
#       2. comp_miscode   - competing deaths miscoded as censoring
#       3. staggered      - staggered entry x administrative censoring
#   * Sensitivity layer: n = 200 at kappa in {0, 0.6} x maturity 0.6 x 4 shapes.
#
# Event process: two arms, piecewise-constant HR shapes, exponential control
# (median 12 months). Analytic cumulative hazards -> exact inversion.
# ------------------------------------------------------------------------------

## --- scenario grids -----------------------------------------------------------

shape_levels <- c("ph", "early", "late", "crossing")

hr_shape <- function(t, shape) {
  switch(shape,
    ph       = rep(0.7, length(t)),
    early    = 1.6 * exp(-t / 6) + 0.4,
    late     = 0.4 + 1.4 * (1 - exp(-t / 6)),
    crossing = ifelse(t < 12, 0.6, 1.5),
    stop("unknown shape: ", shape))
}

# Cumulative HR integral (analytic), shape x.
cumHR <- function(t, shape) {
  switch(shape,
    ph       = 0.7 * t,
    early    = 9.6 * (1 - exp(-t / 6)) + 0.4 * t,
    late     = 1.8 * t + 8.4 * (exp(-t / 6) - 1),
    crossing = ifelse(t < 12, 0.6 * t, 0.6 * 12 + 1.5 * (t - 12)),
    stop("unknown shape: ", shape))
}

make_main_grid <- function() {
  g <- expand.grid(shape = shape_levels, kappa = c(0, 0.3, 0.6),
                   maturity = c(0.4, 0.6, 0.8))
  g$injection <- "gauss_copula"
  g$scenario_id <- sprintf("G_%s_k%.1f_m%.1f", g$shape, g$kappa, g$maturity)
  # Truth labels (draft; finalised in SAP before unblinding)
  g$truth_B <- g$kappa > 0        # dependent censoring present
  g$truth_C <- g$maturity == 0.4  # immature follow-up vs 36-month horizon
  g$truth_A <- FALSE
  g
}

make_orthogonal_grid <- function(maturity = 0.6) {
  mech <- c("logit_dropout", "comp_miscode", "staggered")
  g <- expand.grid(shape = shape_levels, mechanism = mech)
  g$kappa <- 0.3          # nominal label (dependence is mechanism-specific)
  g$maturity <- maturity
  g$injection <- as.character(g$mechanism)
  g$scenario_id <- sprintf("O_%s_%s", g$mechanism, g$shape)
  g$truth_B <- TRUE       # all three induce dependent censoring
  g$truth_C <- FALSE
  g$truth_A <- g$mechanism == "comp_miscode"  # breaks IPCW assumptions
  g$mechanism <- NULL                          # align columns with main grid
  g
}

make_sensitivity_grid <- function() {
  g <- expand.grid(shape = shape_levels, kappa = c(0, 0.6))
  g$maturity <- 0.6
  g$injection <- "gauss_copula"
  g$scenario_id <- sprintf("S_%s_k%.1f_m%.1f", g$shape, g$kappa, g$maturity)
  g$truth_B <- g$kappa > 0
  g$truth_C <- FALSE
  g$truth_A <- FALSE
  g
}

## --- data generation ----------------------------------------------------------

.lambda0 <- log(2) / 12   # control-arm hazard; median survival 12 months
.t_grid   <- seq(0, 240, by = 0.05)   # inversion grid (months)

# Vectorised event-time draw for one arm via exact inversion.
.draw_event <- function(u, arm, shape) {
  if (!arm) return(-log(u) / .lambda0)
  H <- .lambda0 * cumHR(.t_grid, shape)
  stats::approx(H, .t_grid, xout = -log(u), rule = 2)$y
}

# Censoring-scale constant: Weibull(shape 1.3, scale 30) gives a light
# (roughly 15-25%) marginal random-censoring load across scenarios.
.cens_shape <- 1.3
.cens_scale <- 30

# Analytic marginal true survival S(t) over the two balanced arms.
# Used ONLY for scoring (bias / coverage) in simulation - never as a
# detection input (detection rules must consume tool output alone).
true_surv_marginal <- function(t, shape) {
  s0 <- exp(-.lambda0 * t)
  s1 <- exp(-.lambda0 * cumHR(t, shape))
  0.5 * (s0 + s1)
}

generate_scenario <- function(sc, n = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  arm <- rep(0:1, each = floor(n / 2))
  if (length(arm) < n) arm <- c(arm, sample(0:1, 1))

  ## -- latent uniforms, coupled per injection mechanism ------------------------
  U1 <- runif(n)                      # event latent
  U2 <- runif(n)                      # censoring latent
  V  <- rnorm(n)                      # unobserved prognostic frailty
  switch(sc$injection,
    gauss_copula = {
      rho <- sin(pi * sc$kappa / 2)   # kappa = Kendall's tau -> normal rho
      z1 <- qnorm(U1)
      U2 <- pnorm(rho * z1 + sqrt(1 - rho^2) * rnorm(n))
    },
    logit_dropout = {
      # Frailty V accelerates failure and drives dropout; NOT a copula family.
      U1 <- pnorm(qnorm(U1) - 0.4 * V)          # frail patients fail earlier
    },
    comp_miscode = NULL,              # handled below via competing hazard
    staggered   = NULL,               # handled below via calendar time
    stop("unknown injection: ", sc$injection))

  ## -- full event times --------------------------------------------------------
  Tfull <- ifelse(arm == 1,
                  .draw_event(U1, TRUE, sc$shape),
                  .draw_event(U1, FALSE, sc$shape))
  if (sc$injection == "logit_dropout")
    Tfull <- Tfull * exp(-0.3 * V)    # AFT frailty component

  ## -- censoring times ---------------------------------------------------------
  Cfull <- stats::qweibull(U2, shape = .cens_shape, scale = .cens_scale)
  if (sc$injection == "logit_dropout")
    Cfull <- stats::rexp(n, rate = 1 / 24) * exp(-0.5 * V)  # frail drop out early

  ## -- administrative truncation ------------------------------------------------
  if (sc$injection == "staggered") {
    entry <- runif(n, 0, 18)          # staggered entry over 18 months
    close <- 36                       # trial closes at calendar month 36
    t_admin <- pmax(close - entry, 0.1)
  } else {
    t_admin <- rep(stats::quantile(Tfull, probs = sc$maturity, type = 1), n)
  }

  ## -- compete and assemble observed data ---------------------------------------
  if (sc$injection == "comp_miscode") {
    Tcomp <- stats::rweibull(n, shape = 1.3, scale = 30)
    time  <- pmin(pmin(Tfull, Cfull), pmin(t_admin, Tcomp))
    event <- as.integer(Tfull <= pmin(Cfull, t_admin) & Tfull <= Tcomp)
  } else {
    time  <- pmin(pmin(Tfull, Cfull), t_admin)
    event <- as.integer(Tfull <= pmin(Cfull, t_admin))
  }

  obs <- data.frame(time = round(time, 6), event = event, arm = arm)
  attr(obs, "info") <- list(
    scenario_id = sc$scenario_id, n = n, shape = sc$shape,
    kappa = sc$kappa, maturity = sc$maturity, injection = sc$injection,
    seed = seed,
    event_rate = mean(event), cens_rate = 1 - mean(event),
    t_admin_median = stats::median(t_admin))
  obs
}
