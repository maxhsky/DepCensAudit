## DepCensAudit example: NCCTG lung cancer cohort (survival::lung)
##
## Reproducible three-dimension audit on a real, publicly available dataset.
## Run from anywhere after: remotes::install_github("maxhsky/DepCensAudit")
##
## Output: inst/examples/output/lung_audit_report.{html,docx,pdf}

library(DepCensAudit)
library(survival)

## ---- 1. Data ---------------------------------------------------------------
d <- lung[complete.cases(lung[, c("time", "status", "age", "sex", "ph.ecog")]), ]
d$event <- as.integer(d$status == 2)   # lung coding: 1 = censored, 2 = dead
d <- d[, c("time", "event", "age", "sex", "ph.ecog")]
n <- nrow(d)
message("Cohort: NCCTG lung, n = ", n, ", events = ", sum(d$event),
        " (", round(100 * mean(d$event), 1), "%), horizon = 3 years")

## ---- 2. Prediction model for dimension A -----------------------------------
## A routine Cox model; DepCensAudit is model-agnostic (pred_fun interface).
fit <- coxph(Surv(time, event) ~ age + sex + ph.ecog, data = d)
bh  <- basehaz(fit, centered = FALSE)
beta <- coef(fit)

pred_fun <- function(newdata, t0) {
  s0 <- exp(-approx(bh$time, bh$hazard, xout = t0,
                    method = "constant", rule = 2)$y)
  mm <- model.matrix(~ age + sex + ph.ecog, data = newdata)[, -1, drop = FALSE]
  s0^exp(as.numeric(mm %*% beta))
}

## ---- 3. Integrated audit (Phase A) ------------------------------------------
res <- depcens_audit(
  d, time, event,
  dimensions = c("evaluation", "positivity", "extrapolation"),
  pred_fun = pred_fun,
  horizon  = 365 * 3,
  seed     = 20260908
)

print(res)
print(summary(res))

## ---- 4. Render reports ------------------------------------------------------
out_dir <- file.path("inst", "examples", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

html_path <- render_audit_report(
  res, "html", output_file = file.path(out_dir, "lung_audit_report.html"))
message("HTML  report: ", html_path)

docx_path <- render_audit_report(
  res, "docx", output_file = file.path(out_dir, "lung_audit_report.docx"))
message("DOCX  report: ", docx_path)

pdf_ok <- requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex()
if (pdf_ok) {
  pdf_path <- render_audit_report(
    res, "pdf", output_file = file.path(out_dir, "lung_audit_report.pdf"))
  message("PDF   report: ", pdf_path)
} else {
  message("PDF report skipped: no LaTeX installation detected (tinytex).")
}
