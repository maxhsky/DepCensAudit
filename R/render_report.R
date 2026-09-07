#' Render an audit report in HTML, PDF, or DOCX
#'
#' Renders an audit object (\code{depcens_audit} or any single-dimension
#' audit) through the bundled rmarkdown template. The three formats serve
#' different destinations (HTA advisor input, 2026-09-07):
#'
#' \itemize{
#'   \item \code{"pdf"}: static, printable, cross-referenced output intended
#'     for dossier appendices (primary submission format);
#'   \item \code{"docx"}: for merging into dossier chapters / internal
#'     circulation;
#'   \item \code{"html"}: interactive, for internal reproduction.
#' }
#'
#' Every report embeds the reproducibility metadata: random seed, session
#' info, input-data fingerprint, threshold set version, and diagnostic
#' (never adjudicative) wording.
#'
#' @param audit An audit object.
#' @param output_format One of \code{"html"}, \code{"pdf"}, \code{"docx"}.
#' @param output_file Target file path. Defaults to
#'   \code{"depcens_audit_report.<ext>"} in the working directory.
#' @param ... Passed to \code{rmarkdown::render}.
#'
#' @return The path to the rendered file, invisibly.
#' @export
render_audit_report <- function(audit,
                                output_format = c("html", "pdf", "docx"),
                                output_file = NULL, ...) {
  stopifnot(inherits(audit, c("depcens_audit",
                              "depcens_extrapolation_audit",
                              "depcens_positivity_audit",
                              "depcens_evaluation_audit")))
  output_format <- match.arg(output_format)
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required to render reports. ",
         "Alternatively use print()/summary() for console output.")

  fmt <- switch(output_format,
                html = rmarkdown::html_document(toc = TRUE),
                pdf  = rmarkdown::pdf_document(toc = TRUE),
                docx = rmarkdown::word_document(toc = TRUE))
  if (is.null(output_file))
    output_file <- paste0("depcens_audit_report.",
                          switch(output_format, html = "html",
                                 pdf = "pdf", docx = "docx"))

  template <- system.file(
    "rmarkdown", "templates", "audit-report", "skeleton", "skeleton.Rmd",
    package = "DepCensAudit")
  if (template == "") {
    # Development-mode fallback: running from a source tree (not installed).
    dev_path <- file.path("inst", "rmarkdown", "templates", "audit-report",
                          "skeleton", "skeleton.Rmd")
    if (file.exists(dev_path)) template <- normalizePath(dev_path)
  }
  if (template == "")
    stop("Bundled report template not found; reinstall the package ",
         "or run from the package source root.")

  tmp <- tempfile(fileext = ".Rmd")
  file.copy(template, tmp, overwrite = TRUE)
  rmarkdown::render(
    input = tmp,
    output_format = fmt,
    output_file = basename(output_file),
    output_dir = dirname(normalizePath(output_file, mustWork = FALSE)),
    params = list(audit = audit),
    envir = new.env(parent = globalenv()),
    quiet = TRUE, ...
  )
  invisible(output_file)
}
