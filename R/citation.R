#' @title Citation for rnaseqExplorer
#' @description Provides the citation information for the rnaseqExplorer package.
#' @export
citation.rnaseqExplorer <- function() {
  cat("Jeongah Lee. (2024). rnaseqExplorer: A comprehensive Shiny application for analyzing RNA sequencing data. R package version 0.1.0. https://github.com/jeongahblairlee/rnaseqExplorer\n")
}

#' @export
citation <- function(package = "rnaseqExplorer") {
  if (package == "rnaseqExplorer") {
    citation.rnaseqExplorer()
  } else {
    methods::citation(package)
  }
}
