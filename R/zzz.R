.onLoad <- function(libname, pkgname) {
  if (requireNamespace("shiny", quietly = TRUE) &&
      requireNamespace("shinyBS", quietly = TRUE)) {
    shiny::addResourcePath("sbs", system.file("www", package = "shinyBS"))
  }
}
