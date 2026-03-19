#' Run the Metabocano Application
#'
#' @importFrom green bold blue
#' @export
run_metabocano <- function(...) {
  cat("\n")
  cat(crayon::green("             +--------------------+\n"))
  app_name <- paste0(
    crayon::green(crayon::bold("    Metabocano    "))
  )
  cat(crayon::green("             | "), app_name, crayon::green(" |\n"), sep = "")
  cat(crayon::green("             +--------------------+\n"))
  cat("\n")
  cat(crayon::blue(crayon::bold("Enhanced interactive volcano plot for metabolomics studies\n")))
  cat("\n")

  flush.console()

  old_opts <- options(shiny.maxRequestSize = 2 * 1024^3)
    on.exit({
    options(old_opts)
    gc()
    }, add = TRUE)
  shiny::addResourcePath("www", system.file("www", package = "metabocano"))
  shiny::shinyApp(ui = app_ui(), server = app_server, ...)
}
