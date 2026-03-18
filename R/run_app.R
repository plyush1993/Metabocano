#' Run the Metabocano Application
#' Online version: https://plyush1993.shinyapps.io/Metabocano/
#'
#' @export
run_metabocano <- function(...) {
  options(shiny.maxRequestSize = 2 * 1024^3)
  shiny::addResourcePath("www", system.file("www", package = "metabocano"))
  shiny::shinyApp(ui = app_ui(), server = app_server, ...)
}
