#' SHINY MLE
#'
#' @returns shiny app
#' @export
#'
#' @examples
#' \dontrun{shinymle}
shinymle <- function() {
  shiny::runApp(system.file("SHINY", package = "MATH4753F25ianh"), launch.browser = TRUE)
}
