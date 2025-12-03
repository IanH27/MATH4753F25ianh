#' Birthday
#'
#' @param x description
#'
#' @returns Odds of duplicate birthdays
#' @export
#'
#' @examples
#' # birthday(20:24)
birthday <- function(x){
  1 - exp(lchoose(365,x) + lfactorial(x) - x*log(365))
}

