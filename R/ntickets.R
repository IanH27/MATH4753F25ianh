#' Airline Overbooking Problem Function
#'
#' Calculates the number of tickets an airline should sell when the number of
#' seats is N, the probability of a show-up is p, and gamma is the maximum
#' acceptable probability of overbooking (too many people showing up).
#'
#' @param N Integer. Number of seats on the plane.
#' @param gamma Numeric. Target overbooking probability (e.g., 0.02 for 2 Percent).
#' @param p Numeric. Probability a passenger shows up (e.g., 0.95).
#'
#' @return A named list with the following components:
#'   \describe{
#'     \item{nd}{Number of tickets using the discrete Binomial method}
#'     \item{nc}{Number using the Normal approximation}
#'     \item{N}{Input number of seats}
#'     \item{p}{Probability of showing up}
#'     \item{gamma}{Overbooking probability}
#'   }
#'
#' @details
#' The function also creates plots showing the objective function for both
#' discrete and normal approximation cases.
#'
#' @examples
#' ntickets(N = 400, gamma = 0.02, p = 0.95)
#'
#' @importFrom stats pbinom pnorm
#' @importFrom graphics abline plot par
#' @export

ntickets <- function(N, gamma, p) {
  nd <- NA
  for (n in seq(N, N + 100)) {
    prob_overbook <- 1 - pbinom(N, size = n, prob = p)
    if (prob_overbook <= gamma) {
      nd <- n
      break
    }
  }
  nc <- NA
  for (n in seq(N, N + 100)) {
    mu <- n * p
    sigma <- sqrt(n * p * (1 - p))
    prob_overbook_norm <- 1 - pnorm(N + 0.5, mean = mu, sd = sigma)
    if (prob_overbook_norm <= gamma) {
      nc <- n
      break
    }
  }
  n_vals <- seq(N, N + 50)
  obj_discrete <- sapply(n_vals, function(n) (1 - gamma) - pbinom(N, size = n, prob = p))
  obj_normal <- sapply(n_vals, function(n) {
    mu <- n * p
    sigma <- sqrt(n * p * (1 - p))
    (1 - gamma) - pnorm(N + 0.5, mean = mu, sd = sigma)
  })
  par(mfrow = c(1, 2))
  plot(n_vals, obj_discrete, type = "l", col = "blue",
       main = "Objective Function (Discrete Binomial)",
       xlab = "n (tickets sold)", ylab = "Objective Value")
  abline(h = 0, col = "red", lty = 2)
  abline(v = nd, col = "darkgreen", lwd = 2)
  plot(n_vals, obj_normal, type = "l", col = "purple",
       main = "Objective Function (Normal Approximation)",
       xlab = "n (tickets sold)", ylab = "Objective Value")
  abline(h = 0, col = "red", lty = 2)
  abline(v = nc, col = "darkgreen", lwd = 2)
  result <- list(nd = nd, nc = nc, N = N, p = p, gamma = gamma)
  return(result)
}
