#' Plot Normal curve and shade P(X <= a)
#'
#' Draws a Normal density curve, shades the region from -Inf to a,
#' displays the probability on the plot, and returns results as a list.
#'
#' @param mu Mean of the normal distribution.
#' @param sigma Standard deviation (>0).
#' @param a Cutoff where the area is computed (P(X <= a)).
#'
#' @return A named list with mu, sigma, a, and prob.
#' @examples
#' myncurve(mu = 10, sigma = 5, a = 6)
#' @export
myncurve <- function(mu, sigma, a) {
  stopifnot(is.finite(mu), is.finite(sigma), sigma > 0, is.finite(a))
  xlim <- c(mu - 4*sigma, mu + 4*sigma)
  curve(dnorm(x, mean = mu, sd = sigma), xlim = xlim, lwd = 2,
        ylab = "Density", main = sprintf("N(%.3f, %.3f^2)", mu, sigma))
  lo_eff <- xlim[1]
  hi_eff <- min(a, xlim[2])
  xcurve <- seq(lo_eff, hi_eff, length.out = 1000)
  ycurve <- dnorm(xcurve, mean = mu, sd = sigma)
  polygon(c(lo_eff, xcurve, hi_eff), c(0, ycurve, 0),
          col = rgb(0, 0, 1, 0.35), border = NA)
  prob <- pnorm(a, mean = mu, sd = sigma)
  text(mean(c(lo_eff, hi_eff)), 0.9 * max(ycurve),
       paste0("P = ", round(prob, 4)))
  list(mu = mu, sigma = sigma, a = a, prob = prob)
}
