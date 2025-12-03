#' Sampling distribution of the mean for Uniform(a,b)
#'
#' Generates `iter` samples of size `n` from Uniform(a,b) and plots the
#' sampling distribution of the sample mean, overlaying a kernel density
#' and the theoretical normal curve from the CLT.
#'
#' @param n Sample size.
#' @param iter Number of iterations (replications).
#' @param a,b Range of the Uniform distribution.
#' @return Invisibly returns a numeric vector of sample means.
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   w <- mycltu(n = 20, iter = 10000, a = 0, b = 10)
#'   hist(w)
#' }
#' @export
mycltu <- function(n, iter, a = 0, b = 10) {
  y <- runif(n * iter, a, b)
  data <- matrix(y, nrow = n, ncol = iter, byrow = TRUE)
  w <- apply(data, 2, mean)
  param <- hist(w, plot = FALSE)
  ymax <- 1.1 * max(param$density)
  hist(w, freq = FALSE, ylim = c(0, ymax),
       main = paste("Histogram of sample mean\nsample size = ", n, sep = ""),
       xlab = "Sample mean")
  lines(density(w), col = "blue", lwd = 3)
  curve(dnorm(x, mean = (a + b)/2, sd = (b - a)/sqrt(12 * n)),
        add = TRUE, col = "red", lty = 2, lwd = 3)
  curve(dunif(x, a, b), add = TRUE, lwd = 4)
  invisible(w)
}
