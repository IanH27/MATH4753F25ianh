#' Bootstrap Confidence Interval Function
#'
#' This function performs bootstrap resampling and constructs
#' confidence intervals for a user-specified statistic.
#'
#' @param iter Number of bootstrap iterations.
#' @param x Numeric vector of data.
#' @param fun Statistic to compute (e.g., "mean", "var", median, IQR).
#' @param alpha Level of significance.
#' @param cx Text scaling for labels.
#' @param ... Additional graphical parameters for hist().
#'
#' @return A list containing:
#' \item{ci}{Bootstrap confidence interval}
#' \item{xstat}{Bootstrap statistics}
#' \item{fun}{Statistic used}
#' \item{x}{Original sample}
#'
#' @examples
#' # myboot2(x = rnorm(20))
#'
#' @export

myboot2 <- function(iter = 10000,
                    x,
                    fun = "mean",
                    alpha = 0.05,
                    cx = 1.5,
                    ...) {  # Notice where the ... is repeated
  n <- length(x)   # sample size

  # Now sample with replacement
  y <- sample(x, n * iter, replace = TRUE)  # A

  # Make a matrix with all the resampled values
  rs.mat <- matrix(y, nr = n, nc = iter, byrow = TRUE)

  # Compute the statistic for each bootstrap resample
  xstat <- apply(rs.mat, 2, fun)

  # Quantile-based confidence interval
  ci <- quantile(xstat, c(alpha/2, 1 - alpha/2)) # B

  # Histogram of bootstrap statistics
  main.title <- paste("Histogram of Bootstrap sample statistics",
                      "\n", "alpha = ", alpha,
                      " iter = ", iter, sep = "")
  para <- hist(xstat,
               freq = FALSE,
               las  = 1,
               main = main.title,
               ...)

  # Matrix version of original data so we can use apply()
  mat <- matrix(x, nr = length(x), nc = 1, byrow = TRUE)

  # Point estimate using the same function
  pte <- apply(mat, 2, fun)

  # Add vertical line for point estimate
  abline(v = pte, lwd = 3, col = "black")

  # Segment for confidence interval
  segments(ci[1], 0, ci[2], 0, lwd = 4)
  text(ci[1], 0,
       paste("(", round(ci[1], 2), sep = ""),
       col = "red", cex = cx)
  text(ci[2], 0,
       paste(round(ci[2], 2), ")", sep = ""),
       col = "red", cex = cx)

  # Plot the point estimate half-way up the density
  text(pte, max(para$density) / 2, round(pte, 2), cex = cx)

  # Return useful output (including xstat for later tasks)
  invisible(list(ci = ci,
                 fun = fun,
                 x = x,
                 xstat = xstat))
}
