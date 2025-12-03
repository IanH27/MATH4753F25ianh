#' Shade an area under a continuous density curve
#'
#' Generic shading helper for Normal, Gamma, Chi-square, Weibull, Beta, etc.
#'
#' @param xlim Numeric range for plotting.
#' @param densfun Density function (e.g., dnorm).
#' @param pfun CDF function (e.g., pnorm).
#' @param lo,hi Interval endpoints.
#' @param label_digits Number of digits in Label
#' @param label_prefix Prefix for label
#' @param ... Additional arguments passed to densfun/pfun.
#' @export
shade_area <- function(xlim, densfun, pfun, lo, hi, label_digits = 4, label_prefix = "Area = ", ...) {
  curve(densfun(x, ...), xlim = xlim, lwd = 2, ylab = "Density")
  lo_eff <- ifelse(is.finite(lo), lo, xlim[1])
  hi_eff <- ifelse(is.finite(hi), hi, xlim[2])
  xcurve <- seq(lo_eff, hi_eff, length.out = 1000)
  ycurve <- densfun(xcurve, ...)
  polygon(c(lo_eff, xcurve, hi_eff), c(0, ycurve, 0),
          col = rgb(1, 0, 0, 0.35), border = NA)
  prob <- (pfun(hi, ...) - pfun(lo, ...))
  text(mean(c(lo_eff, hi_eff)), 0.9 * max(ycurve),
       paste0(label_prefix, round(prob, label_digits)))
  invisible(prob)
}
