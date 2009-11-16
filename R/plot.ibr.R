plot.ibr <- function(x,...) {
  r <- residuals(x)
  yh <- predict(x)
  plot(yh, r, xlab = "Fit", ylab = "Residuals",...)
  invisible()
}
