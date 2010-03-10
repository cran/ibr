plot.ibr <- function(x,...) {
  r <- residuals(x)
  yh <- predict(x)
  plot(yh, r, xlab = "Index", ylab = "Residuals",...)
  invisible()
}
