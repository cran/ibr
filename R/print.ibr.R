print.ibr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nInitial df:", format(round(x$initialdf,2)), "; Final df:", format(round(x$finaldf,2)), "\n")
  cat("Number of iterations:", x$iter, "chosen by", x$call$criterion, "\n")
}
