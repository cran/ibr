print.ibr <- function(x, digits = max(2, getOption("digits") - 4), ...) {
  cat("\nInitial df:", format(round(x$initialdf,digits)), "; Final df:", format(round(x$finaldf,digits)), "\n")
  cat("Number of iterations:", x$iter, "chosen by", x$call$criterion, "\n")
}
