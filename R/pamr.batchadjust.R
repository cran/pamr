
pamr.batchadjust <- function(data) {
  if (is.null(data$batchlabels)) {
    stop("batch labels are not in data object")
  }
  lab <- data$batchlabels
  dd <- model.matrix( ~ factor(lab) - 1)
  data$x <- data$x - misreg.simple(dd, data$x)
  data
}


misreg.simple <- function(Y, x) {
###Y is a indicator response matrix
  nax <- is.na(x)
  nsamples <- (!nax)%*%Y
  x[nax] <- 0
  xsum <- x%*%Y
  xbar <- xsum/nsamples
  xbar %*% t(Y)
}
