#' A function to plot the cross-validated sample probabilities from the nearest
#' shrunken centroid classifier
#' 
#' A function to plot the cross-validated sample probabilities from the nearest
#' shrunken centroid classifier
#' 
#' \code{pamr.plotcvprob} plots the cross-validated sample probabilities the
#' from nearest shrunken centroid classifier, stratified by the true classses.
#' 
#' @param fit The result of a call to pamr.cv
#' @param data A list with at least two components: x- an expression genes in
#' the rows, samples in the columns), and y- a vector of the class labels for
#' each sample. Same form as data object used by pamr.train.
#' @param threshold Threshold value to be used
#' @author Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' mydata <- list(x=x,y=y)
#' mytrain <-   pamr.train(mydata)
#' mycv <-  pamr.cv(mytrain,mydata)
#' pamr.plotcvprob(mycv,mydata,threshold=1.6)
#' 
#' 
#' 
#' 
#' @export pamr.plotcvprob
pamr.plotcvprob <- function(fit, data, threshold) {
  par(pch = 1)
  ii <- (1:length(fit$threshold))[fit$threshold > threshold]
  ii <- ii[1]
  ss <- data$samplelabels
  pp <- fit$prob[,  , ii]
  if(is.null(fit$newy)) {
    y <- fit$y[fit$sample.subset]
  }
  if(!is.null(fit$newy)) {
    y <- fit$newy[fit$sample.subset]
  }
  o <- order(y)
  y <- y[o]
  if(!is.null(ss)) {
    ss <- ss[o]
  }
  ppp <- pp[o,  ]
  n <- nrow(ppp)
  nc <- length(unique(y))
  par(cex = 1)
  plot(1:n, ppp[, 2], type = "n", xlab = "sample", ylab = 
       "cross-validated probabilities", ylim = c(0, 1.2), axes = FALSE)
  axis(1)
  axis(2, at=seq(0, 1.2, by=0.2), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0", ""))
  axis(4)
  for(j in 1:nc) {
    points(1:n, ppp[, j], col = j + 1)
  }
  for(j in 1:(nc - 1)) {
    abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
  }
  h <- c(0, table(y))
  for(j in 2:(nc + 1)) {
    text(sum(h[1:(j - 1)]) + 0.5 * h[j], 1.02, label = levels(y)[j - 
                                                 1], col = j)
  }
  abline(h = 1)
  if(!is.null(ss)) {
    text(1:length(ss), 1.1, labels = ss, srt = 90, cex = 0.7)
  }
  ##if(!is.null(ss)){axis(3,labels=ss,at=1:length(ss),srt=90)}
}

