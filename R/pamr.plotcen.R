#' A function to plot the shrunken class centroids, from the nearest shrunken
#' centroid classifier
#' 
#' A function to plot the shrunken class centroids, from the nearest shrunken
#' centroid classifier produced by pamr.train
#' 
#' \code{pamr.plotcen} plots the shrunken class centroids for each class, for
#' genes surviving the threshold for at least once class. If genenames are
#' included in "data", they are added to the plot. Note: for many classes and
#' long gene names, this plot may need some manual prettying.
#' 
#' @param fit The result of a call to pamr.train
#' @param data The input data, in the same form as that used by pamr.train
#' @param threshold The desired threshold value
#' @author Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' mydata <- list(x=x,y=y,genenames=as.character(1:1000))
#' mytrain <-   pamr.train(mydata)
#' mycv <- pamr.cv(mytrain,mydata)
#' pamr.plotcen(mytrain, mydata,threshold=1.6)
#'  
#' 
#' @export pamr.plotcen
pamr.plotcen <- function(fit, data, threshold) {
  genenames <- data$genenames[fit$gene.subset]
  x <- data$x[fit$gene.subset, fit$sample.subset]
  clabs <- colnames(fit$centroids)
  scen <- pamr.predict(fit, data$x, threshold = threshold, type = "cent")
  dif <- (scen - fit$centroid.overall)/fit$sd
  if(!is.null(fit$y)){
       nc <- length(unique(fit$y))
  }
   if(is.null(fit$y)){
      nc <- ncol(fit$proby)
}
  o <- drop(abs(dif) %*% rep(1, nc)) > 0
  d <- dif[o,  ]
  nd <- sum(o)
  genenames <- genenames[o]
  xx <- x[o,  ]
  oo <- order(apply(abs(d), 1, max))
  d <- d[oo,  ]
  genenames <- genenames[oo]
  par(mar = c(1, 5, 1, 1), col = 1)
  plot(rep(2, nd) + d[, 1], 1:nd, xlim = c(0, 2*nc+1), ylim = c(1, nd + 3), 
       type = "n", xlab = "", ylab = "", axes = FALSE)
  box()
  abline(h = seq(nd), lty = 3, col = 7)
  jj <- rep(0, nd)
  for(j in 1:nc) {
    segments(jj + 2 * j, seq(nd), jj + 2 * j + d[, j], seq(nd), col
             = j + 1, lwd = 4)
    lines(c(2 * j, 2 * j), c(1, nd), col = j + 1)
    text(2 * j, nd + 2, label = clabs[j], col = j + 1)
  }
  g <- substring(genenames, 1, 20)
  text(rep(0, nd), seq(nd), label = g, cex = 0.4, adj = 0, col = 1)
}
