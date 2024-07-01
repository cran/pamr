#' A function to plot the genes that surive the thresholding from the nearest
#' shrunken centroid classifier
#' 
#' A function to plot the genes that survive the thresholding, from the nearest
#' shrunken centroid classifier produced by pamr.train
#' 
#' \code{pamr.geneplot} Plots the raw gene expression for genes that survive
#' the specified threshold. Plot is stratified by class.  Plot is set up to
#' display only up to about 20 or 25 genes, otherwise it gets too crowded.
#' Hence threshold should be chosen to yield at most about 20 or 25 genes.
#' 
#' @param fit The result of a call to pamr.train
#' @param data The input data.  In the same format as the input data for
#' pamr.train
#' @param threshold The desired threshold value
#' @author Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' mydata <- list(x=x,y=y)
#' mytrain <-   pamr.train(mydata)
#' pamr.geneplot(mytrain, mydata, threshold=1.6)
#'  
#' 
#' @export pamr.geneplot
pamr.geneplot <- function(fit, data, threshold) {
  par(pch = 1, col = 1)
  geneid <- data$geneid
  if(is.null(geneid)) {
    geneid <- as.character(1:nrow(data$x))
  }
  if(is.null(fit$newy)) {
    y <- factor(data$y[fit$sample.subset])
  }
  else {
    y <- factor(fit$newy[fit$sample.subset])
  }
  x <- data$x[fit$gene.subset, fit$sample.subset]
  geneid <- geneid[fit$gene.subset]
  nc <- length(unique(y))
  aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
  cen <- pamr.predict(fit, x, threshold = threshold, type = "cen")
  d <- (cen - fit$centroid.overall)[aa,  ]/fit$sd[aa]
  oo <- order( - apply(abs(d), 1, max))
  aa <- aa[oo]
  ngenes <- length(aa)
  o <- order(y)
  xx <- x[aa, o]
  geneid <- geneid[aa]
  nc <- length(unique(y))
  nn <- c(0, cumsum(table(y)))
  nrow <- trunc(sqrt(ngenes)) + 1
  ncol <- trunc(sqrt(ngenes)) + 1
  if(nrow * (ncol - 1) >= ngenes) {
    ncol <- ncol - 1
  }
  par(mfrow = c(nrow, ncol))
  for(i in 1:ngenes) {
    plot(1:ncol(xx), xx[i,  ], type = "n", xlab = "sample", ylab = 
         "expression", axes = FALSE)
    box()
    axis(2)
    for(j in 1:nc) {
      j1 <- nn[j] + 1
      j2 <- nn[j] + table(y)[j]
      points(j1:j2, xx[i, j1:j2], col = j + 1)
    }
    title(main = as.character(geneid[i]))
    for(j in 1:(nc - 1)) {
      abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
    }
    if(i == 1) {
      h <- c(0, table(y))
      for(j in 2:(nc + 1)) {
        text(sum(h[1:(j - 1)]) + 0.5 * h[j], max(xx[i,  
                                                    ]), label = levels(y)[j - 1], col = j)
      }
    }
  }
  par(mfrow = c(1, 1))
}


