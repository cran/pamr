#' A function to list the genes that survive the thresholding, from the nearest
#' shrunken centroid classifier
#' 
#' A function to list the genes that survive the thresholding, from the nearest
#' shrunken centroid classifier produced by pamr.train
#' 
#' \code{pamr.listgenes} List the geneids, and standardized centroids for each
#' class, for genes surviving at the given threshold. If fitcv is provided, the
#' function also reports the average rank of the gene in the cross-validation
#' folds, and the proportion of times that the gene is chosen (at the given
#' threshold) in the cross-validation folds.
#' 
#' @param fit The result of a call to pamr.train
#' @param data The input data.  In the same format as the input data for
#' pamr.train
#' @param threshold The desired threshold value
#' @param fitcv Optional object, result of a call to pamr.cv
#' @param genenames Include genenames in the list? If yes, they are taken from
#' "data". Default is false (geneid is always included in the list).
#' @author Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' 
#' #generate some data
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' 
#' mydata <- list(x=x,y=factor(y), geneid=as.character(1:nrow(x)),
#'                genenames=paste("g",as.character(1:nrow(x)),sep=""))
#' 
#' 
#' #train classifier
#' mytrain<-   pamr.train(mydata)
#' 
#' pamr.listgenes(mytrain, mydata, threshold=1.6)
#'  
#' 
#' @export pamr.listgenes
pamr.listgenes <- function (fit,   data, threshold, fitcv=NULL, genenames = FALSE)  {
  x <- data$x[fit$gene.subset, fit$sample.subset]
if (genenames) {
    gnames <- data$genenames[fit$gene.subset]
  }
  if (!genenames) {
    gnames <- NULL
  }
  geneid <- data$geneid[fit$gene.subset]
  if(!is.null(fit$y)){
       nc <- length(fit$y)
      }
 if(is.null(fit$y)){
       nc <- ncol(fit$proby)
      }
 clabs <- colnames(fit$centroids)

  aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
  cen <- pamr.predict(fit, x, threshold = threshold, type = "centroid")
  d <- (cen - fit$centroid.overall)[aa,, drop=FALSE]/fit$sd[aa]
  
  gene.order <- order(-apply(abs(d), 1, max))
  d <- round(d, 4)
  g <- gnames[aa]
  g1 <- geneid[aa]
  if (is.null(gnames)) {
    gnhdr <- NULL
  }
  if (!is.null(gnames)) {
    gnhdr <- "name"
  }

if(!is.null(fitcv)){
nfold=length(fitcv$cv.objects)

ind=matrix(F,nrow=nrow(x),ncol=nfold)
ranks=NULL
for( ii in 1:nfold){
	cen=pamr.predict(fitcv$cv.objects[[ii]], x[,-fitcv$folds[[ii]]],threshold=0, type="centroid")
	 dtemp <- (cen - fitcv$cv.objects[[ii]]$centroid.overall)[,, drop=FALSE]/fitcv$cv.objects[[ii]]$sd
	  r <- apply(abs(dtemp), 1, max)
	ranks=cbind(ranks,rank(-abs(r)))

	junk=pamr.predict(fitcv$cv.objects[[ii]], x[,-fitcv$folds[[ii]]],threshold=threshold, type="nonzero")
	ind[junk,ii]=T
}

av.rank=apply(ranks,1,mean)
av.rank=round(av.rank[aa],2)
prop=apply(ind[aa,,drop=F],1,sum)/nfold
}

  options(width = 500)
  schdr <- paste(clabs, "score", sep = "-")

if(is.null(fitcv)){
res <- cbind(as.character(g1), g, d)[gene.order,,drop=F]
  dimnames(res) <- list(NULL, c("id", gnhdr, schdr))

}
if(!is.null(fitcv)){
  res <- cbind(as.character(g1), g, d, av.rank, prop)[gene.order,,drop=F]
  dimnames(res) <- list(NULL, c("id", gnhdr, schdr, "av-rank-in-CV", "prop-selected-in-CV"))
}
  print(res, quote = FALSE)
}

