#' A function giving prediction information for many threshold values, from a
#' nearest shrunken centroid fit.
#' 
#' A function giving prediction information for many threshold values, from a
#' nearest shrunken centroid fit
#' 
#' 
#' @param fit The result of a call to pamr.train
#' @param newx Matrix of features at which predictions are to be made
#' @param threshold The desired threshold values
#' @param prior Prior probabilities for each class. Default is that specified
#' in "fit"
#' @param threshold.scale Additional scaling factors to be applied to the
#' thresholds. Vector of length equal to the number of classes.  Default is
#' that specified in "fit".
#' @param ... Additional arguments to be passed to pamr.predict
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
#' 
#' pamr.predictmany(mytrain, mydata$x)
#'  
#' 
#' @export pamr.predictmany
pamr.predictmany <- function(fit, newx, threshold=fit$threshold,
                             prior = fit$prior,  threshold.scale = fit$threshold.scale,
                             ...) {
  prob <-array(NA,c(length(prior),ncol(newx),length(threshold)))
  predclass <-matrix(NA,nrow=ncol(newx),ncol=length(threshold))
  
  for(i in 1:length(threshold)){
    prob[,,i] <-pamr.predict(fit,newx,threshold=threshold[i],type="posterior",...)
    predclass[,i] <-pamr.predict(fit,newx,threshold=threshold[i],type="class",...)
  }
  
  predclass <-matrix(levels(fit$y)[predclass],ncol=length((threshold)))

  return(list(prob=prob,predclass=predclass))
}













