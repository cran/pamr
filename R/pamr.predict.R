#' A function giving prediction information, from a nearest shrunken centroid
#' fit.
#' 
#' A function giving prediction information, from a nearest shrunken centroid
#' fit
#' 
#' \code{pamr.predict} Give a cross-tabulation of true versus predicted classes
#' for the fit returned by pamr.train or pamr.cv, at the specified threshold
#' 
#' @param fit The result of a call to pamr.train
#' @param newx Matrix of features at which predictions are to be made
#' @param threshold The desired threshold value
#' @param type Type of prediction desired: class predictions, posterior
#' probabilities, (unshrunken) class centroids, vector of genes surviving the
#' threshold
#' @param prior Prior probabilities for each class. Default is that specified
#' in "fit"
#' @param threshold.scale Additional scaling factors to be applied to the
#' thresholds. Vector of length equal to the number of classes.  Default is
#' that specified in "fit".
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
#' mycv <- pamr.cv(mytrain,mydata)
#' pamr.predict(mytrain, mydata$x , threshold=1)
#'  
#' 
#' @export pamr.predict
pamr.predict <-  function(fit, newx, threshold, type = c("class", "posterior", "centroid", "nonzero"), 
                          prior = fit$prior,  threshold.scale = fit$
                          threshold.scale) {
  norm.cen <- fit$norm.cen
  if(!is.null(norm.cen)) {
    newx <- abs(t(scale(t(newx), center = norm.cen, scale = FALSE)))
  }
  type <- match.arg(type)
  sd <- fit$sd
  centroid.overall <- fit$centroid.overall
  centroids <- fit$centroids
  se.scale <- fit$se.scale
  delta <- scale((centroids - centroid.overall)/sd, FALSE, threshold.scale * 
                 se.scale)

  if(fit$sign.contrast=="positive"){delta <- delta*(delta>0)}
  if(fit$sign.contrast=="negative"){delta <- delta*(delta<0)}


  delta.shrunk <- scale(soft.shrink(delta, threshold), FALSE,
                        1/(  threshold.scale * se.scale))
  posid <- drop(abs(delta.shrunk) %*% rep(1, length(prior))) > 0
                
  if(!match(type, c("centroid", "nonzero"), FALSE))
    dd <- diag.disc((newx - centroid.overall)/sd, delta.shrunk, 
                    prior, posid)
  switch(type,
         class = softmax(dd),
         posterior = {
           dd <- safe.exp(dd)
           dd/drop(dd %*% rep(1, length(prior)))
         }
         ,
         centroid = centroid.overall + delta.shrunk * sd,
         nonzero = {
           nz <- drop(abs(delta.shrunk) %*% rep(1, ncol(centroids)
                                                )) > 0
           seq(nz)[nz]
         }
         )
}

safe.exp=function(x){
 xx=sign(x)*pmin(abs(x),500)
 return(exp(xx))
}

