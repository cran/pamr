#' A function to decorrelate (adjust) the feature matrix with respect to some
#' additional predictors
#' 
#' A function to decorrelate (adjust) the feature matrix with respect to some
#' additional predictors
#' 
#' \code{pamr.decorrelate} Does a least squares regression of each row of x on
#' the adjusting predictors, and returns the residuals. If xtest is provided,
#' it also returns the adjusted version of xtest, using the training set least
#' squares regression model for adjustment
#' 
#' @param x Matrix of training set feature values, with genes in the rows,
#' samples in the columns
#' @param adjusting.predictors List of training set predictors to be used for
#' adjustment
#' @param xtest Optional matrix of test set feature values, to be adjusted in
#' the same way as the training set
#' @param adjusting.predictors.test Optional list of test set predictors to be
#' used for adjustment
#' @return A list with components \item{x.adj}{Adjusted x matrix}
#' \item{xtest.adj}{Adjusted xtest matrix, if xtest we provided}
#' @author Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @references Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan,
#' and Gilbert Chu Diagnosis of multiple cancer types by shrunken centroids of
#' gene expression PNAS 99: 6567-6572.  Available at www.pnas.org
#' @examples
#' 
#' #generate some data
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' 
#' x<-matrix(rnorm(1000*20),ncol=20)
#' y<-c(rep(1,10),rep(2,10))
#' adjusting.predictors=list(pred1=rnorm(20), pred2=as.factor(sample(c(1,2),replace
#' =TRUE,size=20)))
#' xtest=matrix(rnorm(1000*10),ncol=10)
#' adjusting.predictors.test=list(pred1=rnorm(10), pred2=as.factor(sample(c(1,2),replace
#' =TRUE,size=10)))
#' 
#' # decorrelate training x wrt adjusting predictors
#' 
#' x.adj=pamr.decorrelate(x,adjusting.predictors)$x.adj
#' # train classifier with adjusted x
#' 
#' d=list(x=x.adj,y=y)
#' a<-pamr.train(d)
#' 
#' # decorrelate training and test x wrt adjusting predictors, then make
#' #predictions for test set
#' 
#' temp <- pamr.decorrelate(x,adjusting.predictors, xtest=xtest,
#'                          adjusting.predictors.test=adjusting.predictors.test)
#' 
#' d=list(x=temp$x.adj,y=y)
#' a<-pamr.train(d)
#' aa<-pamr.predict(a,temp$xtest.adj, threshold=.5)
#' 
#' @export pamr.decorrelate
pamr.decorrelate<-  function (x, adjusting.predictors, xtest=NULL, adjusting.predictors.test=NULL){

foo<- lm(t(x)~., adjusting.predictors)
x.adj=t(foo$res)
xtest.adj=NULL

if(!is.null(adjusting.predictors.test)){
   if(is.null(xtest)){
   stop("xtest must be non-null if adjusting.predictors.test is non-null")
  }
    temp=t(predict(foo,adjusting.predictors.test))
    xtest.adj=xtest-temp
}
return(list(x.adj=x.adj,xtest.adj=xtest.adj))
}
