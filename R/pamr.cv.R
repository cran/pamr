#' A function to cross-validate the nearest shrunken centroid classifier
#' 
#' A function to cross-validate the nearest shrunken centroid classifier
#' produced by pamr.train
#' 
#' \code{pamr.cv} carries out cross-validation for a nearest shrunken centroid
#' classifier.
#' 
#' @param fit The result of a call to pamr.train
#' @param data A list with at least two components: x- an expression genes in
#' the rows, samples in the columns), and y- a vector of the class labels for
#' each sample. Same form as data object used by pamr.train.
#' @param nfold Number of cross-validation folds. Default is the smallest class
#' size
#' @param folds A list with nfold components, each component a vector of
#' indices of the samples in that fold. By default a (random) balanced
#' cross-validation is used
#' @param ... Any additional arguments that are to be passed to pamr.train
#' @return A list with components \item{threshold}{A vector of the thresholds
#' tried in the shrinkage} \item{errors}{The number of cross-validation errors
#' for each threshold value} \item{loglik}{The cross-validated multinomial
#' log-likelihood value for each threshold value} \item{size}{A vector of the
#' number of genes that survived the thresholding, for each threshold value
#' tried.}.  \item{yhat}{A matrix of size n by nthreshold, containing the
#' cross-validated class predictions for each threshold value, in each column}
#' \item{prob}{A matrix of size n by nthreshold, containing the cross-validated
#' class probabilities for each threshold value, in each column} \item{folds}{
#' The cross-validation folds used} \item{cv.objects}{Train objects (output of
#' pamr.train), from each of the CV folds} \item{call}{The calling sequence
#' used}
#' @author Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' 
#' mydata <- list(x=x,y=factor(y), geneid=as.character(1:nrow(x)),
#'  genenames=paste("g",as.character(1:nrow(x)),sep=""))
#' 
#' mytrain <-   pamr.train(mydata)
#' mycv <- pamr.cv(mytrain,mydata)
#' 
#' @export pamr.cv
pamr.cv <-
function(fit, data, nfold = NULL, folds = NULL ,...)
{
        x <- data$x[fit$gene.subset, fit$sample.subset]

        if( !is.null(data$y) & !is.null(data$proby)){
           stop("Must have exactly one of y and  proby  present in the data object")
         }
        
        y <- NULL
        proby <- NULL
        
        if(!is.null(fit$y)){
           y<-  factor(fit$y[fit$sample.subset])
         }
        
        if(!is.null(fit$proby)){
           proby<-  fit$proby[fit$sample.subset,]
         }
        
        this.call <- match.call()
        
# three possibilities, 
# problem.type= class: y are class labels, proby=NULL
#               surv.km: y=NULL, proby are soft class probabilities from KM
#               surv.latent: y are latent class labels,
#                   proby are soft class probabilities from KM
# note; problem type is in fit$problem.type
        
        junk <- nsccv(x, y=y, proby=proby, object = fit, nfold=nfold, folds=folds, 
survival.time=data$survival.time, censoring.status = data$censoring.status, 
ngroup.survival=fit$ngroup.survival, problem.type=fit$problem.type, ...)

        junk$call <- this.call
        
        junk$sample.subset <- fit$sample.subset
        class(junk)="pamrcved"
        junk
}

