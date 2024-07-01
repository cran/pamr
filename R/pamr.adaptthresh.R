#' A function to adaptive choose threshold scales, for use in pamr.train
#' 
#' A function to adaptive choose threshold scales, for use in pamr.train
#' 
#' \code{pamr.adaptthresh} Adaptively searches for set of good threshold
#' scales.  The baseline (default) scale is 1 for each class. The idea is that
#' for easy to classify classes, the threshold scale can be increased without
#' increasing the error rate for that class, and resulting in fewer genes
#' needed for the classification rule. The scalings from pamr.adaptthresh are
#' then used in pamr.train, and pamr.cv. The results may be better than those
#' obtained with the default values of threshold.scale.
#' 
#' @param object The result of a call to pamr.train
#' @param ntries Number of iterations to use in algorithm
#' @param reduction.factor Amount by which a scaling is reduced in one step of
#' the algorithm
#' @param full.out Should full output be returned? Default FALSE
#' @author Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @references
#' 
#' Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert
#' Chu. "Diagnosis of multiple cancer types by shrunken centroids of gene
#' expression" PNAS 2002 99:6567-6572 (May 14).
#' 
#' Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert
#' Chu (2002).  Class prediction by nearest shrunken centroids,with
#' applications to DNA microarrays. Stanford tech report.
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' mydata <- list(x=x,y=y)
#' mytrain <-   pamr.train(mydata)
#' new.scales <- pamr.adaptthresh(mytrain)
#' 
#'  
#' mytrain2 <- pamr.train(mydata, threshold.scale=new.scales)
#' 
#' myresults2 <- pamr.cv(mytrain2, mydata)
#' 
#' 
#' @export pamr.adaptthresh
pamr.adaptthresh <- function(object, ntries = 10, reduction.factor = 0.9, full.out = FALSE) {
  errors <- error.nsc(object)
  threshold <- object$threshold   
### Remove all but the first leading zero errors
  ifirst <- match(TRUE, object$errors > 0, FALSE)
  if (!ifirst)
    stop("Zero training error throughout!")
  else {
    ifirst <- max(ifirst, 1)
    threshold <- threshold[seq(ifirst, length(threshold))]
  }
### initialization
  tscales <- object$threshold.scale
  all.errors <- matrix(0, ntries + 1, length(tscales),
                       dimnames = list(NULL, names(tscales)))
  all.scales <- all.errors
  all.objects <- as.list(seq(ntries + 1))
  rocs <- double(ntries + 1)
  all.scales[1,  ] <- tscales
  all.errors[1,  ] <- errors
  rocs[1] <- roc.nsc(object)      # integrated size^(1/4)*error
  cat("Initial errors:", format(round(errors, 5)), "Roc",
      format(round(rocs[1], 5)), "\n")
  for (i in seq(ntries)) {
    cat("Update", i, "\n")
    j <- rev(order(errors))[1]      # identify the largest error
    tscales[j] <- tscales[j] * reduction.factor     
                                        # and reduce its scale
    all.scales[i + 1,  ] <- tscales/min(tscales)    # and renormalize
    iobject <- update(object, threshold = threshold, 
                      threshold.scale = all.scales[i + 1,  ], remove.zeros = 
                      FALSE)
    all.errors[i + 1,  ] <- errors <- error.nsc(iobject)
    rocs[i + 1] <- roc.nsc(iobject)
    cat("\nErrors", format(round(errors, 5)), "Roc",
        format(round(rocs[i + 1], 5)), "\n")
  }
  j <- order(rocs)[1]     # identify the scales with the smallest "roc"
  opt.scale <- all.scales[j,  ]
  if (full.out)
    list(errors = all.errors, scales = all.scales, rocs = rocs, 
         opt.scale = opt.scale)
  else
    opt.scale
}



#' A function to mean-adjust microarray data by batches
#' 
#' A function to mean-adjust microarray data by batches
#' 
#' \code{pamr.batchadjust} does a genewise one-way ANOVA adjustment for
#' expression values.  Let \eqn{x(i,j)} be the expression for gene \eqn{i} in sample \eqn{j}.
#' Suppose sample \eqn{j} in in batch \eqn{b}, and let \eqn{B} be the set of all samples in batch
#' \eqn{b}. Then \code{pamr.batchadjust} adjusts \eqn{x(i,j)} to \eqn{x(i,j) - mean[x(i,j)]}
#' where the mean is taken over all samples \eqn{j} in \eqn{B}.
#' 
#' @param data The input data. A list with components: x- an expression genes
#' in the rows, samples in the columns, and y- a vector of the class labels for
#' each sample, and batchlabels- a vector of batch labels for each sample.
#' @return A data object of the same form as the input data, with x replaced by
#' the adjusted x
#' @author Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' #generate some data
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' batchlabels <- sample(c(1:5),size=20,replace=TRUE)
#' mydata <- list(x=x,y=factor(y),batchlabels=factor(batchlabels))
#' 
#' mydata2 <- pamr.batchadjust(mydata)
#' 
#' @export pamr.batchadjust
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
