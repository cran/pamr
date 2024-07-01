#' A function to interactively define classes from a clustering tree
#' 
#' function to interactively define classes from a clustering tree
#' 
#' \code{pamr.makeclasses} Using this function the user interactively defines a
#' new set of classes, to be used in pamr.train, pamr.cv etc.  After invoking
#' pamr.makeclasses, a clustering tree is drawn.  This callss the R function
#' \code{hclust}, and any arguments for \code{hclust} can be passed to it.
#' Using the left button, the user clicks at the junction point defining the
#' subgroup 1. More groups can be added to class 1 by clicking on further
#' junction points. The user ends the definition of class 1 by clicking on the
#' rightmost button (in Windows, an additional menu appears and he chooses
#' Stop). This process is continued for classes 2,3 etc.  Note that some
#' sample may be left out of the new classes.  Two consecutive clicks of the
#' right button ends the definition for all classes.
#' 
#' At the end, the clustering is redrawn, with the new class labels shown.
#' 
#' Note: this function is "fragile". The user must click close to the junction
#' point, to avoid confusion with other junction points. Classes 1,2,3..
#' cannot have samples in common (if they do, an Error message will appear).
#' If the function is confused about the desired choices, it will complain and
#' ask the user to rerun pamr.makeclasses. The user should also check that the
#' labels on the final redrawn cluster tree agrees with the desired classes.
#' 
#' @param data The input data. A list with components: x- an expression genes
#' in the rows, samples in the columns, and y- a vector of the class labels for
#' each sample, and batchlabels- a vector of batch labels for each sample.
#' This object if the same form as that produced by pamr.from.excel.
#' @param sort.by.class Optional argument. If true, the clustering tree is
#' forced to put all samples in the same class (as defined by the class labels
#' y in `data') together in the tree. This is useful if a regrouping of classes
#' is desired. Eg: given classes 1,2,3,4 you want to define new classes (1,3)
#' vs (2,4) or 2 vs (1,3)
#' @param ... Any additional arguments to be passed to hclust
#' @return A vector of class labels 1,2,3...  If a component is NA (missing),
#' then the sample is not assigned to any class.  This vector should be
#' assigned to the newy component of data, for use in pamr.train etc.  Note
#' that pamr.train uses the class labels in the component ``newy'' if it is
#' present. Otherwise it uses the data labels ``y''.
#' @author Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' #generate some data
#' x <- matrix(rnorm(1000*40),ncol=40)
#' y <- sample(c(1:4),size=40,replace=TRUE)
#' batchlabels <- sample(c(1:5),size=40,replace=TRUE)
#' 
#' mydata <- list(x=x,y=factor(y),batchlabels=factor(batchlabels),
#'                geneid=as.character(1:nrow(x)),
#'                genenames=paste("g",as.character(1:nrow(x)),sep=""))
#' 
#' # mydata$newy <- pamr.makeclasses(mydata) Run this and define some new classes
#' 
#' train <- pamr.train(mydata)
#' 
#' @export pamr.makeclasses
pamr.makeclasses <- function(data,  sort.by.class = FALSE, ...) {
#  require(cluster)
  as.matrix.dist <- function (x)  {
    size <- attr(x, "Size")
    df <- matrix(0, size, size)
    df[row(df) > col(df)] <- x
    df <- df + t(df)
    labels <- attr(x, "Labels")
    dimnames(df) <- if (is.null(labels)) 
      list(1:size, 1:size)
    else list(labels, labels)
    df
  }
  as.dist <- function (m, diag = FALSE, upper = FALSE) {
    m <- as.matrix(m)
    retval <- m[row(m) > col(m)]
    attributes(retval) <- NULL
    if (!is.null(rownames(m))) 
      attr(retval, "Labels") <- rownames(m)
    else if (!is.null(colnames(m))) 
      attr(retval, "Labels") <- colnames(m)
    attr(retval, "Size") <- nrow(m)
    attr(retval, "Diag") <- diag
    attr(retval, "Upper") <- upper
    attr(retval, "call") <- match.call()
    class(retval) <- "dist"
    retval
  }
  
  if(!is.null(data$samplelabels)) {
    labs <- data$samplelabels
  }
  if(!is.null(data$samplelabels) & !is.null(data$y)) {
    labs <- paste(data$y, labs)
  }
  if(is.null(data$samplelabels)) {
    labs <- 1:ncol(data$x)
  }
  par(col = 1, cex = 1)
  d <- dist(t(data$x))
  dd <- as.matrix.dist(d)
  if(sort.by.class) {
    tt <- table(data$y)
    nc <- length(tt)
    for(i in 1:nc) {
      o <- data$y == names(tt[i])
      d1 <- max(dd[o, o])
      d2 <- min(dd[o, !o])
      fac <- ((0.2 + (0.7 * i)/nc) * d2)/d1
      dd[o, o] <- dd[o, o] * fac
    }
  }
  hc <- hclust(as.dist(dd), ...)
  plot(hc, labels = labs)
  aa <- vector("list", 100)
  go <- TRUE
  i <- 0
  while(go & i < 100) {
    go <- FALSE
    i <- i + 1
    print(c("Identify class", i))
    par(pch = as.character(i), col = 4)
    aa[[i]] <- locator(type = "p")
    if(!is.null(aa[[i]])) {
      go <- TRUE
    }
  }
  nclus <- i - 1
  res <- vector("list", nclus)
  for(i in 1:nclus) {
    res[i] <- aa[i]
  }
  hdelta <- 1
  clus <- vector("list", nclus)
  for(j in 1:nclus) {
    for(jj in 1:length(res[[j]]$x)) {
      r <- c(res[[j]]$x[jj], res[[j]]$y[jj])
      d <- abs(hc$hei - r[2])
      o <- rank(d)
      ncomp <- 5
      oo <- (1:length(o))[o < ncomp + 1 & d < hdelta]
      if(length(oo) == 0) {
        stop(
             "1 Ambigious selection; try pamr.makeclasses again"
             )
      }
      ncomp2 <- length(oo)
      good <- rep(FALSE, ncomp2)
      ordpos <- match(1:length(hc$ord), hc$ord)
      nodes <- vector("list", ncomp2)
      for(ii in 1:ncomp2) {
        ooo <- descendants(hc$mer, oo[ii])[[2]]
        o4 <- as.vector(hc$mer[ooo,  ])
        nodes[[ii]] <- -1 * o4[o4 < 0]
        op <- ordpos[nodes[[ii]]]
        if(r[1] > min(op) & r[1] < max(op)) {
          good[ii] <- TRUE
        }
      }
                                        #browser()
      if(sum(good) != 1) {
        stop(
             "2 Ambigious selection; try pamr.makeclasses again"
             )
      }
                                        #browser()
      ii2 <- (1:ncomp2)[good]
      clus[[j]] <- c(clus[[j]], nodes[[ii2]])
    }
  }
  newy <- rep(NA, ncol(data$x))
  temp <- NULL
  for(i in 1:nclus) {
    clus[[i]] <- unique(clus[[i]])
  }
  for(i in 1:nclus) {
    temp <- c(temp, clus[[i]])
  }
  if(length(unique(temp)) < length(temp)) {
    stop("Clusters overlap; try pamr.makeclasses again")
  }
  for(i in 1:nclus) {
    newy[clus[[i]]] <- i
  }
  labs2 <- as.character(newy)
  labs2[labs2 == "NA"] <- ""
  par(col = 1, cex = 1)
  plot(hc, labels = labs2)
  return(as.factor(newy))
}

