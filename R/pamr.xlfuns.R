pamr.xl.compute.offset <- function(data, offset.percent=50, prior=prior){
  x <- data$x
  y <- data$y
  n.class <- table(y)
  if(min(n.class)==1){stop("Error: each class must have >1 sample")}
  norm.cent <-NULL
  n <- sum(n.class)
  xtest <- x
  ntest <- ncol(xtest)
  K <- length(prior)
  p <- nrow(x)
  Y <- model.matrix( ~ factor(y) - 1, data = list(y = y))
  dimnames(Y) <- list(NULL, names(n.class))
  centroids <- scale(x %*% Y, FALSE, n.class)
  sd <- rep(1, p)
  xdif <- x - centroids %*% t(Y)
  sd <- (xdif^2) %*% rep(1/(n - K), n)
  sd <- drop(sqrt(sd))
  offset  <- quantile(sd, offset.percent/100)
  return(offset)
}

pamr.xl.get.offset  <- function() {
  if (exists("x.train")) {
    return (x.train$offset)
  } else {
    return (pamr.xl.compute.offset(pamr.xl.data, offset.percent=pamr.xl.training.parameters$offset.percent,
                                   prior=pamr.xl.training.parameters$prior))
  }
}

pamr.xl.derive.adjusted.prior  <- function(prior, data) {
  s  <- pamr.xl.get.sample.prior(data)
  temp <- prior - s
  if (sum(temp*temp) < pamr.xl.training.parameters$epsilon) {
    return (list (prior=s, prior.name="Sample Prior"))
  } else {
    s  <-  pamr.xl.get.uniform.prior(data)
    temp  <- prior - s
    if (sum(temp*temp) < pamr.xl.training.parameters$epsilon) {
      return (list (prior=s, prior.name="Uniform Prior"))
    } else {
      return (list (prior=prior, prior.name="Custom Prior"))      
    }
  }
}

pamr.xl.get.default.training.parameters <- function(data) {
  return (list(offset.percent=50, prior=pamr.xl.get.sample.prior(data), prior.name="Sample Prior", sign.contrast="both", epsilon=1e-7))
}

## Return the uniform prior on class labels
pamr.xl.get.uniform.prior  <- function(data) {
  w <- table(data$y)
  n  <- length(w)
  return(rep(1.0/n, n))
}

## Return the sample proportion prior on class labels
pamr.xl.get.sample.prior  <- function(data) {
  w <- table(data$y)
  return(w/sum(w))
}


pamr.xl.process.data <- function() {
  res <- list(x=pamr.xl.raw.data, y=pamr.xl.class.labels, genenames=pamr.xl.gene.names, 
              geneid=pamr.xl.gene.ids, samplelabels=pamr.xl.sample.labels,
              batchlabels=pamr.xl.batch.labels)
  
  if (pamr.xl.data.has.missing.values) {
    res <- pamr.knnimpute(res, k = pamr.xl.knn.neighbors)
  }
  return(res)
}

pamr.xl.compute.cv.confusion  <- function (fit, cv.results, threshold) {
  threshold.rank  <- which(rank(abs(cv.results$threshold - threshold))==1)
  t.threshold  <- cv.results$threshold[threshold.rank]
  true  <- fit$y
  predicted  <- cv.results$yhat[, threshold.rank]
  tt <- table(true, predicted)
  tt1 <- tt
   diag(tt1) <- 0
  tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
  dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
  overall.err  <- round(sum(tt1)/sum(tt), 3)
  return(list(confusion.matrix=tt, overall.error=overall.err, threshold=round(t.threshold, 5)))
 }
pamr.xl.compute.confusion  <- function (fit, threshold) {
  ii <- (1:length(fit$threshold))[fit$threshold > threshold]
  ii <- ii[1]
  predicted <- fit$yhat[, ii]
  if(is.null(fit$newy)) {
    true <- fit$y[fit$sample.subset]
  }
  else {
    true <- fit$newy[fit$sample.subset]
  }
  tt <- table(true, predicted)
  tt1 <- tt
  diag(tt1) <- 0
  tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
  dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
  overall.err  <- round(sum(tt1)/sum(tt), 3)
  return(list(confusion.matrix=tt, overall.error=overall.err))
}
pamr.xl.is.a.subset  <- function(x, y) {
  if (nlevels(factor(x)) == nlevels(factor(c(x, y[!is.na(y)])))) {
    return (1)  # True
  } else {
    return (0)  # False
  }
}
pamr.xl.listgenes.compute  <- function (fit, data, threshold, genenames = FALSE) {
  if (is.null(fit$newy)) {
    y <- factor(data$y[fit$sample.subset])
  }
  if (!is.null(fit$newy)) {
    y <- factor(fit$newy[fit$sample.subset])
  }
  x <- data$x[fit$gene.subset, fit$sample.subset]
  if (genenames) {
    gnames <- data$genenames[fit$gene.subset]
  }
  if (!genenames) {
    gnames <- NULL
  }
  geneid <- data$geneid[fit$gene.subset]
  nc <- length(unique(y))
  aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
  cen <- pamr.predict(fit, x, threshold = threshold, type = "cen")
  d <- (cen - fit$centroid.overall)[aa, ]/fit$sd[aa]
  oo <- order(-apply(abs(d), 1, max))
  d <- round(d, 4)
  g <- gnames[aa]
  g1 <- geneid[aa]
  if (is.null(gnames)) {
    gnhdr <- NULL
  }
  if (!is.null(gnames)) {
    gnhdr <- "name"
  }
  options(width = 500)
  schdr <- paste(dimnames(table(y))$y, "score", sep = " ")
  res <- cbind(as.character(g1), g, d)[oo, ]
  dimnames(res) <- list(NULL, c("id", gnhdr, schdr))
  return(list(gene.headings = dimnames(res)[[2]],
              gene.ids = res[ , 2],   # This was switched with gene.names. 
              gene.names = res[ , 1],
              gene.scores = res[ , -(1:2)]))
  ##print(res, quote = FALSE)
}
pamr.xl.plot.test.probs.compute  <- function(fit, new.x, newx.classes, missing.class.label, 
	threshold, sample.labels=NULL) {
  predicted.probs  <- pamr.xl.predict.test.probs(x.train, new.x, threshold=threshold)
  training.classes  <- levels(factor(fit$y))
  py  <- pamr.xl.predict.test.class.only(x.train, new.x, threshold=threshold)
  order.classes  <- order(newx.classes)
  pp  <- predicted.probs[, order.classes]
  actual.classes <- newx.classes[order.classes]
  actual.classes[is.na(actual.classes)] <- missing.class.label
  ny  <- py$predicted[order.classes]
  n  <- length(ny)
  ss  <- sample.labels
  if (!is.null(ss)) {
    ss  <- ss[order.classes]
  }
  
  return (list(x = 1:n,
               y = t(pp),
               x.label = "Sample",
               y.label = "Predicted Test Probabilities",
               y.names = levels(factor(fit$y)),
               y.lines = cumsum(table(actual.classes)) + 0.5,
               x.dummy = vector(length=2, mode="numeric"),
               y.dummy = vector(length=2, mode="numeric"),
               panel.names = levels(factor(actual.classes)),
               x.names = ss))
}  



pamr.xl.plot.training.error.compute  <- function(trained.object) {
  return (list(x = trained.object$threshold,
               y = trained.object$errors/length(trained.object$y),
               y.ytop = trained.object$nonzero,
               x.label = "Threshold",
               y.label = "Training Error"))
}
pamr.xl.plotcen.compute  <- function(fit, data, threshold) {
  genenames <- data$genenames[fit$gene.subset]
  x <- data$x[fit$gene.subset, fit$sample.subset]
  clabs <- colnames(fit$centroids)
  scen <- pamr.predict(fit, data$x, threshold = threshold, type = "cent")
  dif <- scen - fit$centroid.overall
  nc <- length(unique(fit$y))
  o <- drop(abs(dif) %*% rep(1, nc)) > 0
  d <- dif[o,  ]
  nd <- sum(o)
  genenames <- genenames[o]
  xx <- x[o,  ]
  oo <- order(apply(abs(d), 1, max))
  d <- d[oo,  ]
  genenames <- genenames[oo]
  win.metafile()
  plot.title=paste("Centroid Plot( Threshold =", threshold, ")")
  par(mar = c(1, 5, 1, 1), col = 1)
  plot(rep(2, nd) + d[, 1], 1:nd, xlim = c(0, 2*nc+1), ylim = c(1, nd + 3), 
       type = "n", xlab = "", ylab = "", axes = FALSE, main=plot.title)
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
  dev.off()
#  pamr.plot.y <<- matrix(d, nrow=dim(d)[1])
#  pamr.plot.x <<- seq(nd)
#  pamr.plot.seriesnames <<- dimnames(d)[[2]]
#  pamr.plot.genenames <<- genenames

  return(TRUE)
}
pamr.xl.plotcv.compute  <- function(aa) {
  n <- nrow(aa$yhat)
  y <- aa$y
  if(!is.null(aa$newy)) {
    y <- aa$newy[aa$sample.subset]
  }
  nc <- length(table(y))
  nfolds <- length(aa$folds)
  err <- matrix(NA, ncol = ncol(aa$yhat), nrow = nfolds)
  temp <- matrix(y, ncol = ncol(aa$yhat), nrow = n)
  ni <- rep(NA, nfolds)
  for(i in 1:nfolds) {
    ii <- aa$folds[[i]]
    ni[i] <- length(aa$folds[[i]])
    err[i,  ] <- apply(temp[ii,  ] != aa$yhat[ii,  ], 2, sum)/ni[i]
  }
  se <- sqrt(apply(err, 2, var)/nfolds)

  err2 <- matrix(NA, nrow = length(unique(y)), ncol = length(aa$threshold)-1)
  for(i in 1:(length(aa$threshold) - 1)) {
    s <- pamr.confusion(aa, aa$threshold[i], extra = FALSE)
    diag(s) <- 0
    err2[, i] <- apply(s, 1, sum)/table(y)
  }

  return (list(x = aa$threshold,
               y = aa$error,
               x.label = "Threshold",
               y.label = "Misclasiffication Error",
               y.se = se,
               y.ytop = aa$size,
               cv.err = t(err2),
               cv.legend = dimnames(table(y))[[1]]))
               
}
pamr.xl.plotcvprob.compute  <- function(aa, data, threshold) {
  ii <- (1:length(aa$threshold))[aa$threshold > threshold]
  ii <- ii[1]
  ss <- data$samplelabels
  pp <- aa$prob[,  , ii]
  if(is.null(aa$newy)) {
    y <- aa$y[aa$sample.subset]
  }
  if(!is.null(aa$newy)) {
    y <- aa$newy[aa$sample.subset]
  }
  o <- order(y)
  y <- y[o]
  if(!is.null(ss)) {
    ss <- ss[o]
  }
  ppp <- pp[o,  ]
  n <- nrow(ppp)
  nc <- length(unique(y))


#  axis(2, labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0", ""))
#  if (!is.null(ss)) {
#    pamr.plot.x.names <<- ss
#  }

  return (list(x = 1:n,
               y = ppp,
               x.label = "Sample",
               y.label = "CV Probabilities",
               y.names = levels(y),
               y.lines = cumsum(table(data$y)),
               x.dummy = vector(length=2, mode="numeric"),
               y.dummy = vector(length=2, mode="numeric"),
               x.names = ss))
  
#   for(j in 1:nc) {
#     points(1:n, ppp[, j], col = j + 1)
#   }
#   for(j in 1:(nc - 1)) {
#     abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
#   }
#   h <- c(0, table(y))
#   for(j in 2:(nc + 1)) {
#     text(sum(h[1:(j - 1)]) + 0.5 * h[j], 1.02, label = levels(y)[j - 
#                                                  1], col = j)
#   }
#   abline(h = 1)
#   if(!is.null(ss)) {
#     text(1:length(ss), 1.1, labels = ss, srt = 90, cex = 0.7)
#   }
  ##if(!is.null(ss)){axis(3,labels=ss,at=1:length(ss),srt=90)}
}
pamr.xl.predict.test.class<- function(fit, newx, threshold, test.class.labels) {
  predicted  <- pamr.predict(fit, newx, threshold, type="class")
  return(list(confusion.matrix=table(test.class.labels, predicted), predicted=as.vector(predicted)))
}

pamr.xl.predict.test.class.only  <- function(fit, newx, threshold) {
  return(list(predicted=as.vector(pamr.predict(fit, newx, threshold, type="class"))))
}
pamr.xl.predict.test.class.only  <- function(fit, newx, threshold) {
  return(list(predicted=as.vector(pamr.predict(fit, newx, threshold, type="class"))))
}
pamr.xl.predict.test.probs  <- function(fit, newx, threshold) {
  predicted  <- pamr.predict(fit, newx, threshold, type="posterior")
  return(t(predicted))
}

pamr.xl.test.data.impute  <- function(x, k) {
  N <- dim(x)
  p <- N[2]
  N <- N[1]
  col.nas  <- apply(x, 2, is.na)
  if ((sum(col.nas) == N) > 0) {
    stop("Error: A column has all missing values!")
  }
  
  nas <- is.na(drop(x %*% rep(1, p)))
  xcomplete <- x[!nas,  ]
  xbad <- x[nas,,drop=FALSE ]
  xnas <- is.na(xbad)
  xbadhat <- xbad
  cat(nrow(xbad), fill = TRUE)
  for(i in seq(nrow(xbad))) {
    cat(i, fill = TRUE)
    xinas <- xnas[i,  ]
    xbadhat[i,  ] <- nnmiss(xcomplete, xbad[i,  ], xinas, K = k)
  }
  x[nas,  ] <- xbadhat
  return(x)
}
pamr.xl.test.errors.compute  <- function(fit, newx, newx.classes, threshold=fit$threshold,
                                         prior = fit$prior,  threshold.scale = fit$threshold.scale,
                                         ...) {
  n  <- length(which(!is.na(newx.classes)))
## Note: n is assumed to be nonzero! Check before calling!
  actual.classes  <- newx.classes
  prediction.errs  <- vector(mode="numeric", length=length(threshold))
  
  for(i in 1:length(threshold)){
    t <- pamr.predict(fit,newx,threshold=threshold[i],type="class",...)
    prediction.errs[i]  <- length(which(t != actual.classes)) / n
  }
  
  return(list(x=threshold, y=prediction.errs, x.label="Threshold", y.label="Test Error", ))
  
}
pamr.xl.transform.class.labels  <- function(x) {
  y  <- x
  y[is.na(y)]  <- " "
  return(y)
}

pamr.xl.transform.data <- function(data) {

  if (pamr.xl.take.cube.root) {
    data$x = pamr.cube.root(data$x)
  }

  if (pamr.xl.batch.labels.present) {
    data <- pamr.batchadjust(data)
  }

  if (pamr.xl.center.columns && pamr.xl.scale.columns) {
    data$x = scale(data$x, center=TRUE, scale=TRUE)
  } else if (pamr.xl.center.columns) {
    data$x = scale(data$x, center=TRUE, scale=FALSE)
  } else if (pamr.xl.scale.columns) {
    data$x = scale(data$x, center=FALSE, scale=TRUE)
  }

  return (data)
}

pamr.xl.transform.test.data <- function(test.x) {
  res <- test.x
  if (pamr.xl.take.cube.root) {
    res = pamr.cube.root(res)
  }

  if (pamr.xl.center.columns && pamr.xl.scale.columns) {
    res = scale(res, center=TRUE, scale=TRUE)
  } else if (pamr.xl.center.columns) {
    res = scale(res, center=TRUE, scale=FALSE)
  } else if (pamr.xl.scale.columns) {
    res = scale(res, center=FALSE, scale=TRUE)
  }

  return (res)
}

