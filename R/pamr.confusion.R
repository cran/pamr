pamr.confusion <- function(fit, threshold, extra = TRUE) {
  ii <- (1:length(fit$threshold))[fit$threshold > threshold]
  ii <- ii[1]
  predicted <- fit$yhat[, ii]
  if (is.null(fit$newy)) {
    true <- fit$y[fit$sample.subset]
  }
  else {
    true <- fit$newy[fit$sample.subset]
  }
  tt <- table(true, predicted)
  if (extra) {
    tt1 <- tt
    diag(tt1) <- 0
    tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
    dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
    print(tt)
    cat(c("Overall error rate=", round(sum(tt1)/sum(tt), 3)),
        fill= TRUE)
  }
  if (!extra) {
    return(tt)
  }
}
