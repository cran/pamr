pamr.listgenes <- function (fit, data, threshold, genenames = FALSE)  {
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
  print(res, quote = FALSE)
}

