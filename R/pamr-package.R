#' Prediction Analysis of Microarrays
#'
#' @description Functions for training and predicting using shrunken
#'   centroids. While originally developed in the context of
#'   microarray data analysis (<doi:10.1073/pnas.082099299>), the
#'   method is broadly applicable.
#' @aliases pamr.cube.root pamr.pairscore pamr.pvalue.survival
#'   pamr.score.to.class1 pamr.score.to.class2 print.nsc print.nsccv
#'   print.pamrcved print.pamrtrained
#' @import survival cluster
#' @importFrom grDevices dev.off savePlot
#' @importFrom graphics abline axis box legend lines locator mtext par plot points segments text title
#' @importFrom stats approx dist hclust kmeans lm median model.matrix pchisq predict quantile runif update var
#' @importFrom utils  menu
#' @keywords internal
"_PACKAGE"



