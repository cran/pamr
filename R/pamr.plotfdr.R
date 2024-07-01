#' A function to plot the FDR curve from the nearest shrunken centroid
#' classifier
#' 
#' A function to plot the FDR curve the nearest shrunken centroid classifier
#' 
#' \code{pamr.plotfdr} plots the FDR curves from nearest shrunken centroid
#' classifier. The median FDR (solid line) and upper 90 percentile (broken
#' line) are shown
#' 
#' @param fdrfit The result of a call to pamr.fdr
#' @param call.win.metafile Used by Excel interface
#' @author Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:2),size=20,replace=TRUE)
#' x[1:50,y==2]=x[1:50,y==2]+3
#' mydata <- list(x=x,y=y)
#' mytrain <-   pamr.train(mydata)
#' myfdr <-  pamr.fdr(mytrain, mydata)
#' pamr.plotfdr(myfdr)
#' 
#' @export pamr.plotfdr
pamr.plotfdr <- function(fdrfit,  call.win.metafile=FALSE){

##if(call.win.metafile){win.metafile()}

om=fdrfit$results[,"Number of significant genes"]>0

na.min=function(x){min(x[!is.na(x)])}
na.max=function(x){max(x[!is.na(x)])}
plot(fdrfit$results[om,"Number of significant genes"],fdrfit$results[om,"Median FDR"],log="x",
xlab="Number of genes called significant",
ylab="False discovery rate (median and 90th percentile)",type="b",  
ylim=c(na.min(fdrfit$results[om,"Median FDR"]), na.max(fdrfit$results[om,"90th percentile of FDR"])))
x=fdrfit$results[om,"Number of significant genes"]
xlim <- range(x)
barw <- abs((log(x)))*1.2
upper=fdrfit$results[om,"90th percentile of FDR"]
lower=fdrfit$results[om,"Median FDR"]
segments(x, upper, x, lower, lty=2)
segments(x - barw, upper, x + barw, upper, lty=2)

axis(3,at=fdrfit$results[om,"Number of significant genes"], labels=round(fdrfit$results[om,"Threshold"],2))

mtext("Threshold", 3, 2, cex = 1.0)

if (call.win.metafile) {
    savePlot("", type="wmf")
}
dev.off()

return()
}
