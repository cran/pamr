#' A function to plot the survival curves in each Kaplan Meier stratum
#' 
#' A function to plot the survival curves in each Kaplan Meier stratum
#' 
#' 
#' @param fit The result of a call to pamr.train
#' @param survival.time Vector of survival times
#' @param censoring.status Vector of censoring status values
#' @author Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' 
#' gendata<-function(n=100, p=2000){
#'   tim <- 3*abs(rnorm(n))
#'   u<-runif(n,min(tim),max(tim))
#'   y<-pmin(tim,u)
#'    ic<-1*(tim<u)
#' m <- median(tim)
#' x<-matrix(rnorm(p*n),ncol=n)
#'   x[1:100, tim>m] <-  x[1:100, tim>m]+3
#'   return(list(x=x,y=y,ic=ic))
#' }
#' 
#' # generate training data; 2000 genes, 100 samples
#' 
#' junk<-gendata(n=100)
#' y<-junk$y
#' ic<-junk$ic
#' x<-junk$x
#' d <- list(x=x,survival.time=y, censoring.status=ic,
#' geneid=as.character(1:nrow(x)), genenames=paste("g",as.character(1:nrow(x)),sep=
#' ""))
#' 
#' # train model
#' a3<- pamr.train(d, ngroup.survival=2)
#' 
#' 
#' pamr.plotstrata(a3, d$survival.time, d$censoring.status)
#' 
#' @export pamr.plotstrata
pamr.plotstrata <-
function (fit, survival.time, censoring.status)
{
    group <-apply(fit$proby,1,which.is.max)
#    require(survival)
    n.class <- length(unique(group))
    junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
    junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
 
  pv <- 1-pchisq(2*(junk2$loglik[2]-junk2$loglik[1]),df=n.class-1)

if(!is.null(fit$cutoffs.survival)){
    labels <- rep(NULL,n.class)
    labels[1] <- paste("(1)   ","<= ", round(fit$cutoffs.survival[1],2),sep="")
    if(n.class>2){
        for(i in 2:(n.class-1)){
          labels[i] <- paste("(",as.character(i),")  ", " > ",
        round(fit$cutoffs.survival[i-1],2), "  & <= ", 
        round(fit$cutoffs.survival[i],2), sep="")
     }}
    labels[n.class] <-  paste("(",as.character(n.class),")  ", " > ",round(fit$cutoffs.survival[n.class-1],2),sep="")
  }

else{labels <- as.character(1:n.class)}

#    win.metafile()
    plot(junk, col = 2:(2 + n.class - 1), xlab = "Time", ylab = "Probability of survival", main="Survival Strata Plot")
 #   legend(0.7 * max(fit$survival.time), 0.9, col = 2:(2 + n.class -
     legend(.01* max(fit$survival.time), 0.2, col = 2:(2 + n.class -
        1), lty = rep(1, n.class), legend = labels)
     text(0.1 * max(fit$survival.time), .25, paste("pvalue=",as.character(round(pv,4))))

#   dev.off()
#   return(TRUE)
  }
