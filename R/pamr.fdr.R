#' A function to estimate false discovery rates for the nearest shrunken
#' centroid classifier
#' 
#' A function to estimate false discovery rates for the nearest shrunken
#' centroid classifier
#' 
#' \code{pamr.fdr} estimates false discovery rates for a nearest shrunken
#' centroid classifier
#' 
#' @param trained.obj The result of a call to pamr.train
#' @param data Data object; same as the one passed to pamr.train
#' @param nperms Number of permutations for estimation of FDRs.  Default is 100
#' @param xl.mode Used by Excel interface
#' @param xl.time Used by Excel interface
#' @param xl.prevfit Used by Excel interface
#' @return A list with components: \item{results}{Matrix of estimates FDRs for
#' various various threshold values. Reported are both the median and 90th
#' percentile of the FDR over permutations} \item{pi0}{The estimated proportion
#' of genes that are null, i.e. not significantly different}
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
#'                genenames=paste("g",as.character(1:nrow(x)),sep=""))
#' 
#' mytrain <-   pamr.train(mydata)
#' myfdr <- pamr.fdr(mytrain, mydata)
#' 
#' @export pamr.fdr
pamr.fdr <- function(trained.obj, data, nperms=100, xl.mode=c("regular","firsttime","onetime","lasttime"),
                        xl.time=NULL, xl.prevfit=NULL){

this.call <- match.call()
 xl.mode=match.arg(xl.mode)


 if(xl.mode=="regular" | xl.mode=="firsttime"){

  y= data$y
  m=nrow(data$x)

 nclass=length(table(y))

  threshold <- trained.obj$threshold
 n.threshold=length(threshold)

  tt <- scale((trained.obj$centroids - trained.obj$centroid.overall)/trained.obj$sd, FALSE,
        trained.obj$threshold.scale * trained.obj$se.scale)


  ttstar <- array(NA,c(m,nperms,nclass))
results=NULL
pi0=NULL

}
  if(xl.mode=="onetime" |  xl.mode=="lasttime"){
 y=xl.prevfit$y
 m=xl.prevfit$m
 nclass=xl.prevfit$nclass
 threshold=xl.prevfit$threshold
 n.threshold=xl.prevfit$n.threshold
tt=xl.prevfit$tt
ttstar=xl.prevfit$ttstar
nperms=xl.prevfit$nperms
results=xl.prevfit$results
pi0=xl.prevfit$pi0
  }


  if(xl.mode=="regular"){
    first=1;last=nperms
  }
  if(xl.mode=="firsttime"){
    first=1;last=1
  }
  if(xl.mode=="onetime"){
    first=xl.time;last=xl.time
  }
  if(xl.mode=="lasttime"){
    first=nperms;last=nperms
  }


  for(i in first:last){
    cat("",fill=T)
     cat(c("perm=",i),fill=T)
    ystar <- sample(y)
    data2 <- data
    data2$y <- ystar
    foo<-pamr.train(data2, threshold=0, scale.sd = trained.obj$scale.sd,
    threshold.scale =  trained.obj$threshold.scale,
    se.scale = trained.obj$se.scale, offset.percent = 50, hetero = trained.obj$hetero,
prior = trained.obj$prior,  sign.contrast = trained.obj$sign.contrast)


   sdstar=foo$sd-foo$offset+trained.obj$offset
    ttstar[,i,] =scale((foo$centroids - foo$centroid.overall)/sdstar, FALSE,
        foo$threshold.scale * foo$se.scale)
}

 if(xl.mode=="regular" | xl.mode=="lasttime"){

fdr=rep(NA,n.threshold)
fdr90=rep(NA,n.threshold)
ngenes=rep(NA,n.threshold)

for(j in 1:n.threshold){
 # nobs=sum( (abs(tt)-threshold[j])%*%rep(1,ncol(tt)) >0)
 nobs=sum(drop((abs(tt)-threshold[j] > 0) %*% rep(1, ncol(tt))) > 0)

 temp=abs(ttstar)-threshold[j] >0
 temp2=rowSums(temp, dims=2)
 nnull=colSums(temp2>0)
  fdr[j]=median(nnull)/nobs
  fdr90[j]=quantile(nnull,.9)/nobs
  ngenes[j]=nobs
}


  q1 <- quantile(ttstar, .25)
  q2 <- quantile(ttstar, .75)

  pi0 <- min(sum( tt> q1 & tt< q2 )/(.5*m*nclass) ,1 )

  fdr <- fdr*pi0
fdr90=fdr90*pi0
fdr=pmin(fdr,1)
fdr90=pmin(fdr90,1)

  results <- cbind(threshold, ngenes, fdr*ngenes, fdr, fdr90)
om=is.na(fdr)
results=results[!om,]


 dimnames(results) <- list(NULL,c("Threshold", "Number of significant genes", "Median number of null genes",
"Median FDR", "90th percentile of FDR"))
# last time through, delete the temp stuff that is used just by Excel interface
y=NULL;x=NULL;m=NULL;threshold=NULL;n.threshold=NULL;tt=NULL;nperms=NULL;
ttstar=NULL
}

  return(list(results=results,pi0=pi0, y=y,m=m,threshold=threshold,n.threshold=n.threshold, tt=tt,ttstar=ttstar, nperms=nperms))
}
