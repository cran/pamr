#' A function that takes estimate class probabilities and produces a class
#' prediction or indeterminate prediction
#' 
#' A function that takes estimate class probabilities and produces a class
#' prediction or indeterminate prediction
#' 
#' 
#' @param prob Estimated class probabilities, from pamr.predict with
#' type="posterior")
#' @param mingap Minimum difference between highest and second highest
#' probability. If difference is < mingap, prediction is set to indeterminate
#' (NA)
#' @author Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' mydata <- list(x=x,y=y)
#' mytrain <-   pamr.train(mydata)
#' prob<- pamr.predict(mytrain, mydata$x , threshold=1, type="posterior")
#' pamr.indeterminate(prob,mingap=.75)
#' 
#' @export pamr.indeterminate
pamr.indeterminate <-  function(prob, mingap=0){
n=nrow(prob)
yhat=rep(NA,n)
for(i in 1:n){

  r=rank(-prob[i,])
if(sum(r==1)==1){
 pr1=prob[i,r==1]
 pr2=prob[i,r==2]
if(pr1-pr2 >= mingap){ yhat[i]=(1:ncol(prob))[r==1]}
}}
yhat=as.factor(dimnames(prob)[[2]][yhat])
return(yhat)
}
