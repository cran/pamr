\name{pamr.surv.to.class2}
\alias{pamr.surv.to.class2}
\title{A function to  assign observations to categories, based on their
survival times.}
\description{A function to  assign observations to categories, based on their
survival times.}


\usage{
pamr.surv.to.class2(y, icens, cutoffs=NULL, n.class=NULL,  class.names=NULL, newy=y, newic=icens)
}

\arguments{
  \item{y}{vector of survival times}
\item{icens}{Vector of censorng status values: 1=died, 0=censored}
\item{cutoffs}{Survival time cutoffs for categories. Default NULL}
\item{n.class}{Number of classes to create: if cutoffs is NULL,
     n.class equal classes are created.}
\item{class.names}{Character names for classes}
\item{newy}{New set of survival times, for which probabilities are computed
(see below). Default is y}
\item{newic}{New set of censoring statuses, for which probabilities are computed
(see below). Default is icens}

}
  
\details{
\code{pamr.pamr.surv.to.class2}  splits observations into categories
based on their survival times and the Kaplan-Meier estimates.
For example if n.class=2, it makes
two categories, one below the median survival, the other above.
For each observation (newy, ic), it then computes the probability
of that observation falling in each category. For an uncensored observation
that probability is just  1 or 0 depending on when the death occurred.
For a censored observation, the probabilities are based on the Kaplan
Meier and are typically between 0 and 1.
}

\value{
  \item{class} {The category labels}
 \item{prob}{The estimates class probabilities}
 \item{cutoffs}{The cutoffs used}
}

\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu}

\examples{

gendata<-function(n=100, p=2000){
  tim <- 3*abs(rnorm(n))
  u<-runif(n,min(tim),max(tim))
  y<-pmin(tim,u)
   ic<-1*(tim<u)
m <- median(tim)
x<-matrix(rnorm(p*n),ncol=n)
  x[1:100, tim>m] <-  x[1:100, tim>m]+3
  return(list(x=x,y=y,ic=ic))
}

# generate training data; 2000 genes, 100 samples

junk<-gendata(n=100)
y<-junk$y
ic<-junk$ic
x<-junk$x
d <- list(x=x,survival.time=y, censoring.status=ic,
geneid=as.character(1:nrow(x)), genenames=paste("g",as.character(1:nrow(x)),sep=
""))

# train model
a3<- pamr.train(d, ngroup.survival=2)

# generate test data
junkk<- gendata(n=500)

dd <- list(x=junkk$x, survival.time=junkk$y, censoring.status=junkk$ic)

# compute soft labels
proby <-  pamr.surv.to.class2(dd$survival.time, dd$censoring.status,
             n.class=a3$ngroup.survival)$prob

}
\keyword{}

