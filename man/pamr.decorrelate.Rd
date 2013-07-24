\name{pamr.decorrelate}
\alias{pamr.decorrelate}
\title{ A function to decorrelate (adjust) the feature matrix with respect
to some additional predictors}
\description{
A function to decorrelate (adjust) the feature matrix with respect
to some additional predictors
}
\usage{
 pamr.decorrelate(x, adjusting.predictors, xtest=NULL, adjusting.predictors.test=NULL)
}


\arguments{
  \item{x}{Matrix of training set feature values,  with
    genes in the rows, samples in the columns}
  \item{adjusting.predictors}{List of training set predictors to be used for adjustment}
   \item{xtest}{Optional matrix of test set feature values, to be adjusted
in the same way as the training set}
 \item{adjusting.predictors.test}{Optional list of test set predictors to be used for adjustment}
}

\details{
  \code{pamr.decorrelate}
 Does a least squares regression of each row of x on the adjusting
predictors, and returns the residuals. If xtest is provided, it also
returns the adjusted version of  xtest, using the
training set  least squares regression model for adjustment
}

\value{
  A list with components
  \item{x.adj}{Adjusted x matrix} 
  \item{xtest.adj}{Adjusted xtest matrix, if xtest we provided} 
}

\references{
  Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu 
  Diagnosis of multiple cancer types by shrunken centroids of gene expression 
  PNAS 99: 6567-6572.   Available at www.pnas.org}




\author{ Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
#generate some data
set.seed(120)

x<-matrix(rnorm(1000*20),ncol=20)
y<-c(rep(1,10),rep(2,10))
adjusting.predictors=list(pred1=rnorm(20), pred2=as.factor(sample(c(1,2),replace
=TRUE,size=20)))
xtest=matrix(rnorm(1000*10),ncol=10)
adjusting.predictors.test=list(pred1=rnorm(10), pred2=as.factor(sample(c(1,2),replace
=TRUE,size=10)))

# decorrelate training x wrt adjusting predictors

x.adj=pamr.decorrelate(x,adjusting.predictors)$x.adj
# train classifier with adjusted x

d=list(x=x.adj,y=y)
a<-pamr.train(d)

# decorrelate training and test x wrt adjusting predictors, then make
#predictions for test set

temp <- pamr.decorrelate(x,adjusting.predictors, xtest=xtest,
                         adjusting.predictors.test=adjusting.predictors.test)

d=list(x=temp$x.adj,y=y)
a<-pamr.train(d)
aa<-pamr.predict(a,temp$xtest.adj, threshold=.5)
}

