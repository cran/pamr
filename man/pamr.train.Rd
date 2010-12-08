\name{pamr.train}
\alias{pamr.train}
\title{ A function to train a nearest shrunken centroid
  classifier}
\description{
  A function that computes a nearest shrunken centroid for gene
  expression
  (microarray) data
}
\usage{
pamr.train(data, gene.subset=NULL, sample.subset=NULL,
         threshold = NULL, n.threshold = 30, 
        scale.sd = TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50,
         hetero=NULL, prior = NULL, remove.zeros = TRUE, sign.contrast="both",
ngroup.survival = 2)
}

\arguments{
  \item{data}{The input data. A list with components: x- an expression
    genes in the rows, samples in the columns), and y-  a vector of
    the class labels for each sample. 
    Optional components-
    genenames,
    a vector of gene names, and geneid- a vector of gene identifiers.}
  \item{gene.subset}{Subset of genes to be used.  Can be either
    a logical vector of length total number of genes, or a list
    of integers of the row numbers of the genes to be used}
  \item{sample.subset}{Subset of samples to be used.  Can be either
    a logical vector of length total number of samples, or a list
    of integers of the column numbers of the samples to be used.}
  \item{threshold}{A vector of threshold values for the centroid
    shrinkage.Default is a set of 30 values chosen by the software}
  \item{n.threshold}{Number of threshold values desired (default 30)}
  \item{scale.sd}{Scale each threshold by the wthin class standard
    deviations? Default: true}
  \item{threshold.scale}{Additional scaling factors to be applied
    to the thresholds. Vector of length equal to the number of
    classes.
    Default- a vectors of ones.}
  \item{se.scale}{Vector of scaling factors for the within class
    standard errors. Default is sqrt(1/n.class-1/n), where n is
    the overall sample size and n.class is the sample sizes in each
    class. This default adjusts for different class sizes.}
  \item{offset.percent}{Fudge factor added to  the denominator of each t-statistic,
    expressed as a percentile of the gene standard deviation values.
    This is a small positive quantity to penalize genes with
    expression
    values near zero, which can result in very large ratios. This
    factor
    is expecially impotant for Affy data. Default
    is the median of the standard deviations of each gene.}
 \item{hetero}{Should a heterogeneity transformation be done?
    If yes, hetero must be set to one of the class labels (see Details below).
    Default is no (hetero=NULL) }
  \item{prior}{Vector of length the number of classes, representing
    prior probabilities for each of the classes. The prior is used
    in Bayes rule for making class prediction. Default is NULL, and prior
is then taken to be  n.class/n,
    where  n is the overall sample size and n.class is the sample
    sizes in   each class.}
\item{remove.zeros}{Remove threshold values yielding zero genes? Default TRUE}
  \item{sign.contrast}{Directions of allowed deviations of class-wise average gene  expression from the 
    overall average  gene expression. Default is ``both'' (positive or negative).
    Can also  be set to ``positive'' or ``negative''.}
\item{ngroup.survival}{Number of groups formed for survival data. Default 2}
}
\details{
  \code{pamr.train} fits a nearest shrunken centroid classifier to gene
  expression data. Details may be found in the PNAS paper referenced
  below. One feature not described there is "heterogeneity analysis".
  Suppose there are two classes labelled "A" and "B".
  CLass "A" is considered a normal class, and "B" an abnormal class.
  Setting hetero="A" transforms   expression values x[i,j] to
  |x[i,j]- mean(x[i,j])| where the mean is taken only over samples in
  class "A". The transformed feature values are then used in Pam.
  This is useful when the abnormal class "B" is heterogeneous, i.e.
  a given gene might have higher expresion than normal for some
  class "B" samples, and lower for others.
  With more than 2 classes, each class is centered on the class specified
  by hetero.
}

\value{
  A list with components
  \item{y}{The outcome classes.} 
  \item{yhat}{A matrix of predicted classes, each column representing
    the results from one threshold.}.
  \item{prob}{A array of predicted class probabilities. of dimension
    n by nclass by n.threshold. n is the number samples, nclass is the number
    of classes, n.threshold is the number of thresholds tried}
  \item{centroids}{A matrix of (unshrunken) class centroids, n by
    nclass}
  \item{hetero}{Value of hetero used in call to pamr.train}
  \item{norm.cent}{Centroid of "normal" group, if hetero was specified}
  \item{centroid.overall}{A vector containing the (unshrunken) overall
    centroid (all classes together)}
  \item{sd}{A vector of the standard deviations for each gene}
  \item{threshold}{A vector of the threshold tried in the shrinkage}
  \item{nonzero}{A vector of the number of genes that survived the
    thresholding, for each threshold value tried}
  \item{threshold.scale}{A vector of threshold scale factors that were
    used}
  \item{se.scale}{A vector of standard error scale factors that were
    used}
  \item{call}{The calling sequence used}
  \item{prior}{The prior probabilities used}
  \item{errors}{The number of trainin errors for each threshold value}
}

\references{
  Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu 
  Diagnosis of multiple cancer types by shrunken centroids of gene expression 
  PNAS 99: 6567-6572.   Available at www.pnas.org}




\author{ Trevor Hastie,Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
#generate some data
set.seed(120)
x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)
mydata <- list(x=x,y=factor(y))

#train classifier
results<-   pamr.train(mydata)

# train classifier on all  data except class 4
results2 <- pamr.train(mydata,sample.subset=(mydata$y!=4))
 
# train classifier on  only the first 500 genes
results3 <- pamr.train(mydata,gene.subset=1:500)

}

