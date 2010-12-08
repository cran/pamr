\name{pamr.makeclasses}
\alias{pamr.makeclasses}
\title{A function to interactively define classes from a clustering tree}
\description{function to interactively define classes from a clustering tree}
\usage{
pamr.makeclasses(data,sort.by.class=FALSE,...)
}

\arguments{
  \item{data}{The input data. A list with components: x- an expression
    genes in the rows, samples in the columns, and y-  a vector of
    the class labels for each sample, and batchlabels- a vector of batch
    labels for each sample.
    This object if the same form as that produced by pamr.from.excel.}
  \item{sort.by.class}{Optional argument. If true, the clustering tree
    is forced to put all samples in the same class  (as defined by the
    class labels y in `data') together in the tree. This is useful
    if a regrouping of classes is desired. Eg: given classes 1,2,3,4
    you want to define new classes (1,3) vs (2,4)  or 2 vs (1,3)}
  \item{...}{Any additional arguments to be passed to hclust}
}

\details{
  \code{pamr.makeclasses} 
  Using this function the user interactively defines a new set of classes,
  to be used in pamr.train, pamr.cv etc. 
  After invoking pamr.makeclasses, a clustering tree is drawn.
  This callss the R function \code{hclust}, and any arguments for
  \code{hclust} can be passed to it.
  Using the left button, the user clicks at the junction point defining
  the subgroup 1. More groups can be added to class 1 by clicking
  on further junction points. The user ends the definition of class 1 by
  clicking on the rightmost button [in Windows,  an additional menu appears
  and he chooses Stop] . This process is continued for classes 2,3 etc.
  Note that some sample may be left out of the new classes.
  Two consecutive clicks of the right button ends the definition for all classes.
  
  At the end, the clustering is redrawn, with the new class labels shown.
  
  Note: this function is "fragile". The user must click close to the junction
  point, to avoid confusion with other junction points. Classes 1,2,3..
  cannot have samples in common (if they do, an Error message will appear).
  If the function is confused about the desired choices, it will
  complain and ask the user to rerun pamr.makeclasses. The user should
  also check that the labels on the  final redrawn cluster tree agrees with the desired 
  classes.
}

\value{A vector of class labels 1,2,3...
  If a component is NA (missing), then the
  sample is not assigned to any class.
  This vector should be assigned to the newy component of data, for use in pamr.train etc.
  Note that pamr.train uses the class labels in the component ``newy'' if
  it is present. Otherwise it uses the data labels ``y''.
}




\author{ Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and Gilbert Chu  }

\examples{
set.seed(120)
#generate some data
x <- matrix(rnorm(1000*40),ncol=40)
y <- sample(c(1:4),size=40,replace=TRUE)
batchlabels <- sample(c(1:5),size=40,replace=TRUE)

mydata <- list(x=x,y=factor(y),batchlabels=factor(batchlabels), geneid=as.character(1:nrow(x)), genenames=paste("g",as.character(1:nrow(x)),sep=""))

# mydata$newy <- pamr.makeclasses(mydata)    Run this and define some new classes

train <- pamr.train(mydata)
}


