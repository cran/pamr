#' A function that interactively leads the user through a PAM analysis
#' 
#' A function that interactively leads the user through a PAM analysis
#' 
#' \code{pamr.menu} provides a menu for training, cross-validating and plotting
#' a nearest shrunken centroid analysis.
#' 
#' @param data A list with at least two components: x- an expression genes in
#' the rows, samples in the columns), and y- a vector of the class labels for
#' each sample. Same form as data object used by pamr.train.
#' @author Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan, and
#' Gilbert Chu
#' @examples
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(120)
#' x <- matrix(rnorm(1000*20),ncol=20)
#' y <- sample(c(1:4),size=20,replace=TRUE)
#' mydata <- list(x=x,y=y)
#' #  pamr.menu(mydata)
#' 
#' @export pamr.menu
pamr.menu <- function(data) {
  done <- FALSE
  junk.train <- NULL
  junk.results <- NULL
  while(!done) {
    cat("", fill = TRUE)
    switch(menu(c("pamr.train", "pamr.cv", "pamr.plotcv", 
                  "pamr.plotcen", "pamr.confusion", 
                  "pamr.plotcvprob", "pamr.geneplot", 
                  "pamr.listgenes", 
                  "pamr.train with heterogeneity analysis", 
                  "Exit")),
           junk.train <- pamr.train(data),
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               junk.results <- pamr.cv(junk.train, data)
             }
           }
           ,
           {
             if(is.null(junk.results)) {
               cat("Error: need to run pamr.cv first", fill
                   = TRUE)
             }
             if(!is.null(junk.results)) {
               pamr.plotcv(junk.results)
             }
           }
           ,
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.plotcen(junk.train, data, threshold = 
                            threshold)
             }
           }
           ,
           {
             if(is.null(junk.results)) {
               cat("Error: need to run pamr.cv first", fill
                   = TRUE)
             }
             if(!is.null(junk.results)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.confusion(junk.results, threshold = 
                              threshold)
             }
           }
           ,
           {
             if(is.null(junk.results)) {
               cat("Error: need to run pamr.cv first", fill
                   = TRUE)
             }
             if(!is.null(junk.results)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.plotcvprob(junk.results, data, threshold
                               = threshold)
             }
           }
           ,
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.geneplot(junk.train, data, threshold = 
                             threshold)
             }
           }
           ,
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.listgenes(junk.train, data, threshold = 
                              threshold)
             }
           }
           ,
           {
             junkk.train <- NULL
             cat("Normal class?", fill = TRUE)
             normal <- scan("", nlines = 1, what = "")
             junk.train <- pamr.train(data, hetero = normal)
           }
           ,
           done <- TRUE)
  }
  cat("Done\n")
}

pamr.pairscore <-function(x, pair.ind=NULL) {
}

