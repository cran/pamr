#this is the wrong place but we need to put the methods stuff
#somewhere
#.initClasses <- function(env) {
#    setClass("uarray", representation(uexpr="character",
#                                      samplenames="character",
#                                      genenames="character",
#                                      samplecluster="dendrogram",
#                                      genecluster="dendrogram",
#                                      phenodata="character") )
#
#     #define a generic for obtaining the data
#    setGeneric("uexpr", function(object) standardGeneric("uexpr"))
#
#    setMethod("uexpr", "uarray", function(object) get(object@uexpr))
#
#    #define a generic for obtaining the phenotypic data
#    setGeneric("phenodata", function(object)
#               standardGeneric("phenodata"))
#
#    setMethod("phenodata", "uarray", function(object) {
#        if( length(object@phenodata) == 1 )
#            get(object@phenodata)
#        else
#            NULL
#    })
#
# #deal with the names
#    setGeneric("samplenames", function(object)
#               standardGeneric("samplenames"))
#    setMethod("samplenames", "uarray", function(object)
#              object@samplenames)
#
#    setGeneric("genenames", function(object)
#               standardGeneric("genenames"))
#    setMethod("genenames", "uarray", function(object) object@genenames )
#
#                                        # deal with the clusters
#    setGeneric("samplecluster", function(object)
#               standardGeneric("samplecluster"))
#  setMethod("samplecluster", "uarray", function(object) object@samplecluster )
#
#   setGeneric("genecluster", function(object)
#              standardGeneric("genecluster"))
#   setMethod("genecluster", "uarray", function(object) object@genecluster )
#
#                                       # plotting
#   if( !isGeneric("plot") )
#       setGeneric("plot")
#
#   setMethod("plot", "uarray", function(x, ...) {
#    expr <- as.matrix(uexpr(x))
#    #scale
#    expr <- sweep(expr, 1, apply(expr, 1, mean, na.rm = TRUE))
#    f <- function(v) {
#        v <- v[!is.na(v)]
#        sqrt(sum(v^2)/max(1, length(v) - 1))
#    }
#    expr <- sweep(expr, 1, apply(expr, 1, f), "/")
#    breaks <- seq(-3,3,by=.2)
#    colors<- GetColor(breaks)
#    breaks <- c(-100,breaks,100)
#    colors <- c(colors[1], colors)
#    opar<-par(mar=c(1,1,4,10))
#    on.exit(par(mar=opar))
#    image(1:ncol(expr), 1:nrow(expr), z = t(expr), axes = F,
#          col=colors, breaks=breaks, xlab="", ylab="")
#    axis(3, at=1:ncol(expr), labels=samplenames(x),tick=FALSE)
#    axis(4, at=1:nrow(expr), labels=genenames(x), tick=FALSE, las=1)
#})
#
#}
#given a set of labels, a set of weights and a set of clusters
#order the labels, by weights within clusters

#order.restricted <- function(labels, weights, clusters) {
#    n <- length(labels)
#    if (length(weights) != n || length(clusters) != n)
#        stop("all arguments must be the same length")
#    which <- split(labels, clusters)
#    wts <- split(weights, clusters)
#    cwts <- sapply(wts, mean)
#    cord <- order(cwts)
#    rval <- NULL
#    for(j in cord )
#        rval <- c(rval, which[[j]][order(wts[[j]])])
#    rval
#}

