pamr.train <-

function(data, gene.subset=1:nrow(data$x), sample.subset=1:ncol(data$x),
         threshold = NULL, n.threshold = 30,
        scale.sd = TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50, hetero=NULL,
         prior = NULL,  remove.zeros = TRUE, sign.contrast="both")

{
        this.call <- match.call()

        if(is.null(prior))
          {prior <- table(data$y[sample.subset])/length(data$y[sample.subset])
           prior <- prior[prior!=0]
        }

        if(!is.null(sample.subset) & !is.null(data$newy)) {
           stop("Can't have both newy present in data object, and sample.subset specified"
                        )
                     }
       if(is.null(sample.subset)){sample.subset <-1:ncol(data$x)}

        if(is.null(data$y) & is.null(data$newy)) {
                stop("must have either y or newy present in data object")
        }
        if(is.null(data$newy)) {
                y <- data$y
        }
        if(!is.null(data$newy)) {
                y <- data$newy
                sample.subset <- (1:ncol(data$x))[!is.na(y)]
                print("Using classes `newy' from data object")
        }
        junk <- nsc(data$x[gene.subset, sample.subset], factor(y[sample.subset]), 
          offset.percent=offset.percent,  threshold = threshold, hetero=hetero,
          n.threshold = n.threshold,  scale.sd= scale.sd, threshold.scale=threshold.scale,
           se.scale= se.scale, prior=prior, remove.zeros=remove.zeros,
            sign.contrast=sign.contrast)

        junk$call <- this.call
        junk$gene.subset <- gene.subset
        junk$sample.subset <- sample.subset
        junk$newy <- data$newy
        return(junk)
}
