pamr.cv <-
function(fit, data, nfold = min(table(data$y)), folds = balanced.folds(data$y),...)
{
        x <- data$x[fit$gene.subset, fit$sample.subset]
        if(is.null(fit$newy)) {
                y <- factor(data$y[fit$sample.subset])
        }
        else {
                y <- factor(data$newy[fit$sample.subset])
        }
        this.call <- match.call()
        junk <- nsccv(x, y, object = fit, ...)
        junk$call <- this.call
        junk$newy <- fit$newy
        junk$sample.subset <- fit$sample.subset
        return(junk)
}

