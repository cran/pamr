nsccv <-
function(x, y, nfold = min(table(y)), folds = balanced.folds(y), threshold =
        NULL, threshold.scale = NULL, prior, object, ...)
{
        this.call <- match.call()
        n <- length(y)
        if(is.null(folds)) {
                folds <- split(sample(1:n), rep(1:nfold, length = n))
        }
        else nfold <- length(folds)
        if(missing(prior)) {
                if(missing(object))
                        prior <- table(y)/n
                else prior <- object$prior
        }
        if(missing(threshold)) {
                if(missing(object))
                        stop("Must either supply threshold argument, or an nsc object"
                                )
                else {
                        threshold <- object$threshold
                        threshold.scale <- object$threshold.scale
                        se.scale <- object$se.scale
                }
        }
        n.threshold <- length(threshold)        ### Set up the data structures
        yhat <- rep(list(y), n.threshold)
        names(yhat) <- paste(seq(n.threshold))
        yhat <- data.frame(yhat)
        n.class <- table(y)
        prob <- array(1, c(n, length(n.class), n.threshold))
        size <- double(n.threshold)
        hetero <-object$hetero
        for(ii in 1:nfold) {
                cat("Fold", ii, ":")
                a <- nsc(x[,  - folds[[ii]]], y[ - folds[[ii]]], x[, folds[[ii
                        ]], drop = FALSE], threshold = threshold, threshold.scale
                         = threshold.scale, se.scale = se.scale, prior = prior,
                          hetero=hetero,
                        ..., remove.zeros = FALSE)
                size <- size + a$nonzero
                prob[folds[[ii]],  ,  ] <- a$prob
                yhat[folds[[ii]],  ] <- a$yhat
                cat("\n")
        }
        if(missing(object))
                size <- round(size/nfold)
        else size <- object$nonzero
        error <- rep(NA, n.threshold)
        loglik <- error
        for(i in 1:n.threshold) {
                error[i] <- sum(yhat[, i] != y)/n
                loglik[i] <- sum(log(prob[,  , i][cbind(seq(1, n), codes(y))]))/                        n
        }
obj<- list(threshold=threshold, error=error, loglik=loglik,size=size, yhat=yhat,y=y,prob=prob,folds=folds,
                call = this.call)
        class(obj) <- "nsccv"
        obj
}

