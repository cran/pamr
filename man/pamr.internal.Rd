\name{pamr-internal}
\title{Internal pamr functions}
\alias{pamr.cube.root}
\alias{pamr.knnimpute.old}
\alias{pamr.pairscore}
\alias{pamr.pvalue.survival}
\alias{pamr.score.to.class1}
\alias{pamr.score.to.class2}
\alias{print.nsc}
\alias{print.nsccv}
\alias{print.pamrcved}
\alias{print.pamrtrained}
\alias{pamr.xl.compute.offset}
\alias{pamr.xl.error.trace}
\alias{pamr.xl.get.offset}
\alias{pamr.xl.derive.adjusted.prior}
\alias{pamr.xl.get.default.training.parameters}
\alias{pamr.xl.get.uniform.prior}
\alias{pamr.xl.get.sample.prior}
\alias{pamr.xl.get.class.names}
\alias{pamr.xl.get.class.labels}
\alias{pamr.xl.get.number.of.classes}
\alias{pamr.xl.process.data}
\alias{pamr.xl.compute.cv.confusion}
\alias{pamr.xl.compute.confusion}
\alias{pamr.xl.is.a.subset}
\alias{pamr.xl.listgenes.compute}
\alias{pamr.xl.plot.test.probs.compute}
\alias{pamr.xl.plot.training.error.compute}
\alias{pamr.xl.plotcen.compute}
\alias{pamr.xl.plotcv.compute}
\alias{pamr.xl.plotcvprob.compute}
\alias{pamr.xl.predict.test.class}
\alias{pamr.xl.predict.test.class.only}
\alias{pamr.xl.predict.test.probs}
\alias{pamr.xl.test.data.impute}
\alias{pamr.xl.test.errors.surv.compute}
\alias{pamr.xl.test.errors.compute}
\alias{pamr.xl.transform.class.labels}
\alias{pamr.xl.transform.data}
\alias{pamr.xl.transform.test.data}
\alias{pamr.xl.plotsurvival}
\alias{pamr.xl.plotsurvival.test}
\alias{pamr.xl.predict.test.surv.class}
\alias{pamr.xl.plotsurvival.strata}
\alias{pamr.xl.test.get.soft.classes}
\alias{pamr.xl.get.threshold.range}
\alias{pamr.xl.get.soft.class.labels}

\description{Internal pamr functions}
\usage{
pamr.pairscore(x, pair.ind=NULL)
pamr.pvalue.survival(group, survival.time, censoring.status, ngroup.survival)
pamr.score.to.class1(x, scores, cutoff=2, n.class=2)
pamr.score.to.class2(x, scores, cutoff=2, n.pc=1, n.class=2)
pamr.knnimpute.old(data, k = 10) 
pamr.cube.root(x)
print.nsc(x, ...)
print.nsccv(x, ...)
\method{print}{pamrtrained}(x, ...)
\method{print}{pamrcved}(x, ...)
pamr.xl.error.trace() 
pamr.xl.get.threshold.range(fit) 
pamr.xl.get.soft.class.labels(fit, survival.times, censoring.status) 
pamr.xl.compute.offset(data, offset.percent=50, prior=prior)
pamr.xl.get.offset() 
pamr.xl.derive.adjusted.prior(prior, data) 
pamr.xl.get.default.training.parameters(data) 
pamr.xl.get.uniform.prior(data, nclasses=NULL) 
pamr.xl.get.sample.prior(data) 
pamr.xl.get.class.names() 
pamr.xl.get.class.labels() 
pamr.xl.get.number.of.classes() 
pamr.xl.process.data(use.old.version=FALSE) 
pamr.xl.compute.cv.confusion (fit, cv.results, threshold) 
pamr.xl.compute.confusion (fit, threshold) 
pamr.xl.is.a.subset(a, y) 
pamr.xl.listgenes.compute (fit, data, threshold, fitcv=NULL,  genenames = FALSE) 
pamr.xl.plot.test.probs.compute(fit, new.x, newx.classes, missing.class.label,
                                threshold, sample.labels=NULL) 
pamr.xl.plot.training.error.compute(trained.object) 
pamr.xl.plotcen.compute(fit, data, threshold) 
pamr.xl.plotcv.compute(aa) 
pamr.xl.plotcvprob.compute(fit, data, threshold) 
pamr.xl.predict.test.class(fit, newx, threshold, test.class.labels) 
pamr.xl.predict.test.surv.class(fit, newx, threshold, survival.times,
                                censoring.status) 
pamr.xl.predict.test.class.only(fit, newx, threshold) 
pamr.xl.predict.test.probs(fit, newx, threshold) 
pamr.xl.test.data.impute(x, k, use.old.version=FALSE) 
pamr.xl.test.errors.surv.compute(fit, newx, threshold=fit$threshold,
                                 survival.times, censoring.status) 
pamr.xl.test.errors.compute(fit, newx, newx.classes, threshold=fit$threshold,
                            prior = fit$prior,
                            threshold.scale = fit$threshold.scale, ...) 
pamr.xl.transform.class.labels(x) 
pamr.xl.transform.data(data) 
pamr.xl.transform.test.data(test.x) 
pamr.xl.plotsurvival(fit, data, threshold) 
pamr.xl.plotsurvival.test(fit, newx, survival.time, censoring.status, threshold) 
pamr.xl.plotsurvival.strata(fit, data) 
pamr.xl.test.get.soft.classes(fit, survival.times, censoring.status)
pamr.xl.get.soft.class.labels(fit, survival.times, censoring.status) 

}
\author{Balasubramanian Narasimhan
 and Rob Tibshirani}
\details{
  These are not to be called by the user.
}
\keyword{internal}

