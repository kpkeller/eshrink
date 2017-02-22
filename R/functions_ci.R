##' @name check_CIbound
##' @title Test function for creating confidence intervals
##
##' @description 'Inverts the test' to determine if a given value should lie in a confidence region
##
##' @details This function is used as part of an 'inverting the test' approach to generate confidence intervals for estimators from \code{\link{festRidge}}. \code{Bstar} datasets are generated from slices of the posterior distribution of the model parameters where beta (or other parameter indicated by \code{ind}) is fixed at the value \code{testBeta}. For each dataset, beta is estimated via \code{\link{festRidge}}, and the resulting distribution of estimators is compared against the estimate from the observed data (\code{obsEst}).
##' 
##' The values of \code{lambdaseq}, \code{X}, \code{nPost}, and \code{loss} are passed to \code{\link{festRidge}} and typically match the values that were used to compute \code{obsEst}.
##'
##' The computational cost of this function is most affected by the values of
##'	\code{nPost} and \code{Bstar}. Large values of the latter are important for 
##'	adequately representing the distribution of parameter estimates. In some
##'	settings, \code{nPost} can be reduced without substantially impacting
##'	the results. However, each dataset is likely to be different.
##
##' @param testBeta Candidate value of beta to test.
##' @param obsEst Estimate of beta from the observed data for which a confidence interval is desired
##' @param posthyper List of hyperparameters for the posterior distribution of model parameters. See \code{\link{samplePosterior}} for expected names.
##' @param type String indicating ridge or LASSO
##' @param lambdaseq Sequence of penalty values
##' @param X deisgn matrix
##' @param nPost Number of posterior samples to use.
##' @param ind Index of parameter to test. Defaults to 1.
##' @param loss Either \code{"fMBV"} or \code{"fMSE"}.
##' @param Bstar Number of estimators to compute for comparison distribution. Larger values improve the precision of the procedure but increase computational cost.
##' @param B Passed to \code{\link{festLASSO}}
##' @param lowerBound Logical indicating that the test is for a lower bound
##' @param reproducible Should the simulated datasets be reproducible?
##' @param alpha Percentile of the distribution to compare against. See details.
##' @param returnDist If TRUE, then distribution of estimates generated is returned
##'			instead of comparison against \code{alpha}
##'
##' @export
##' @importFrom stats uniroot
##' @author Joshua Keller
check_CIbound <- function(testBeta, obsEst, type=c("ridge", "lasso"), posthyper,  lambdaseq, X, nPost, ind=1, Bstar=100,  B=500, loss="fMBV", lowerBound=TRUE, reproducible=TRUE, alpha=0.025, returnDist=FALSE, ...) {
	type <- match.arg(type)
	p <- ncol(X)
	XtX <- crossprod(X)
	gammaMu <- posthyper $postMu[-ind] + posthyper $postV[-ind, ind]/posthyper $postV[ind, ind]*(testBeta - posthyper $postMu[ind])
	gammaV <- posthyper $postV[-ind, -ind] - tcrossprod(posthyper $postV[-ind, ind], posthyper $postV[ind, -ind])/posthyper $postV[ind, ind]
	if(reproducible) set.seed(ind[1])
	sigma2star <- 1/rgamma(n= Bstar, shape= posthyper $an, rate= posthyper $bn)
	gammastar <- MASS::mvrnorm(n= Bstar, mu=rep(0, p-1), Sigma= gammaV)
	gammastar <- t(t(gammastar*sqrt(sigma2star)) + as.vector(gammaMu))
	betafMBVstar <- numeric(Bstar)
	XtXlamIinv <- lapply(lambdaseq, getXtXlamIinv, XtX=XtX)
	temp_paramvector <- numeric(p)
	temp_paramvector[ind] <- testBeta
	for (j in 1:Bstar){
		temp_paramvector[-ind] <- gammastar[j,]
		if(reproducible) set.seed(j*7)
		Xbetastar <- X %*% temp_paramvector
		eps <- rnorm(n=nrow(X), sd=sqrt(sigma2star[j]))
		ystar <- Xbetastar + eps
		ystar <- ystar - mean(ystar)
		if (type=="ridge"){
			bhatstar <- festRidge(X=X, y=ystar, loss=loss, lseq= lambdaseq, nPost=nPost, se.version="none", XtXlamIinv= XtXlamIinv, ind=ind, ...)[[loss]]
		} else if (type=="lasso"){
			bhatstar <- festLASSO(X=X, y=ystar, loss=loss, lseq= lambdaseq, nPost=nPost, se.version="none", ind=ind, B=B,...)[[loss]]			
		}
		betafMBVstar[j] <- bhatstar$beta[ind] 
	}
	if (returnDist){
		return(betafMBVstar)
	} else {
		if (lowerBound) {
			mean(obsEst <betafMBVstar) - alpha
		} else {
			mean(obsEst> betafMBVstar) - alpha
		}
	}
}

##' @name invertTest
##'	@rdname check_CIbound
##
##' @param interval Interval to check. Used for both upper and lower bound, if they are not provided
##' @param lower.interval,upper.interval Bounding intervals over which to check for lower and upper endpoints of CI
##' @param ... In \code{invertTest }, these are passed to \code{check_CIbound}. In 
##'		\code{check_CIbound}, these arguments are passed to \code{\link{samplePosterior}}.
##' @param tol Passed to \code{\link{uniroot}}
##' @param fulldetail If TRUE, then output from \code{\link{uniroot}} is included.
##
## @description Confidence intervals for festRidge and festLASSO estimators
## @details A wrapper for \code{\link{check_CIbound}} that uses a root-finding
## search to find boundaries.
##' @export
## @author Joshua Keller
invertTest <- function(interval, type="ridge", lower.interval=interval, upper.interval=interval, ..., tol=0.005, fulldetail=FALSE){
	
	ui_lower <- uniroot(f= check_CIbound, lowerBound=TRUE, type= type, ...,  interval= lower.interval, tol=tol)
lowerBound <- ui_lower $root
	ui_upper <- uniroot(f= check_CIbound, lowerBound=FALSE, type= type,   ...,  interval=upper.interval, tol=tol)
upperBound <- ui_upper $root
out <- c(lower=lowerBound, upper=upperBound)
if (fulldetail){
	suppressWarnings(out$lower_output <- ui_lower)
	out$upper_output <- ui_upper
}
out
}





