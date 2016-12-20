#####################################
#####################################
#####################################
# Functions for the Ridge Estimator
#
# Included in this file:
#
#	ridgeEst()
#	ridgeBias()
#	ridgeVar()
#	ridgeMSE()
#	vcovRidgefEst()
#	ridgefEst()
#


##' @name ridgeEst
##' @title Compute Ridge Estimate
##
##' @description Computes a ridge estimate for a given regression problem
##
##' @details Computes the ridge estimator from the regression of y on X, with penalty (vector!) penalize.
##'
##' The input \code{penalize} provides a vector of penalty factors, 
##'	such that the penalty for covariate j is lambda*penalize[j].
##' Primary purpose is for indicating which variable to penalize (1)
##' and which to not (0), although more flexible penalties that include
##’ decimal values are possible. Defaults to c(0, rep(1, p-1)), where
##’ p is number of columns in X (this penalizes all coefficients but 
##’ the first).  
##
##  Input:
##' @param lambda Penalty factor to be applied
##' @param X Design matrix of regression problem
##' @param y outcome vector. Typically centered.
##' @param penalize Vector giving penalty structure. Values must be in [0, 1].
##’	See Details for more information.
##' @param XtX the cross product of the design matrix. Will be computed from
##'		\code{X} if not provided.
##
##' @export
##' @seealso \code{\link{ridgefEst}}, \code{\link{ridgeMSE}}
ridgeEst <- function(lambda, X, y,  penalize, XtX=crossprod(X)){
	if(missing(penalize)) penalize <- c(0, rep(1, ncol(XtX)-1))
	penalize <- as.numeric(penalize)
	if (any(penalize>1) | any(penalize<0)) stop("Element of penalize must be between 0 and 1.")
	penalty <- lambda*penalize
	as.vector(solve(XtX + diag(penalty), crossprod(X, y)))
}


##' @name ridgeBias
##' @title Get Bias of Ridge Estimator
##
##' @description Computes the (analytic) bias of Ridge estimator
##
##' @details If \code{beta} is provided as a matrix, this will
##’	treat each column of \code{beta} as a different true
##’	parameter vector and return a matrix of bias vectors.
##
##  Input:
##' @param lambda Penalty parameter (should have length 1).
##' @param beta True parameter vector. Can also be p x n matrix.
##' @param XtX cross product of design matrix
##' @param penalize vector giving what to penalize. See \code{\link{ridgeEst}} for more information.
##' @param ind Numerical or logical vector indicating which elements
##'		should be returned. Defaults to the first element.
##' @param XtXlamIinv explicit value of (X'X + diag(penalty))^{-1}.  Useful
##'		for simulations to save computation.
##
##' @export
##' @author Joshua Keller
##'	@seealso \code{\link{ridgeMSE}, \link{ridgeVar}}
ridgeBias <- function (lambda, beta, XtX, penalize=c(0, rep(1, ncol(XtX)-1)), ind = 1, XtXlamIinv=NULL) {
    if (is.null(XtXlamIinv)){
    	p <- ncol(XtX)
	    if(nrow(XtX)!=p) stop("XtX must be square.")
	    penalize <- as.numeric(penalize)
	    if (max(abs(penalize)) > 1) 
	        stop("Element of penalize must be <= 1.")
	    if (length(lambda) > 1) {
	        lambda <- lambda[1]
	        warning("lambda should be of length 1. Only using first value.")
	    }
	    penalty <- lambda * penalize
	    XtXlamIinv <- chol2inv(chol(XtX + diag(penalty)))
    } else {
	    p <- ncol(XtXlamIinv)	
    }
    if (is.matrix(beta) && nrow(beta) != p) 
        stop("Beta should be vector, or matrix with ncol(XtX) rows.")
#    XtXlamIinvB <- as.matrix(solve(XtX + diag(penalty), penalty * beta))
#	cholXtXlamInv <- chol(XtX + diag(penalty))
#	XtXlamIinvB  <- backsolve(cholXtXlamInv, backsolve(cholXtXlamInv, penalty*beta, transpose=T))
#	XtXlamIinvB <- chol2inv(chol(XtX + diag(penalty))) %*%  (penalty * beta)
#    return(XtXlamIinvB[ind, ])
	XtXlamIinvB <- crossprod(penalty*t(XtXlamIinv[ind, , drop=F]),  beta)
    return(XtXlamIinvB)
}


##' @name ridgeVar
##' @title Get Ridge variance
##
##' @description To be added..
##
##' @details See \code{\link{ridgeEst}} for description of \code{penalize} parameter.
##
##  Input:
##' @param XtX cross product of design matrix
##' @param lambda penalty parameter
##'	@param penalize vector giving what to penalize
##' @param ind Numerical or logical vector indicating which elements
##'		of the variance matrix should be returned. Defaults to the
##'		(1,1) element
##' @param sigma2 The true variance parameter
##' @param XtXlamIinv explicit value of (X'X + diag(penalty))^{-1}.  Useful
##'		for simulations to save computation.
##' @export
##' @author Joshua Keller
##' @seealso \code{\link{ridgeMSE}, \link{ridgeBias}}
ridgeVar <- function(lambda, XtX, penalize, sigma2=1, ind=1, XtXlamIinv=NULL){
	if (is.null(XtXlamIinv)) {
		if(missing(penalize)) penalize <- c(0, rep(1, ncol(XtX)-1))
		penalize <- as.numeric(penalize)
		if (max(abs(penalize))>1) stop("Element of penalize must be <= 1.")
		penalty <- lambda*penalize	
		XtXlamIinv <- chol2inv(chol(XtX + diag(penalty)))
	}
	rVar <-  sigma2 * XtXlamIinv %*% XtX %*% XtXlamIinv
	return(rVar[ind, ind])
}







##' @name ridgeMSE
##' @title Compute MSE for Ridge Estimator
##
##' @description Computes the MSE for ridge estimators given different
##'		values of the true \code{beta} and \code{sigma2} coefficients.
##
##' @details The bias and variance are computed using the functions
##'		\code{\link{ridgeBias}} and \code{\link{ridgeVar}}. 
##'
##'		The \code{XtXlamIinvB} and \code{returmMSE} arguments were designed
##'		for simulations, where pre-computation of matrices and reducing
##'		output provides computational benefits.  These are typically
##'		not needed for analysis of a single dataset.
##
##  Input:
##' @param lambda Shrinkage penalty parameters (vector).
##'	@param XtX the cross product of the design matrix.
##' @param beta true value of the regression parameters
##' @param sigma2 true value of the variance parameter
##' @param ind numeric or logical vector indicating which
##' 		element for which MSE should be returned
##' @param XtXlamIinv explicit value of (X'X + diag(penalty))^{-1}.  Useful
##'		for simulations to save computation.
##'
##' @export
##
ridgeMSE <- function(lambda, XtX, beta, sigma2, ind=1, XtXlamIinv=NULL){
	if (is.logical(ind)){
		if (sum(ind)>1) stop("'ind' should be a single integer or logical vector with sum 1.")
	} else {
		if (length(ind)>1)stop("'ind' should be a single integer or logical vector with sum 1.")
	}
	
	if (is.matrix(beta) && nrow(beta)!=ncol(XtX)) stop("Beta should be vector, or matrix with ncol(XtX) rows.")

	
	if (length(lambda)==1){
		rv <- ridgeVar(lambda= lambda, XtX=XtX, ind=ind)
		rv <- as.vector(tcrossprod(sigma2, rv))
		rb <- as.vector(ridgeBias(lambda, XtX, beta=beta, ind=ind))
		
#		cat("dim rb = ", length(rb), "\n")
#		cat("dim rv = ", length(rv), "\n")
		out <- list(var=rv, bias=rb, mse=rv+rb^2)
	} else {
		if (is.null(XtXlamIinv)){
			XtXlamIinv <- lapply(lambda, getXtXlamIinv, XtX=XtX)
		}
		rv <- sapply(XtXlamIinv, function(w) (w %*% XtX %*% w )[ind, ind])
#		rv <- sapply(lambda, ridgeVar, XtX=XtX, ind=ind)
		rv <- tcrossprod(sigma2, rv)
		lamDiag <- lapply(lambda, function(w) diag(w*c(0, rep(1, ncol(XtX)-1))))
		XtXlamIinvMX  <-  mapply(function(x, y) (x %*% y)[ind, ], XtXlamIinv, lamDiag)
#		XtXlamIinvMX <- simplify2array(XtXlamIinv)[ind, , ]
		rb <- crossprod(as.matrix(beta), XtXlamIinvMX)
#		cat("dim rb = ", dim(rb), "\n")
#		cat("dim rv = ", dim(rv), "\n")
#		rb <- sapply(lambda, ridgeBias, XtX=XtX, beta=beta, ind=ind)
		out <- list(var=rv, bias=rb, mse=rv+rb^2)
	}
	return(out)
}



getXtXlamIinv <- function (lambda, XtX, penalize) {
    p <- ncol(XtX)
    if (missing(penalize)) 
        penalize <- c(0, rep(1, p - 1))
    penalize <- as.numeric(penalize)
    if (max(abs(penalize)) > 1) 
        stop("Element of penalize must be <= 1.")
    if (length(lambda) > 1) {
        lambda <- lambda[1]
        warning("lambda should be of length 1. Only using first value.")
    }
    penalty <- lambda * penalize
	XtXlamIinv <- chol2inv(chol(XtX + diag(penalty))) 
    return(XtXlamIinv)
}




##' @name vcovRidgefEst
##' @title Standard Error Estimate
##
##' @description Computes an estimate of the variance-covariance
##'		matrix for an 'fLoss' ridge estimator
##
##' @details Computes a standard error estimate for an 'fLoss'
##'		estimator, where 'fLoss' is typically fMSE or fMBV.
##'		Approximates the variance of the estimator using the 
##'		the variance conditional on the observed data (i.e. 
##'		using the posterior distribution of parameters).
##'     Currently, two different versions are available.
##
##  Input:
##' @param fLoss A matrix of loss function values to be minimized.
##'		Assumed structure is the columns correspond to different
##'		values of penalty parameter and rows correspond to points
##'		in a posterior sample of (beta, sigma).
##'	@param lambda The sequence of penalty parameter values
##'		corresponding to the columns of \code{fLoss}.
##'	@param XtX Cross product of the design matrix.
##' @param postBeta Matrix containing the posterior sample of 
##'		beta values. Assumed to be n-by-p, where n is number of
##'		samples (and equal to number of rows in fLoss) and p is
##'		number of regression parameters in model.
##' @param postSigma2 Vector containing the posterior sample
##'		of variance parameters.  Should have same length as 
##'		postBeta.
##'	@param penalize Vector indicating which variables are
##'		penalized in the regression model.  See details for
##'		\code{\link{ridgeEst}} for further details.
##' @param ind Numerical or logical vector indicating which elements
##'		of the variance matrix should be returned. Defaults to the
##'		(1,1) element
##'	@param version Character string indicating which version of
##'		standard error to compute. 'varExp' or 'full', corresponding
##'		to the variance of the conditional mean of the estimator or
##'		or that plus the expected value of the conditional variance,
##'		which gives the full posterior variance of the estimator.
##' @importFrom stats var
##' @export
##' @seealso \code{\link{ridgefEst}, \link{samplePosterior}}
vcovRidgefEst <- function(fLoss, lambda, XtX, postBeta, postSigma2, penalize, ind=1, version=c("varExp", "full")){
	version <- match.arg(version)
	p <- ncol(XtX)
	if(missing(penalize)) penalize <- c(0, rep(1,p-1))
	penalize <- as.numeric(penalize)
	if (max(abs(penalize))>1) stop("Element of penalize must be <= 1.")
	
	# Dimension checks
	n <- nrow(fLoss)
	nlambda <- ncol(fLoss)
	if (length(lambda)!=nlambda) stop("Length of lambda sequence must be equal to ncol(fLoss).")
	if (is.vector(postBeta)) postBeta <- matrix(postBeta, ncol=1)
	if (nrow(postBeta)!=n) stop("Number of beta samples should equal nrow(fLoss).")
	
	lmintemp <- apply(fLoss, 1, which.min)
	condMean <- apply(cbind(lambda[lmintemp], postBeta), 1, function(w) (solve(XtX + diag(w[1]*penalize), XtX %*%  w[-1]))[ind])
	if (version=="varExp"){
		if (is.matrix(condMean)) {
			varmx <- var(t(condMean))
		} else {
			varmx <- var(condMean)
		}
	} else if (version=="full") {
		condVar <- mapply(ridgeVar, lambda=lambda[lmintemp], sigma2=postSigma2, XtX=list(XtX), ind=list(ind))
		if (is.matrix(condMean)) {
		varmx <- var(t(condMean)) + matrix(rowMeans(condVar), ncol=sqrt(nrow(condVar)))
		} else {
		varmx <- var(condMean) + mean(condVar)
		}
	}
	varmx
}



##' @name ridgefEst
##' @title Compute `Future MSE' Ridge Estimate
##
##' @description Computes a ridge estimate for a given regression problem, with penalty parameter chosen to minimize future MSE.
##
##' @details To be added.....
##'	
##'		To balance the influence of covariates, it is recommended
##'		that the design matrix be standardized.  This can be done by
##'		the user or via the argument \code{scale}.  If \code{scale=TRUE},
##'		then coefficent and standard error estimates are back-transformed.
##'
##  Input:
##'	@param X Design matrix of regression problem
##' @param y outcome vector. Typically centered.
##' @param type The type of estimator. See details.
##' @param lseq Sequence of penalty values to consider.
##' @param scale Logical indicating whether the design matrix X be scaled. See details.
##'	@param returnMSE Logical indicating whether mse object should be returned.
##' @param postsamp List containing posterior sample (from \code{samplePosterior}). If
##'		missing, then a posterior sample is drawn.  Currently checks on the provided
##'		\code{postsamp} are limited, so use with caution.  Designed to facilitate
##'		simualtions or other scenarios where it may be pre-computed.
##' @param returnPS logical indicating whether or not the full posterior sample should
##'		be included in output.
##' @param nPost Size of posterior sample to compute
##' @param se logical indicating whether or not standard errors should be returned.
##' @param XtXlamIinv explicit value of (X'X + diag(penalty))^{-1}.  Useful
##'		for simulations to save computation.
##' @param ... Other arguments passed to \code{samplePosterior}
##' @export
##' @seealso \code{\link{lassofEst}}, \code{\link{ridgeMSE}}
ridgefEst <- function(X, y, type=c("fMSE", "fMBV", "all"), lseq, scale=FALSE, returnMSE=FALSE, postsamp, returnPS=FALSE, nPost=1000, se=TRUE, XtXlamIinv =NULL, ...){
	type <- match.arg(type)

	if (scale) X <- scale(X)	
	if(missing(lseq)) lseq <-  exp(seq(log(0.01), log(1000), length=500))
	
	if (missing(postsamp)) postsamp <- samplePosterior(X, y, n=nPost,...)

	XtX <- crossprod(X)
	rmse <- ridgeMSE(lseq, XtX, beta= t(postsamp$beta), sigma2= postsamp $sigma2, XtXlamIinv = XtXlamIinv)

	# Get estimator
	if(type=="fMSE"){
		# Compute posterior means of MSE
		pmse <- colMeans(rmse$mse)
		lmin <- which.min(pmse)
		beta <- ridgeEst(lseq[lmin], X=X, y=y)
		if (se) {
			se <- sqrt(diag(vcovRidgefEst(fLoss=rmse$mse, lambda=lseq, XtX=XtX, postBeta=postsamp$beta, postSigma2=postsamp$sigma2, ind=rep(TRUE, length(beta)))))
		} else {
			se <- rep(0, length(beta))
		}
		if (scale) {
			beta <- beta/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"
		}
		out <- list(beta=beta, lambda=lseq[lmin], lmin=lmin, se=se)	
	} else if (type=="fMBV"){
#		pmbv <- colMeans(pmax(rmse$bias^2, rmse$var))
		rb2 <- rmse$bias^2
		pmbv <- colMeans( (rb2 +  rmse$var + abs(rb2- rmse$var))/2)
		lmin <- which.min(pmbv)
		beta <- ridgeEst(lseq[lmin], X=X, y=y)
		if (se){
			se <- sqrt(diag(vcovRidgefEst(fLoss=pmax(rmse$bias^2, rmse$var), lambda=lseq, XtX=XtX, postBeta=postsamp$beta, postSigma2=postsamp$sigma2, ind=rep(TRUE, length(beta)))))		
		} else {
			se <- rep(0, length(beta))
		}
		if (scale) {
			beta <- beta/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"
		}
		out <- list(beta=beta, lambda=lseq[lmin], lmin=lmin,  se=se)
	} else if (type=="all"){
#		pmbv <- colMeans(pmax(rmse$bias^2, rmse$var))
		rb2 <- rmse$bias^2
		pmbv <- colMeans( (rb2 +  rmse$var + abs(rb2- rmse$var))/2)
		lminfMBV <- which.min(pmbv)
		betafMBV <- ridgeEst(lseq[lminfMBV], X=X, y=y)
		pmse <- colMeans(rmse$mse)
		lminfMSE <- which.min(pmse)
		betafMSE <- ridgeEst(lseq[lminfMSE], X=X, y=y)
		if (se){
			sefMBV <- sqrt(diag(vcovRidgefEst(fLoss=pmax(rmse$bias^2, rmse$var), lambda=lseq, XtX=XtX, postBeta=postsamp$beta, postSigma2=postsamp$sigma2, ind=rep(TRUE, length(betafMBV)))))
			sefMSE <- sqrt(diag(vcovRidgefEst(fLoss=rmse$mse, lambda=lseq, XtX=XtX, postBeta=postsamp$beta, postSigma2=postsamp$sigma2, ind=rep(TRUE, length(betafMSE)))))
		} else {
			sefMBV <- rep(0, length(betafMBV))		
			sefMSE <- rep(0, length(betafMSE))		
		}
		if (scale) {
			betafMSE <- betafMSE/attributes(X)$"scaled:scale"
			betafMBV <- betafMBV/attributes(X)$"scaled:scale"
			sefMBV <- sefMBV/attributes(X)$"scaled:scale"
			sefMSE <- sefMSE/attributes(X)$"scaled:scale"
		}
		se <- cbind(fMSE= sefMSE, fMBV= sefMBV) 

		out <- list(beta=cbind(fMSE= betafMSE, fMBV= betafMBV), lambda=c(fMSE=lseq[lminfMSE], fMBV=lseq[lminfMBV]), lmin=c(fMSE = lminfMSE, fMBV = lminfMBV), se=se)
	}

	if (returnMSE){
		out <- list(out, MSE=rmse$mse, Bias2=rmse$bias^2, Var=rmse$var)
	}
	if (returnPS){
		out$postsamp <- postsamp
	}
	return(out)
}


