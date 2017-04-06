#####################################
#####################################
#####################################
# Functions for the Ridge Estimator
#
# Included in this file:
#
#	estRidge()
#	biasRidge()
#	varRidge()
#	mseRidge()
#	vcovfestRidge()
#	festRidge()
#


##' @name estRidge
##' @title Estimate Coefficients for Ridge Regression 
##
##' @description Computes a vector of regression coefficients for a provided ridge penalty.
##
##
##' @details 
##' The input \code{penalize} is a vector of ridge penalty factors, 
##'	such that the penalty for covariate j is \code{lambda*penalize[j]}.
##' Although its primary purpose is for indicating which variables to penalize (1)
##' and which to not penalize (0), fractional values between 0 and 1 are accepted.
##' Defaults to c(0, rep(1, p-1)), where
##' p is number of columns in X (this penalizes all coefficients but 
##' the first).  
##' 
##' The design matrix \code{X} is assumed to contain only numeric values, so 
##'			any factors should be coded according to desired contrast (e.g., via \code{\link{model.matrix}})
##'
##' 
##
##  Input:
##' @param lambda ridge penalty factor 
##' @param X design matrix for the regression. 
##' @param y outcome vector. Unless \code{X} contains an intercept column, this should typically be centered.
##' @param penalize vector giving penalty structure. Values must be in [0, 1].
##'	See Details for more information.
##' @param XtX (optional) cross product of the design matrix. If running simulations or 
##'		other procedure for identical \code{X}, providing a pre-computed value 
##'		can reduce computational cost.
##
##' @export
##' @author Joshua Keller
##' @seealso \code{\link{festRidge}}, \code{\link{mseRidge}}
estRidge <- function(lambda, X, y,  penalize, XtX=crossprod(X)){
	if(missing(penalize)) penalize <- c(0, rep(1, ncol(XtX)-1))
	penalize <- as.numeric(penalize)
	if (any(penalize>1) | any(penalize<0)) stop("Element of 'penalize' must be between 0 and 1.")
	penalty <- lambda*penalize
	as.vector(solve(XtX + diag(penalty), crossprod(X, y)))
}



##' @name mseRidge
##' @title Compute MSE, Bias, and Variance for Ridge Estimator
##
##' @description Computes the analytic mean-squared error (MSE), bias, and 
##'		variance for ridge regression estimators given different
##'		values of the true \code{beta} and \code{sigma2} parameters.
##
##' @details 
##' 	The computations assume that all covariates are correctly included in the
##'		mean model and bias is due only to penalization. The bias is given by:
##'
##'			\eqn{-(X'X + \Lambda)^{-1}\Lambda\beta}
##'
##'		where \eqn{\Lambda = diag(\lambda*penalize)}. The variance is given by: 
##'
##'			\eqn{\sigma^2(X'X + \Lambda)^{-1}X'X(X'X + \Lambda)^{-1}}
##'		
##' 	If \code{beta} is provided as a matrix, this will treat
##' 	each column of \code{beta} as a different true parameter vector 
##' 	and return a matrix of bias values (or a vector, if \code{ind} has length 1).
##'
##' 	Providing a pre-computed value of \code{XtXlamIinv} can reduce the computational
##'		cost in simulations. However, the user is responsible for assuring that the value
##' 	of \code{lambda} provided matches the value used to compute \code{XtXlamIinv}.
##
##
##  Input:
##' @param lambda penalty parameter value. For \code{biasRidge} and \code{varRidge}, this should be 
##'		a single value. For \code{mseRidge}, either a single value of a list of values.
##' @param XtX Cross product of design matrix. Not needed if \code{XtXlamIinv} is provided.
##' @param penalize Vector of penalty factors. See \code{\link{estRidge}} for more information.
##' @param beta True parameter values. Either a vector of length \code{p} or a 
##'			\code{p} x \code{d} matrix.
##' @param sigma2 Value of the variance parameter
##' @param ind Numerical or logical vector indicating which elements of the bias vector and
##'		 variance matrix should be returned. Defaults to the first element.
##' @param XtXlamIinv Optional explicit value of \code{(XtX + diag(lambda*penalize))^(-1)}.
##	
## Output:
##'	@return  For \code{mseRidge}, a list containing the variance, bias, and MSE. For \code{biasRidge} and \code{varRidge}, a matrix is returned.
##
##'
##' @export
##' @author Joshua Keller
mseRidge <- function(lambda, XtX, beta, sigma2, penalize, ind=1, XtXlamIinv=NULL){
	if (is.logical(ind)){
		if (sum(ind)>1) stop("'ind' should be a single integer or logical vector with sum 1.")
	} else {
		if (length(ind)>1)stop("'ind' should be a single integer or logical vector with sum 1.")
	}
	if(missing(penalize)) penalize <- c(0, rep(1, ncol(XtX)-1))
	penalize <- as.numeric(penalize)
	if (max(abs(penalize))>1) stop("Element of penalize must be <= 1.")

	if (is.matrix(beta) && nrow(beta)!=ncol(XtX)) stop("Beta should be vector, or matrix with ncol(XtX) rows.")

	if (length(lambda)==1){
		rv <- varRidge(lambda= lambda, XtX=XtX, penalize=penalize,  ind=ind)
		rv <- as.vector(tcrossprod(sigma2, rv))
		rb <- as.vector(biasRidge(lambda, XtX, beta=beta, penalize=penalize, ind=ind))
		out <- list(var=rv, bias=rb, mse=rv+rb^2)
	} else {
		if (is.null(XtXlamIinv)){
			XtXlamIinv <- lapply(lambda, getXtXlamIinv, XtX=XtX, penalize=penalize)
		}
		rv <- sapply(XtXlamIinv, function(w) sum(diag(w %*% XtX %*% w )[ind]))
#		rv <- sapply(lambda, varRidge, XtX=XtX, ind=ind)
		rv <- tcrossprod(sigma2, rv)
		lamDiag <- lapply(lambda, function(w) diag(w*penalize))
		XtXlamIinvMX  <-  mapply(function(x, y) crossprod(x,y)[ind, ], XtXlamIinv, lamDiag)
#		XtXlamIinvMX <- simplify2array(XtXlamIinv)[ind, , ]
		rb <- crossprod(as.matrix(beta), XtXlamIinvMX)
#		rb <- sapply(lambda, biasRidge, XtX=XtX, beta=beta, ind=ind)
		out <- list(var=rv, bias=rb, mse=rv+rb^2)
	}
	return(out)
}



##' @name biasRidge
##' @rdname mseRidge
##
## NOTE: Merged into mseRidge
##
##
##' @export
## @author Joshua Keller
##	@seealso \code{\link{estRidge}}
biasRidge <- function (lambda, XtX, beta, penalize, ind = 1, XtXlamIinv=NULL) {
    if (is.null(XtXlamIinv)){
    		    	p <- ncol(XtX)
    } else {
    	  	    p <- ncol(XtXlamIinv)	
    }
    
    	if(missing(penalize)) penalize <- c(0, rep(1, p-1))
	
	penalize <- as.numeric(penalize)
    if (max(abs(penalize)) > 1) 
        stop("Element of penalize must be <= 1.")
    if (length(lambda) > 1) {
        lambda <- lambda[1]
        warning("lambda should be of length 1. Only using first value.")
    }
    penalty <- lambda * penalize
    
    if (is.null(XtXlamIinv)){
	    if(nrow(XtX)!=p) stop("XtX must be square.")
	    XtXlamIinv <- chol2inv(chol(XtX + diag(penalty)))
    }
    if (is.matrix(beta) && nrow(beta) != p) 
        stop("Beta should be vector, or matrix with ncol(XtX) rows.")
	XtXlamIinvB <- -crossprod(penalty*t(XtXlamIinv[ind, , drop=F]),  beta)
    return(XtXlamIinvB)
}


##' @name varRidge
##' @rdname mseRidge
##
## @description Computes the variance matrix of a ridge estimator, conditional on the penalty and global variance \code{sigma2}. Specifically, this function computes:
##		where \eqn{\Lambda = diag(\lambda*penalize)}.
##
## @details See \code{\link{estRidge}} for description of \code{penalize} parameter.
##
##  Input:
## @param XtX Cross product of design matrix
## @param lambda Penalty parameter.
##	@param penalize Vector giving what to penalize
## @param ind Numerical or logical vector indicating which elements
##		of the variance matrix diagonal should be returned. Defaults to the
##		(1,1) element
## @param XtXlamIinv Explicit value of \eqn{(X'X + diag(penalty))^(-1)},
##		where \code{penalty=lambda*penalize}. Useful
##		for simulations to save computation.
##' @export
## @author Joshua Keller
## @seealso \code{\link{mseRidge}, \link{biasRidge}}
varRidge <- function(lambda, XtX, sigma2=1, penalize, ind=1, XtXlamIinv=NULL){
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


# Helper function
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




##' @name vcovfestRidge
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
##'		\code{\link{estRidge}} for further details.
##' @param ind Numerical or logical vector indicating which elements
##'		of the variance matrix should be returned. Defaults to the
##'		(1,1) element
##'	@param version Character string indicating which version of
##'		standard error to compute. 'varExp' or 'full', corresponding
##'		to the variance of the conditional mean of the estimator or
##'		that plus the expected value of the conditional variance. In
##'		practice, the latter is often too large.
##
##' @importFrom stats var
##' @export
##' @author Joshua Keller
##' @seealso \code{\link{festRidge}, \link{samplePosterior}}
vcovfestRidge <- function(fLoss, lambda, XtX, postBeta, postSigma2, penalize, ind=1, version=c("varExp", "full")){
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
		condVar <- mapply(varRidge, lambda=lambda[lmintemp], sigma2=postSigma2, XtX=list(XtX), ind=list(ind))
		if (is.matrix(condMean)) {
		varmx <- var(t(condMean)) + matrix(rowMeans(condVar), ncol=sqrt(nrow(condVar)))
		} else {
		varmx <- var(condMean) + mean(condVar)
		}
	}
	varmx
}



##' @name festRidge
##' @rdname festRidge
##' @title Compute `Future Loss' Ridge or LASSO Estimates 
##
##' @description Computes a ridge or LASSO estimate for a given regression problem, with penalty parameter chosen to minimize bias and variance.
##
##' @details 
##'		The value of the ridge or LASSO penalty is selected by minimizing the 
##'		posterior expectation of the loss function, which is chosen by the argument
##'		\code{loss}. Possible options are \code{fMBV}, which uses the loss function
##'		\eqn{fMBV = max(Bias(\beta)^2, Var(\beta))} and \code{fMSE}, which uses the loss function
##'		\eqn{fMSE = Bias(\beta)^2 +  Var(\beta)}. 
##'	
##'		To balance the influence of covariates, it is recommended
##'		that the design matrix be standardized.  This can be done by
##'		the user or via the argument \code{scale}.  If \code{scale=TRUE},
##'		then coefficient and standard error estimates are back-transformed.
##'
##'		Use the \code{XtXlamIinv} argument with caution. No checks are done on the provided
##'		value. Note that \code{lseq} is re-ordered to be decreasing, and provided values
##'		of \code{XtXlamIinv} must account for this.
##'
##  Input:
##' @param X Design matrix for the regression. Assumed to contain only numeric values, so 
##'			any factors should be coded according to desired contrast (e.g., via \code{\link{model.matrix}})
##' @param y Outcome vector. Unless \code{X} contains an intercept column, this should typically be centered.
##' @param loss Loss function for choosing the penalty parameter. See details.
##'	@param ind Vector of integers or logicals indicating which coefficients the loss is to be computed on.
##' @param lseq Sequence of penalty values to consider.
##' @param penalize See \code{\link{estRidge}}
##' @param scale Logical indicating whether the design matrix X be scaled. See details.
##'	@param returnMSE Logical indicating whether mse object should be returned.
##' @param postsamp List containing posterior sample (from \code{samplePosterior}). If
##'		missing, then a posterior sample is drawn.  Currently checks on the provided
##'		\code{postsamp} are limited, so use with caution.  Designed to facilitate
##'		simulations or other scenarios where it may be pre-computed.
##' @param returnPS logical indicating whether or not the full posterior sample should
##'		be included in output.
##' @param nPost Size of posterior sample to compute
##' @param se.version String indicating which version of standard errors to use. See \code{\link{vcovfestRidge}}.
##' @param XtXlamIinv explicit value of (X'X + diag(penalty))^{-1}.  Useful
##'		for simulations to save computation. 
##' @param ... Other arguments passed to \code{samplePosterior}
##' @export
##' @seealso  \code{\link{mseRidge},\link{vcovfestRidge}, \link{simLASSO}, \link{check_CIbound}}
## \code{\link{festLASSO}}
festRidge <- function(X, y, loss=c("fMSE", "fMBV", "both"), ind=1, lseq, penalize, scale=FALSE, returnMSE=FALSE, postsamp, returnPS=FALSE, nPost=1000, se.version =c("varExp", "full", "none"), XtXlamIinv=NULL, ...){

	fncall <- match.call()
	loss <- match.arg(loss)
	se.version <- match.arg(se.version)	

	if (scale) X <- scale(X)	
	if (missing(lseq)) lseq <-  exp(seq(log(1000), log(0.01), length=500))

	# Sort lseq in decreasing order
	lseq <- rev(sort(lseq))

	XtX <- crossprod(X)
	p <- ncol(X)

	if(missing(penalize)) {
		penalize <- rep(1, p)
		penalize[ind] <- 0
		#c(0, rep(1,p-1))
	}
	penalize <- as.numeric(penalize)
	if (any(penalize>1) | any(penalize<0)) stop("Element of penalize must be between 0 and 1.")
	
	# Make Posterior Draw, if missing
	if (missing(postsamp)) postsamp <- samplePosterior(X, y, n=nPost,...)

	# Compute Bias2 and Variance
	rmse <- mseRidge(lseq, XtX, beta= t(postsamp$beta), sigma2= postsamp $sigma2, XtXlamIinv = XtXlamIinv, penalize=penalize, ind=ind)

	out <- list()
	# Get estimator
	if(loss=="fMSE" | loss =="both"){
		# Compute posterior means of MSE
		pmse <- colMeans(rmse$mse)
		lmin <- which.min(pmse)
		beta <- estRidge(lseq[lmin], X=X, y=y, penalize=penalize)
		if (se.version =="none"){
			se <- rep(0, length(beta))
		} else {
			se <- sqrt(diag(vcovfestRidge(fLoss=rmse$mse, lambda=lseq, XtX=XtX, postBeta=postsamp$beta, postSigma2=postsamp$sigma2, ind=rep(TRUE, length(beta)),version= se.version)))
		}
		if (scale) {
			beta <- beta/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"
		}
		out$fMSE <- list(beta=beta, lambda=lseq[lmin], lmin=lmin, se=se)	
	}
	if (loss=="fMBV" | loss =="both"){
#		pmbv <- colMeans(pmax(rmse$bias^2, rmse$var))
		rb2 <- rmse$bias^2
		mbv <- (rb2 +  rmse$var + abs(rb2- rmse$var))/2
		pmbv <- colMeans(mbv )
		lmin <- which.min(pmbv)
		beta <- estRidge(lseq[lmin], X=X, y=y, penalize=penalize)
		if (se.version =="none"){
			se <- rep(0, length(beta))		
		} else {
			se <- sqrt(diag(vcovfestRidge(fLoss= mbv, lambda=lseq, XtX=XtX, postBeta=postsamp$beta, postSigma2=postsamp$sigma2, ind=rep(TRUE, length(beta)),version= se.version)))
		}
		if (scale) {
			beta <- beta/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"
		}
		out$fMBV <- list(beta=beta, lambda=lseq[lmin], lmin=lmin,  se=se)
	}
	
	if (returnMSE){
		out$MSE <- rmse$mse
		out$Bias2 <- rmse$bias^2
		out$Var <- rmse$var
	}
	if (returnPS){
		out$postsamp <- postsamp
	}
	return(out)
}


