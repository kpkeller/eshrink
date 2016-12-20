#####################################
#####################################
# Functions for the Ridge Estimator
#
# Included in this file:
#
#	simLassoEst()
#	lassofutureEst()
#
#####################################




##' @name simLassoEst
##' @title Compute Lasso Estimator for simulated Data
##
##' @description Simulates data from a regression model
##'		and computes the lasso estimate for this data.
##
##' @details Simulates data from a regression model with true
##'		coefficient parameter \code{beta} and normal errors with
##'		standard deviation \code{sigma}.  Computes the LASSO
##'		estimate for the coefficient vector using the \code{glmnet}
##'		function from the package of the same name.
##
##  Input:
##' @param lambda Penalty factor to be applied
##'	@param X Design matrix of regression problem
##' @param beta true value of parameter vector to simulate from
##' @param sigma true value of square root of variance parameter for simulating.
##' @param penalize Vector giving penalty structure. Supplied to glmnet as `\code{penalty.factor}'. By default, all coefficients except first are penalized.
##'	@param rescale.lambda Should lambda be rescaled to account for the 
##'		default re-scaling done by glmnet?
##'	@param ind Index of coefficient to be returned.  Value of 0 implies
##'		all coefficients (i.e. the full parameter vector estimate)
##
##'	@author Joshua Keller
##' @export
##' @import glmnet
##' @importFrom stats rnorm
##
# Gets (sequence of) Lasso estimator for a given (sequence of) lambda
# values, by simulating new observations
#
# Note: Current hard coding of lack of intercept
# and that first column is variable of interest
# Also hard coded right now is that all but
# exposure of interest are penalized
simLassoEst <- function(lambda, X,  beta, sigma, penalize, rescale.lambda=TRUE, ind=1){
	n <- nrow(X)
	p <- ncol(X)
	
	if (ind==0){
		ind <- 1:length(beta)
	}
	
	if(missing(penalize)) penalize <- c(0, rep(1,p-1))
	penalize <- as.numeric(penalize)
	if (any(penalize>1) | any(penalize<0)) stop("Element of penalize must be between 0 and 1.")

	eps <- stats::rnorm(n, sd=sigma)
	y <- X %*% beta + eps	
	
	# Rescale lambda to account for re-scaling done by glmnet
	if (rescale.lambda){
		if(is.logical(rescale.lambda)){
			plam <- sum(penalize==0)
			lambda <- lambda*(p-plam)/p
		} else {
			lambda <- lambda*rescale.lambda
		}
	}
	gfit <- glmnet::glmnet(x= X, y=y, lambda= lambda,  penalty.factor=penalize, standardize=FALSE, intercept=FALSE, family="gaussian")
	lassoEsts <- gfit$beta[ind,]
	return(lassoEsts)
}


	

##' @name lassofEst
##' @title Compute Lasso Estimator for simulated Data
##
##' @description To be added...
##
##' @details To be added...
##
##  Input:
##' @param lambda Penalty factor to be applied
##'	@param X Design matrix of regression problem
##' @param y outcome vector. Typically centered.
##' @param type Which type of estimate should be computed? See details.
##' @param penalize Vector giving penalty structure. Supplied to glmnet as `\code{penalty.factor}'.
##'	@param rescale.lambda Should lambda be rescaled to account for the 
##'		default re-scaling done by glmnet?
##' @param scale Logical indicating whether the design matrix X be scaled. See details.
##'	@param returnMSE Logical indicating whether mse object should be returned.
##' @param lseq Sequence of penalty values to consider.
##' @param B Number of future datasets to simulate for each point in posterior sample.
##'	@param se.version Which version of Standard errors to report?
##' @param postsamp List containing posterior sample (from \code{samplePosterior}). If
##'		missing, then a posterior sample is drawn.  Currently checks on the provided
##'		\code{postsamp} are limited, so use with caution.  Designed to facilitate
##'		simualtions or other scenarios where it may be pre-computed.
##' @param returnPS logical indicating whether or not the full posterior sample should
##'		be included in output.
##' @param nPost Size of posterior sample to compute
##' @param ... Other arguments passed to \code{samplePosterior}
##'
##' @export
##' @import glmnet
##
##' @seealso \code{\link{ridgefEst}}, \code{\link{simLassoEst}}
##
# Need to make penalty factor code
# more general
lassofEst <- function(X, y, type=c("fMSE", "fMBV", "all"), lseq, B=500, penalize=NULL, rescale.lambda=TRUE, scale=FALSE, returnMSE=FALSE, se.version=c("varExp", "full", "both"), postsamp, returnPS=FALSE, nPost=1000, ...){
	
	fncall <- match.call()
	type <- match.arg(type)
	se.version <- match.arg(se.version)	
		
	if(scale) X <- scale(X)
	if(missing(lseq)) lseq <-  exp(seq(log(100), log(0.01), length=500))

	# Sort lseq in decreasing order
	lseq <- rev(sort(lseq))

	XtX <- crossprod(X)
	p <- ncol(X)
	
	if(missing(penalize)) penalize <- c(0, rep(1,p-1))
	penalize <- as.numeric(penalize)
	if (any(penalize>1) | any(penalize<0)) stop("Element of penalize must be between 0 and 1.")
	
	# Rescale lambda to account for re-scaling done by glmnet
	if (rescale.lambda){
		if(is.logical(rescale.lambda)){
			plam <- sum(penalize==0)
			lseq <- lseq*(p-plam)/p
		} else {
			lseq <- lseq*rescale.lambda
		}
	}
	
	# Make Posterior Draw
	if(missing(postsamp)) postsamp <- samplePosterior(X, y, n= nPost, ...)

	Bias2 <- matrix(NA, nrow=nPost, ncol=length(lseq))
	MSE <- Var <- Bias2
	lassoEsts <- array(NA, dim=c(B, length(lseq)))
	lassoEstVarMBV <- lassoEstVarMSE <- lassoEstMeanMBV  <- lassoEstMeanMSE <- numeric(nPost)
	for (bg in 1:nPost){
	for (i in 1:B){
		set.seed(50005 + i)
		lassoEsts[i, ] <- simLassoEst(lambda=lseq, X=X, beta= postsamp $beta[bg,], sigma=sqrt(postsamp $sigma2[bg]), penalize= penalize, rescale.lambda=FALSE) # DOn't rescale, because already rescaled
	}
	# Note this hard coding of 'X' here.
	# should make more general
	Bias2[bg,] <-  apply(lassoEsts, 2, getBias2, truth= postsamp $beta[bg,"X"])
	Var[bg,] <- apply(lassoEsts, 2, getPopVar)
	MSE[bg, ] <- Bias2[bg,] + Var[bg,]
	MBV <- pmax(Bias2[bg,], Var[bg,])
	if (type =="fMSE") {
		lassoEstMeanMSE[bg] <- mean(lassoEsts[,which.min(MSE[bg, ])])
		lassoEstVarMSE[bg] <- var(lassoEsts[,which.min(MSE[bg, ])])
	} else if (type =="fMBV"){
		lassoEstMeanMBV[bg] <- mean(lassoEsts[,which.min(MBV)])
		lassoEstVarMBV[bg] <- var(lassoEsts[,which.min(MBV)])		
	} else if (type =="all"){
		lassoEstMeanMSE[bg] <- mean(lassoEsts[,which.min(MSE[bg, ])])
		lassoEstMeanMBV[bg] <- mean(lassoEsts[,which.min(MBV)])
		lassoEstVarMSE[bg] <- var(lassoEsts[,which.min(MSE[bg, ])])
		lassoEstVarMBV[bg] <- var(lassoEsts[,which.min(MBV)])
	}
	}

	# Get estimator
	if(type=="fMSE"){
		pmse <- colMeans(MSE)
		lmin <- which.min(pmse)
		gfit <- glmnet::glmnet(x= X, y= y, lambda= lseq,  penalty.factor= penalize, standardize=FALSE, intercept=FALSE, family="gaussian")
		beta <- gfit$beta[,lmin]
		if (se.version=="varExp") {
			se <- sqrt(var(lassoEstMeanMSE))
		} else if (se.version=="full") {
			se <- sqrt(var(lassoEstMeanMSE) + mean(lassoEstVarMSE))
		} else if (se.version=="both") {
			se <- c(varExp= sqrt(var(lassoEstMeanMSE)), full=sqrt(var(lassoEstMeanMSE) + mean(lassoEstVarMSE)))
		}
		if (scale) {
			beta <- beta/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"[1]
		}
		out <- list(beta=beta, lambda=lseq[lmin], lmin=lmin, se=se)
	} else if (type=="fMBV"){
		pmaxBias2Var <- colMeans(pmax(Bias2, Var))
		lmin <- which.min(pmaxBias2Var)
		gfit <- glmnet::glmnet(x= X, y= y, lambda= lseq,  penalty.factor= penalize, standardize=FALSE, intercept=FALSE, family="gaussian")
		beta <- gfit$beta[,lmin]
		if (se.version=="varExp") {
			se <- sqrt(var(lassoEstMeanMBV))
		} else if (se.version=="full") {
			se <- sqrt(var(lassoEstMeanMBV) + mean(lassoEstVarMBV))
		} else if (se.version=="both") {
			se <- c(varExp= sqrt(var(lassoEstMeanMBV)), full=sqrt(var(lassoEstMeanMBV) + mean(lassoEstVarMBV)))
		}
		if (scale) {
			beta <- beta/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"[1]
		}
		out <- list(beta=beta, lambda=lseq[lmin], lmin=lmin, se=se)
	} else if (type=="all"){
		pmaxBias2Var <- colMeans(pmax(Bias2, Var))
		lminBias2Var <- which.min(pmaxBias2Var)
		gfit <- glmnet::glmnet(x= X, y= y, lambda= lseq,  penalty.factor= penalize, standardize=FALSE, intercept=FALSE, family="gaussian")
		betaMaxBias2var <- gfit$beta[, lminBias2Var]
		pmse <- colMeans(MSE)
		lmin <- which.min(pmse)
		betaMSE <- gfit$beta[, lmin]
		if (se.version=="varExp") {
			se <- c(fMSE=sqrt(var(lassoEstMeanMSE)), fMBV=sqrt(var(lassoEstMeanMBV)))
		} else if (se.version=="full") {
			se <- c(fMSE=sqrt(var(lassoEstMeanMSE) + mean(lassoEstVarMSE)), fMBV=sqrt(var(lassoEstMeanMBV) + mean(lassoEstVarMBV)))
		} else if (se.version=="both") {
		se <- c(fMSE_varExp=sqrt(var(lassoEstMeanMSE)), fMSE_full=sqrt(var(lassoEstMeanMSE) + mean(lassoEstVarMSE)), fMBV_varExp=sqrt(var(lassoEstMeanMBV)), fMBV_full=sqrt(var(lassoEstMeanMBV) + mean(lassoEstVarMBV)))
		}				
		if (scale) {
			betaMSE <- betaMSE/attributes(X)$"scaled:scale"
			betaMaxBias2var <- betaMaxBias2var/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"[1]
		}
		out <- list(beta=cbind(fMSE=betaMSE, fMBV= betaMaxBias2var), lambda=c(fMSE =lseq[lmin], fMBV=lseq[lminBias2Var]), lmin=c(fMSE =lmin, fMBV =lminBias2Var), se=se, call=fncall)
	}

	if (returnMSE){
		out <- list(out, MSE=MSE, Bias2=Bias2, Var=Var)
	}
	if (returnPS){
		out$postsamp <- postsamp
	}
	return(out)
}


