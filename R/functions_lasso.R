#####################################
#####################################
# Functions for the Ridge Estimator
#
# Included in this file:
#
#	simLASSO()
#	festLASSO()
#
#####################################




##' @name simLASSO
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
simLASSO <- function(lambda, X,  beta, sigma, penalize, rescale.lambda=TRUE, ind=1){
	n <- nrow(X)
	p <- ncol(X)
	
	if (sum(abs(ind))==0){
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



mseLASSO <- function(lambda, X, postsamp, nPost, B, penalize, ind){
	Bias2 <- matrix(NA, nrow=nPost, ncol=length(lambda))
	MBV <- MSE <- Var <- Bias2
	coefsLASSO <- array(NA, dim=c(B, ncol(X), length(lambda)))
	estMeanMSE <- matrix(NA, nrow=nPost, ncol=ncol(X))
	estMeanMBV <- estMeanMSE
	estVarMSE <- array(NA, dim=list(nPost, ncol(X), ncol(X)))
	estVarMBV <- estVarMSE
	for (bg in 1:nPost){
	for (i in 1:B){
		set.seed(50005 + i)
		coefsLASSO[i, , ]  <- as.matrix(simLASSO(lambda=lambda, X=X, beta= postsamp$beta[bg,], sigma=sqrt(postsamp$sigma2[bg]), penalize=penalize, rescale.lambda=FALSE, ind=1:ncol(X))) # Don't rescale, because already rescaled
	}
	Bias2[bg,] <-  colSums((apply(coefsLASSO[, ind, , drop=FALSE], 2:3, mean) - postsamp$beta[bg,ind])^2)
	Var[bg,] <- colSums(apply(coefsLASSO[, ind, , drop=FALSE], 2:3, getPopVar))
	MSE[bg, ] <- Bias2[bg,] + Var[bg,]
	MBV[bg,] <- pmax(Bias2[bg,], Var[bg,])	
	MSEmin <- which.min(MSE[bg, ])
	MBVmin <- which.min(MBV[bg, ])
	
	estMeanMSE[bg,] <- apply(coefsLASSO[, , MSEmin], 2, mean)
	estMeanMBV[bg,] <- apply(coefsLASSO[, , MBVmin], 2, mean)
	estVarMSE[bg, , ] <- var(coefsLASSO[,, MSEmin])
	estVarMBV[bg, , ] <- var(coefsLASSO[,, MBVmin])
	}
	
	out <- list(Bias2=Bias2, Var=Var, MSE=MSE, MBV = MBV , estMeanMSE= estMeanMSE, estMeanMBV= estMeanMBV, estVarMSE= estVarMSE, estVarMBV= estVarMBV )
	out
}	

## @name festLASSO
##' @rdname festRidge
# ##' @title Compute Lasso Estimator for simulated Data
# ##
# ##' @description To be added...
# ##
# ##' @details To be added...
# ##
# ##  Input:
# ##'	@param X Design matrix of regression problem
# ##' @param y outcome vector. Typically centered.
# ##' @param loss Loss function for choosing the penaly parameter. See details.
# ##' @param penalize Vector giving penalty structure. Supplied to glmnet as `\code{penalty.factor}'.
##'	@param rescale.lambda If \code{TRUE}, then lambda is rescaled to account for the 
##'		default re-scaling done by \code{glmnet}. Can also be a scalar scaling factor.
# ##' @param scale Logical indicating whether the design matrix X be scaled. See details.
# ##'	@param returnMSE Logical indicating whether mse object should be returned.
# ##'	@param ind Vector of ntegers or logicals indicating which coefficients the loss is to be computed on.
# ##' @param lseq Sequence of penalty values to consider. Sorted in decreasing order.
##' @param B Number of future datasets to simulate for each point in posterior sample.
# ##'	@param se.version Which version of Standard errors to report?
# ##' @param postsamp List containing posterior sample (from \code{samplePosterior}). If
# ##'		missing, then a posterior sample is drawn.  Currently checks on the provided
# ##'		\code{postsamp} are limited, so use with caution.  Designed to facilitate
# ##'		simualtions or other scenarios where it may be pre-computed.
# ##' @param returnPS logical indicating whether or not the full posterior sample should
# ##'		be included in output.
# ##' @param nPost Size of posterior sample to compute
# ##' @param ... Other arguments passed to \code{samplePosterior}
# ##'
##' @export
# ##' @import glmnet
# ##
# ##' @seealso \code{\link{festRidge}}, \code{\link{simLASSO}}
# ##
festLASSO <- function(X, y, loss=c("fMSE", "fMBV", "both"), ind=1, lseq, B=500, penalize, rescale.lambda=TRUE, scale=FALSE, returnMSE=FALSE, postsamp, returnPS=FALSE, nPost=1000, se.version=c("varExp", "full", "none"), ...){
	
	fncall <- match.call()
	loss <- match.arg(loss)
	se.version <- match.arg(se.version)	
		
	if(scale) X <- scale(X)
	if(missing(lseq)) lseq <-  exp(seq(log(100), log(0.01), length=500))

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
	
	# Rescale lambda to account for re-scaling done by glmnet
	if (rescale.lambda){
		if(is.logical(rescale.lambda)){
			plam <- sum(penalize==0)
			lseq <- lseq*(p-plam)/p
		} else {
			lseq <- lseq*rescale.lambda
		}
	}
	
	# Make Posterior Draw, if missing
	if(missing(postsamp)) {
		postsamp <- samplePosterior(X, y, n= nPost, ...)
	} else {
		nPost <- length(postsamp$sigma2)
	}
	

	# Compute Bias2 and Variance, using simulated values
	lmse <- mseLASSO(lambda=lseq, X=X, postsamp=postsamp, nPost=nPost, B=B, penalize=penalize, ind=ind)

	out <- list()
	# Get estimator
	if(loss=="fMSE" | loss=="both"){
		pmse <- colMeans(lmse$MSE)
		lmin <- which.min(pmse)
		gfit <- glmnet::glmnet(x= X, y= y, lambda= lseq,  penalty.factor= penalize, standardize=FALSE, intercept=FALSE, family="gaussian")
		beta <- gfit$beta[,lmin]
		if (se.version=="varExp") {
			se <- sqrt(diag(var(lmse$estMeanMSE)))
		} else if (se.version=="full") {
			se <- sqrt(diag(var(lmse$estMeanMSE)  + apply(lmse$estVarMSE, 2:3, mean)))
		}  else if (se.version=="none") {
			se <- rep(0, length(beta))
		}
		if (scale) {
			beta <- beta/attributes(X)$"scaled:scale"
			se <- se/attributes(X)$"scaled:scale"[1]
		}
		out$fMSE <- list(beta=beta, lambda=lseq[lmin], lmin=lmin, se=se)
	}
	if (loss=="fMBV" | loss=="both"){
		pmaxBias2Var <- colMeans(lmse$MBV)
		lmin <- which.min(pmaxBias2Var)
		gfit <- glmnet::glmnet(x= X, y= y, lambda= lseq,  penalty.factor= penalize, standardize=FALSE, intercept=FALSE, family="gaussian")
		beta <- gfit$beta[,lmin]
if (se.version=="varExp") {
			se <- sqrt(diag(var(lmse$estMeanMBV)))
		} else if (se.version=="full") {
			se <- sqrt(diag(var(lmse$estMeanMBV)  + apply(lmse$estVarMBV, 2:3, mean)))
		}  else if (se.version=="none") {
			se <- rep(0, length(beta))
		}
		out$fMBV <- list(beta=beta, lambda=lseq[lmin], lmin=lmin, se=se)
	}

	if (returnMSE){
		out$MSE <- lmse$MSE
		out$Bias2 <- lmse$Bias2
		out$Var <- lmse$Var
	}
	if (returnPS){
		out$postsamp <- postsamp
	}
	return(out)
}


