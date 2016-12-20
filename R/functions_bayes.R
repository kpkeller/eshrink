

##' @name samplePosterior
##' @title Posterior Sample for Bayesian Linear Regression
##
##' @description Computes and provides a sample from the posterior
##'		distribution of parameters from a Bayesian Linear regression
##'		model. 
##
##' @details Assumes a normal-inverse-gamma model, where the distribution
##'		of the regression coefficients is conditional on \eqn{sigma^2}.
##
##  Input:
##' @param X A data matrix of size \code{n} by \code{p}.
##' @param y Outcome variable
##'	@param a0,b0 Hyperparameters for inverse gamma
##'		distribution of the error varaince.
##'	@param v0inv Prior precision for the error term.  Either a single value
##'		to be repeated in a diagonal precision matrix, or a \code{p} by \code{p}
##'		matrix.
##' @param mu0 Prior mean.  Defaults to zero vector.
##' @param n Size of posterior sample
##' @param hyperPost Logical indicating whether the posterior hyperparameters are returned.
##' @param intercept Logical indicating whether an intercept is included in the model.
##'		If \code{FALSE}, then \code{y} is centered prior to computations.
##
## Output:
##'	@return  A list containing the following elements:
##' \item{sigma2}{A vector containing the posterior sample of \eqn{sigma^2} values.}
##' \item{beta}{Matrix containing the posterior sample of \eqn{beta} values.}
##' \item{postMu}{Vector containing the posterior mean (if \code{moments=TRUE}).}
##' \item{postV}{Matrix giving the posterior mean (if \code{moments=TRUE}).}
##' \item{wSS}{Mean within-cluster sum-of-squares}
##
##' @export
##' @importFrom MASS mvrnorm
##' @importFrom stats rgamma
##' @author Joshua Keller
##
samplePosterior <- function(X, y, n, a0=1, b0=5e-5, v0inv=1/1000, mu0=0, hyperPost=FALSE, intercept=FALSE){
	
	# Sample size
	nx <- nrow(X)
	# Number of parameters
	p <- ncol(X)
	
	if (!intercept){
		# Center the outcome
		ycen <- y - mean(y)	
	}

	
	if (is.matrix(v0inv)){
		V0inv <- v0inv
	} else {
		V0inv <- diag(v0inv, nrow=p)  # Prior Precision matrix for coefs		
	}

	if(length(mu0)!=p){
		mu0 <- rep(mu0[1], p)
	}


	postV <- solve(crossprod(X) + V0inv) # Posterior Variance for coefs
	postMu <- postV %*% (V0inv %*% mu0 + crossprod(X, ycen))
	an <- a0 + nx/2 # Posterior shape for (sigma2)^-1
	bn <- b0 + 1/2*(crossprod(ycen) + crossprod(mu0, V0inv %*% mu0) - crossprod(postMu, crossprod(X, ycen)))  # Posterior rate for (sigma2)^-1
	# Make the posterior draw:
	postSigma <- 1/stats::rgamma(n=n, shape=an, rate=bn)
	postBeta <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma=postV)
	postBeta <- t(t(postBeta*sqrt(postSigma)) + as.vector(postMu))
	colnames(postBeta) <-  c("X", paste0("Z", 1:(p-1)))
	if (hyperPost){
		out <- list(sigma2=postSigma, beta=postBeta, postMu=postMu, postV=postV, an=an, bn=bn)
	} else {
		out <- list(sigma2=postSigma, beta=postBeta)
	}
	return(out)
}

