

##' @name samplePosterior
##' @title Posterior Sample for Bayesian Linear Regression
##
##' @description Draws a sample from the posterior
##'		distribution of parameters from a Bayesian Linear regression
##'		model. 
##
##' @details This function draws a sample from the posterior distributions of the coefficient parameter (\eqn{\beta}) and error variance parameter (\eqn{\sigma^2}) from a Bayesian linear regression model. The variance parameter is assumed to follow an inverse-gamma distribution. Conditional on the error variance parameter and a specified precision matrix, the coefficient parameters (\eqn{\beta}) are multivariate normal. 
##
##  Input:
##' @param X Design matrix of size \code{n} by \code{p}.
##' @param y Outcome variable
##'	@param a0,b0 Hyperparameters (shape, rate) for inverse gamma
##'		distribution of the error variance.
##'	@param v0inv Prior precision for the error term.  Either a single value
##'		to be repeated in a diagonal precision matrix, or a \code{p} by \code{p}
##'		matrix.
##' @param mu0 Prior mean. Either a single value that will be repeated,
##'		or a vector of length \code{p}.   Defaults to zero vector.
##' @param n Size of posterior sample to be computed. A value of 0 is accepted.
##' @param returnParams Logical indicating whether the parameters of the posterior distribution are returned.
##' @param intercept Logical indicating whether an intercept is included in the model.
##'		If \code{FALSE}, then \code{y} is centered.
##
## Output:
##'	@return  A list containing the following elements:
##' \item{sigma2}{A vector containing the posterior sample of \eqn{\sigma^2} values.}
##' \item{beta}{Matrix containing the posterior sample of \eqn{\beta} values.}
##' \item{postMu}{Vector containing the posterior mean (if \code{returnParams =TRUE}).}
##' \item{postV}{Matrix giving the posterior mean (if \code{returnParams =TRUE}).}
##' \item{an,bn}{Posterior hyperparameters for the  inverse gamma
##'		distribution of the error variance (if \code{returnParams =TRUE}).}
##
##' @export
##' @importFrom MASS mvrnorm
##' @importFrom stats rgamma
##' @author Joshua Keller
##
##' @examples
##' x <- rnorm(40, mean=2, sd=2)
##' y <- x + rnorm(40, sd=1)
##' samplePosterior(X=x, y=y, n=10)
##' samplePosterior(X=cbind(1, x), y=y, n=10, intercept=TRUE)
##'samplePosterior(X=cbind(1, x), y=y, n=0, mu=c(3, 3), intercept=TRUE)
##
##
samplePosterior <- function(X, y, n, a0=1, b0=5e-5, v0inv=1/1000, mu0=0, returnParams=TRUE, intercept=FALSE){
	
	X <- as.matrix(X)

	# Sample size
	nx <- nrow(X)
	# Number of parameters
	p <- ncol(X)
	
	
	if (!intercept){
		# Center the outcome
		ycen <- y - mean(y)	
	} else {
		ycen <- y
	}

	if (is.matrix(v0inv)){
		V0inv <- v0inv
	} else {
		V0inv <- diag(v0inv, nrow=p)  # Prior Precision matrix for coefs		
	}

	if(length(mu0)==1){
		mu0 <- rep(mu0, p)
	}
	if(length(mu0)!=p){
		warning("Prior mean (mu0) length not equal to 'p'. Replacing with rep(mu0[1], p).")
		mu0 <- rep(mu0[1], p)
	}
	
	if (n<0){
		stop("'n' must be greater than or equal to zero.")
	}

	# Posterior Variance for coefficients
	postV <- solve(crossprod(X) + V0inv) 
	postMu <- postV %*% (V0inv %*% mu0 + crossprod(X, ycen))
	# Posterior shape for (sigma2)^-1:
	an <- a0 + nx/2 
	# Posterior rate for (sigma2)^-1:
	bn <- b0 + 1/2*(crossprod(ycen) + crossprod(mu0, V0inv %*% mu0) - crossprod(postMu, crossprod(X, ycen))) 
	bn <- as.numeric(bn) 
	if (n==0){
		postSigma <- NULL
		postBeta <- NULL
	} else {
		# Make the posterior draw:
		postSigma <- 1/stats::rgamma(n=n, shape=an, rate=bn)
		postBeta <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma=postV)
		postBeta <- t(t(postBeta*sqrt(postSigma)) + as.vector(postMu))
		colnames(postBeta) <-  colnames(X)
	}
	out <- list(sigma2=postSigma, beta=postBeta)
	if (returnParams) out <- c(out, list(postMu=postMu, postV=postV, an=an, bn=bn))
	return(out)
}

