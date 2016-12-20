##' @name getMSE
##' @title Mean Squared Error of a Vector
##
##' @description Compute the Mean Squared Error of a vector
##
##' @details Details here!
##
##' @param x vector of estimates
##' @param truth vector of true parameters
getMSE <- function(x, truth){
	if(any(is.na(x))) warning("Missing values in x.")
	mean((x - truth)^2, na.rm=TRUE)
}

getPopVar <- function(x) {
	if(any(is.na(x))) warning("Missing values in x.")
	mean((x - mean(x, na.rm=TRUE))^2, na.rm=TRUE)
}

getBias2 <- function(x, truth) {
	if(any(is.na(x))) warning("Missing values in x.")
	(mean(x, na.rm=TRUE) - truth)^2
}





