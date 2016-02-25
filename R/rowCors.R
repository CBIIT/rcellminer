#' Row-wise correlations
#' 
#' Correlation between ith row of x and ith row of y for all i
#' 
#' @param X a matrix
#' @param Y a matrix
#' @return a list of two vectors: cor (correlation values) and pval (correlation p-values)
#' 
#' @examples
#' a <- matrix(runif(100), nrow=10, ncol=10)
#' b <- matrix(runif(100), nrow=10, ncol=10)
#' c <- rowCors(a, b)
#' 
#' @author Sudhir Varma, NCI-LMP
#' 
#' @concept rcellminer
#' @export
rowCors <- function(X,Y) {
	na.vals=which(is.na(X) | is.na(Y))
	X[na.vals]=NA
	Y[na.vals]=NA
	
	X=sweep(X, 1, rowMeans(X, na.rm=TRUE), "-")
	Y=sweep(Y, 1, rowMeans(Y, na.rm=TRUE), "-")
	X=sweep(X, 1, sqrt(rowSums(X*X, na.rm=TRUE)), "/")
	Y=sweep(Y, 1, sqrt(rowSums(Y*Y, na.rm=TRUE)), "/")
	z=X*Y
	r=rowSums(z, na.rm=TRUE)
	qw=which(abs(r)>1)
	r[qw]=sign(r[qw])
	n=rowSums(!is.na(z))
	df=n-2
	qw=which(df>=0)
	t=array(data=NA, dim=length(df))
	t[qw]=sqrt(df[qw]) * r[qw] / sqrt(1 - r[qw]^2)
	pval=2*pt(-abs(t), df)
	
	return(list(cor=r, pval=pval))
}