#' Calculate cross-correlations with between rows of input matrices
#' 
#' @param X a matrix or data.frame 
#' @param Y a matrix or data.frame 
#' @param method a string specifying the type of correlation, chosen from pearson 
#' (default) or spearman.
#' @return a list containing matrices of pairwise correlations and their p-values 
#' between rows of the input matrices or dataframes.
#' 
#' @author Sudhir Varma, NCI-LMP, with input checks, support for Spearman's correlation
#' added by VNR.
#' 
#' @examples
#' drugActData <- exprs(getAct(rcellminerData::drugData))
#' crossCors(drugActData[c("94600"), ], drugActData[c("727625", "670655"), ])
#' crossCors(drugActData[c("94600"), ], drugActData[c("727625", "670655"), ], method="spearman")
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom stats cor.test pt
crossCors <- function(X, Y = NULL, method = "pearson") {
	if (!(method %in% c("pearson", "spearman"))){
		stop("Parameter 'method' must be set to either 'pearson' or 'spearman'.")
	}

	if (is.null(Y)){
		Y <- X
	}
	if (is.data.frame(X)){
		X <- as.matrix(X)
	}
	if (is.data.frame(Y)){
		Y <- as.matrix(Y)
	}
	if (is.vector(X)){
		X <- matrix(data=X, nrow=1, ncol=length(X))
	}
	if (is.vector(Y)){
		Y <- matrix(data=Y, nrow=1, ncol=length(Y))
	}
	
	if (method == "spearman"){
		return(crossCorsSpearman(X, Y))
	}
	
	# fix for error when nrow(Y) == 1
	if ((nrow(X) == 1) && (nrow(Y) == 1)){
		tmp <- cor.test(as.numeric(X), as.numeric(Y))
		corMat <- matrix(tmp$estimate, nrow = 1, ncol = 1)
		rownames(corMat) <- rownames(X)
		colnames(corMat) <- rownames(Y)
		pvalMat <- matrix(tmp$p.value, nrow = 1, ncol = 1)
		rownames(pvalMat) <- rownames(X)
		colnames(pvalMat) <- rownames(Y)
		return(list(cor=corMat, pval=pvalMat))
	}
	
	# fix for error when nrow(Y) == 1
	if ((nrow(X) > 1) && (nrow(Y) == 1)){
		return(lapply(crossCors(Y, X), FUN = t))
	}
	
	r=array(data=NA, dim=c(nrow(X), nrow(Y)))
	pval=array(data=NA, dim=c(nrow(X), nrow(Y)))
	for(i in 1:nrow(X))
	{
		x=X[rep(i, nrow(Y)),]
		y=Y
		na.vals=which(is.na(x) | is.na(y))
		x[na.vals]=NA
		y[na.vals]=NA
		
		x=sweep(x, 1, rowMeans(x, na.rm=TRUE), "-")
		y=sweep(y, 1, rowMeans(y, na.rm=TRUE), "-")
		x=sweep(x, 1, sqrt(rowSums(x*x, na.rm=TRUE)), "/")
		y=sweep(y, 1, sqrt(rowSums(y*y, na.rm=TRUE)), "/")
		z=x*y
		r[i,]=rowSums(z, na.rm=TRUE)
		qw=which(abs(r[i,])>1)
		r[i,qw]=sign(r[i,qw])
		n=rowSums(!is.na(z))
		df=n-2
		qw=which(df>=0)
		t=array(data=NA, dim=length(df))
		t[qw]=sqrt(df[qw]) * r[i,qw] / sqrt(1 - r[i,qw]^2)
		pval[i,]=2*pt(-abs(t), df)
		
	}
	
	if(!is.null(rownames(X)))
		rownames(r)=rownames(pval)=rownames(X)
	if(!is.null(rownames(Y)))
		colnames(r)=colnames(pval)=rownames(Y)
	
	return(list(cor=r, pval=pval))
}

#' Calculate Spearman's correlations with between rows of input matrices
#' 
#' @param X a matrix or data.frame 
#' @param Y a matrix or data.frame 
#' @return a list containing matrices of pairwise Spearman's correlations and 
#' their p-values between rows of the input matrices or dataframes.
#' 
#' 
#' @examples
#' \dontrun{
#' crossCorsSpearman(drugActData[c("94600"), ], drugActData[c("727625", "670655"), ])
#' }
#' 
#' @concept rcellminer
crossCorsSpearman <- function(X, Y = NULL){
	if (is.null(Y)){
		Y <- X
	}
	if (ncol(X) != ncol(Y)){
		stop("X and Y must have the same number of columns.")
	}
	output <- list()
	output$cor <- matrix(NA, nrow = nrow(X), ncol = nrow(Y))
	rownames(output$cor) <- rownames(X)
	colnames(output$cor) <- rownames(Y)
	output$pval <- matrix(NA, nrow = nrow(X), ncol = nrow(Y))
	rownames(output$pval) <- rownames(X)
	colnames(output$pval) <- rownames(Y)
	
	# Note: look for faster ways to compute correlations.
	for (i in seq_len(nrow(X))){
		for (j in seq_len(nrow(Y))){
			xvec <- as.numeric(X[i, ])
			yvec <- as.numeric(Y[j, ])
			corResults <- suppressWarnings(cor.test(xvec, yvec, method = "spearman", ))
			output$cor[i, j] <- corResults$estimate
			output$pval[i, j] <- corResults$p.value
		}
	}
	
	return(output)
}

