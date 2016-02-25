#' Restricts a feature matrix to only include features associated with a specified gene set.
#' 
#' @param geneSet a character vector of gene names.
#' @param featureMat a matrix or data frame with feature vectors along rows and feature
#' names specified in rownames(featureMat).
#' @param prefixSet a set of feature name prefixes to be prepended to each element of geneSet to 
#' obtain a collection of geneSet-associated features.
#' @return a matrix containing the features in the intersection of rownames(featureMat) and
#' the set of geneSet-derived features (obtained by prepending each element of prefixSet to 
#' each gene in geneSet).
#' 
#' #' @examples 
#' X <- matrix(1:25, nrow=5)
#' rownames(X) <- c("expA", "expB", "copC", "mutC", "expD")
#' restrictFeatureMat(geneSet = c("B", "C"), X)
#' 
#' @concept rcellminer
#' @export
restrictFeatureMat <- function(geneSet, featureMat, prefixSet=c("cop", "exp", "mut")){
	names(prefixSet) <- prefixSet
	geneSetFeatures <- c(lapply(prefixSet, FUN=function(s) paste(s, geneSet, sep="")), recursive=TRUE)
	
	return(featureMat[intersect(rownames(featureMat), geneSetFeatures), ])
}

