
#' Compare an input pattern against a set of patterns.
#' 
#' @param pattern An N element input pattern specified as either a named vector or an
#' 1 x N matrix or data frame.  Names (or column names) must match the column names of each element
#' of profileMatrixList.
#' @param profileMatrixList A single matrix (or data frame) or a list of matrices (or data frames).
#' Each matrix (data frame) must be k x N - that is the k patterns for comparison with the input
#' pattern must be specified along the rows, with rownames set appropriately.
#' @param method a string specifying the type of correlation, chosen from pearson 
#' (default) or spearman.
#' 
#' @return A data frame with pattern comparison results. Specifically, if M is the total number
#' patterns in profileMatrixList elements, an M x 2 matrix is returned with sorted Pearson's 
#' correlations in the first column and corresponding p-values in the second column. Comparison
#' pattern names are indicated in the row names.
#'  
#' 
#' @examples 
#' drugAct <- exprs(getAct(rcellminerData::drugData))
#' molDataMats <- getMolDataMatrices()[c("exp", "mut")]
#' molDataMats <- lapply(molDataMats, function(X) X[1:10, ])
#' pcResults <- patternComparison(drugAct["609699", ], molDataMats)
#' pcResults <- patternComparison(drugAct["609699", ], molDataMats, method="spearman")
#' pcResults <- patternComparison(drugAct["609699", ], molDataMats$exp, method="spearman")
#' 
#' @concept rcellminer
#' @export
patternComparison <- function(pattern, profileMatrixList, method = "pearson"){
	
	# ------[INPUT CHECKS AND POSSIBLE DATA TYPE TRANSFORMATIONS]------------------------------
	# (1) If input pattern is a vector, convert to a 1 x length(pattern) data.frame,
	#     so that cross-correlation computation below works properly.
	if (is.vector(pattern)){
		patternNames <- names(pattern)
		pattern <- as.data.frame(matrix(pattern, nrow=1, ncol=length(pattern)))
		rownames(pattern) <- "pattern"
		colnames(pattern) <- patternNames
	}
	
	# (2) If input profileMatrixList is actually a single matrix (or data.frame),
	#     place this in a single-element list, so that the code below (which expects
	#     a list of matrices) works properly.
	if (is.matrix(profileMatrixList) || is.data.frame(profileMatrixList)){
		tmp <- profileMatrixList
		profileMatrixList <- list()
		profileMatrixList[["1"]] <- tmp
	}
	
	# (3) Make sure that the (e.g. cell line) names in the pattern match the names associated
	#     with profiles organized along the rows of the matrices in profileMatrixList.
	for (i in (1:length(profileMatrixList))){
		if (!identical(colnames(pattern), colnames(profileMatrixList[[i]]))){
			stop(paste("names(pattern) is not identical to colnames(profileMatrixList[[", i, "]]."))
		}
	}
	# -----------------------------------------------------------------------------------------------
	
	profileNames <- NULL
	for (i in (1:length(profileMatrixList))){
		profileNames <- c(profileNames, rownames(profileMatrixList[[i]]))
	}
	
	output <- as.data.frame(matrix(0, nrow=length(profileNames), ncol=2))
	rownames(output) <- profileNames
	colnames(output) <- c("COR", "PVAL")
	
	for (i in (1:length(profileMatrixList))){
		xcorOutput <- crossCors(X=pattern, Y=profileMatrixList[[i]], method)
		output[colnames(xcorOutput$cor),  "COR"]  <- as.vector(xcorOutput$cor)
		output[colnames(xcorOutput$pval), "PVAL"] <- as.vector(xcorOutput$pval)
	}
	
	output <- output[order(output$COR, decreasing=TRUE, na.last=TRUE), ]
	
	return(output)
}