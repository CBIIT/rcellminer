#' Compare an input pattern against a set of patterns, excluding the effect
#' of a fixed pattern.
#' 
#' @param x An N element input pattern specified as either a vector or a
#' 1 x N matrix or data frame. 
#' @param Y A k x N matrix with k patterns for comparison with the input pattern x
#' specified along the rows, with rownames set appropriately.
#' @param z An N element  pattern specified as either a vector or a
#' 1 x N matrix or data frame. This is the pattern whose effect (with respect to
#' a linear model) is to be excluded when comparing x with each row entry of Y.
#' 
#' @return A data frame with pattern comparison results (ordered by PARCOR):
#' NAME: Name of entry in Y being compared.
#' PARCOR: Partial correlation between x and the entry in Y with respect to z.
#' PVAL: p-value.
#'  
#' @examples 
#' x <- exprs(getAct(rcellminerData::drugData))["609699", ]
#' Y <- rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]][1:100, ]
#' z <- rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]]["SLFN11", ]
#' results <- parCorPatternComparison(x, Y, z)
#' 
#' @concept rcellminer
#' @export
parCorPatternComparison <- function(x, Y, z){
	if (!is.vector(x)){
		x <- as.numeric(x)
	}
	if (!is.vector(z)){
		z <- as.numeric(z)
	}
	if (is.null(rownames(Y))){
		rownames(Y) <- 1:nrow(Y)
	}
	if (any(duplicated(rownames(Y)))){
		stop("Remove duplicate rownames in input matrix Y.")
	}
	
	results <- data.frame(NAME = rownames(Y), PARCOR = NA, PVAL = NA,
												stringsAsFactors = FALSE)
	rownames(results) <- results$NAME
	
	xzNotNa <- (!is.na(x)) & (!is.na(z))
	for (name in rownames(results)){
		notNa <- xzNotNa & (!is.na(Y[name, , drop=TRUE]))

		tmp <- cor.test(
			residuals(lm(x[notNa] ~ z[notNa])), 
			residuals(lm(Y[name, notNa, drop=TRUE] ~ z[notNa])))
		
		results[name, "PARCOR"] <- tmp$estimate
		results[name, "PVAL"]   <- tmp$p.value
	}
	
	results <- results[order(results$PARCOR, decreasing = TRUE), ]
	
	return(results)
}