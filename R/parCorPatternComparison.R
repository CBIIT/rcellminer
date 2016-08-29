#' Compare an input pattern against a set of patterns, excluding the predictive
#' effect of a fixed pattern or set of patterns.
#' 
#' @param x An N element input pattern specified as either a vector or a
#' 1 x N matrix or data frame. 
#' @param Y An N element  pattern specified as a vector for comparison with the input
#' pattern x or a k x N matrix with k patterns for comparison with the input pattern x
#' specified along the rows, with rownames set appropriately.
#' @param Z An N element  pattern specified as a vector or a k x N matrix of
#' patterns specified along the rows. These are the patterns whose effect (with respect 
#' to a linear model) is to be excluded when comparing x with Y or each row entry of Y.
#' @param updateProgress A optional function to be invoked with each computed
#' partial correlation to indicate progress.
#' 
#' @return A data frame with pattern comparison results (ordered by PARCOR):
#' NAME: Name of entry in Y being compared.
#' PARCOR: Partial correlation between x and the entry in Y with respect to Z.
#' PVAL: p-value.
#'  
#' @examples 
#' x <- exprs(getAct(rcellminerData::drugData))["609699", ]
#' Y <- rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]][1:100, ]
#' Z <- rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]][c("SLFN11", "JAG1"), ]
#' results <- parCorPatternComparison(x, Y, Z)
#' Y <- rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]][1, , drop=TRUE]
#' Z <- rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]]["SLFN11", , drop=TRUE]
#' results <- parCorPatternComparison(x, Y, Z)
#' 
#' @concept rcellminer
#' @export
parCorPatternComparison <- function(x, Y, Z, updateProgress = NULL){
	if (!is.vector(x)){
		x <- as.numeric(x)
	}
	if (is.vector(Y)){
		Y <- matrix(Y, nrow = 1, ncol = length(Y))
	}
	if (is.vector(Z)){
		Z <- matrix(Z, nrow = 1, ncol = length(Z))
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
	
	# Compute logical vector indicating Z matrix columns w/o missing data.
	zNotNa <- apply(Z, MARGIN = 2, FUN = function(w){
		!any(is.na(w))
	})
	xzNotNa <- (!is.na(x)) & zNotNa
	
	lmData <- data.frame(var_x = x, var_Y_n = NA)
	for (j in seq_len(nrow(Z))){
		lmData[, paste0("var_Z_", j)] <- Z[j, , drop = TRUE]
	}
	
	for (name in rownames(results)){
		lmData[, "var_Y_n"] <- Y[name, , drop = TRUE]
		notNa <- xzNotNa & (!is.na(Y[name, , drop=TRUE]))

		tmp <- cor.test(
			residuals(lm(formula = var_x ~ .,   data = lmData[notNa, -2])), 
			residuals(lm(formula = var_Y_n ~ ., data = lmData[notNa, -1])))
		
		results[name, "PARCOR"] <- tmp$estimate
		results[name, "PVAL"]   <- tmp$p.value
		
		lmData[, "var_Y_n"] <- NA
		if (is.function(updateProgress)){
			updateProgress()
		}
	}
	
	results <- results[order(results$PARCOR, decreasing = TRUE), ]
	
	return(results)
}