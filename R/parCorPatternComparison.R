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
#' Note that for the partial correlation to be value, the pattern(s) in Z should 
#' not overlap with those in x or Y.
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
		rownames(Y) <- paste0("Y_", 1:nrow(Y))
	}
	if (is.null(rownames(Z))){
		rownames(Z) <- paste0("Z_", 1:nrow(Z))
	}
	
	if (any(duplicated(rownames(Y)))){
		stop("Remove duplicate rownames in input matrix Y.")
	}
	if (any(duplicated(rownames(Z)))){
		stop("Remove duplicate rownames in input matrix Z.")
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
	
	numProcessed <- 0
	N <- nrow(Y)
	for (name in rownames(results)){
		lmData[, "var_Y_n"] <- Y[name, , drop = TRUE]
		notNa <- xzNotNa & (!is.na(Y[name, , drop=TRUE]))
		
		#----[check]-----------------------------------------------------------------
		# Note: For the partial correlation to be valid, The variables being compared, 
		# x and the current row of Y, should not be identical to the any of the Z 
		# variables. We are removing the effect of the Z variables on the x and Y 
		# variables by fitting linear models for these with respect to the Z variables. 
		# The residual vectors from the above model fits are then correlated to obtain 
		# the partial correlation. 
		# If some of the Z variables are identical to either x or 
		# the current Y, in theory, a zero residual vector should be obtained,
		# making the above correlation computation invalid. In pratice, numerical
		# imprecision will produce a residual vector with small, but non-zero,
		# entries. If this were to be used in the mentioned correlation, a spurious
		# partial correlation value would result.
		xzMatch <- vapply(lmData[, c(-1, -2), drop = FALSE], function(z) {
			identical(lmData[, "var_x"], z)
		}, logical(1))
		yzMatch <- vapply(lmData[, c(-1, -2), drop = FALSE], function(z) {
			identical(lmData[, "var_Y_n"], z)
		}, logical(1))
		if (any(xzMatch) || any(yzMatch)){
			next
		}
		#----------------------------------------------------------------------------

		tmp <- cor.test(
			residuals(lm(formula = var_x ~ .,   data = lmData[notNa, -2])), 
			residuals(lm(formula = var_Y_n ~ ., data = lmData[notNa, -1])))
		
		results[name, "PARCOR"] <- tmp$estimate
		results[name, "PVAL"]   <- tmp$p.value
		
		lmData[, "var_Y_n"] <- NA
		if (is.function(updateProgress)){
			numProcessed <- numProcessed + 1
			percentCompleted <- floor((numProcessed / N) * 100)
			text <- paste0(percentCompleted, "% completed.")
			updateProgress(text)
		}
	}
	
	results <- na.exclude(results)
	results <- results[order(results$PARCOR, decreasing = TRUE), ]
	
	return(results)
}