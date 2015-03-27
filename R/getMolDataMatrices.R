#' Returns a list of molecular data type matrices, with rownames in 
#' each matrix prefixed with a data type abbreviation.
#' 
#' @param molDataMats A named list of molecular data type matrices with feature data 
#' specified along the rows, and feature names indicated in the row names.
#' 
#' @return a list containing molecular data type matrices, with rownames 
#' in each matrix prefixed with a data type abbreviation, e.g., 'exp' for
#' mRNA expression, etc.  The matrix-specific data type abbreviations are
#' derived from the names of molDataMats.
#' 
#' 
#' @examples
#' molDataMats <- getMolDataMatrices()
#' 
#' @concept rcellminer
#' @export
getMolDataMatrices <- function(molDataMats = NULL){
	if (is.null(molDataMats)){
		molDataMats <- getAllFeatureData(rcellminerData::molData)
	}
	
	for (molDataType in names(molDataMats)){
		rownames(molDataMats[[molDataType]]) <- paste0(molDataType, 
			rownames(molDataMats[[molDataType]]))
	}
	
	return(molDataMats)
}