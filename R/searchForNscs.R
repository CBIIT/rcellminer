#' Search for NSCs 
#' 
#' @param pattern a search pattern. This string will be treated as a regular expression 
#'   with the case ignored.  
#' @return A vector of matching NSCs
#' 
#' @details Use this function with caution. Not all compounds have names and 
#'   compounds can have many synonyms not included in CellMiner. 
#'   
#' @examples 
#' searchForNscs("nib$")  
#' 
#' @concept rcellminer
#' @export
searchForNscs <- function(pattern) {
	df <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	idx <- grep(pattern, df[, "NAME"], ignore.case = TRUE)
	results <- df[idx, "NSC"]
	names(results) <- df[idx, "NAME"]
	return(results)
}
