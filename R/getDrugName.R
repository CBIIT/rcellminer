#' Get the drug names for a set of NSC identifiers.
#' 
#' @param nscSet A character vector of NSC strings 
#' 
#' @return A named character vector indicating the compound names
#' for each NSC in nscSet (with an empty string returned if no
#' such information is available, and an NA returned if the
#' NSC is not included in the CellMiner database).
#' 
#' @examples 
#' nscSet <- c("609699", "94600")
#' getDrugName(nscSet)
#' 
#' @concept rcellminer
#' @export
getDrugName <- function(nscSet){
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	if (!identical(rownames(drugAnnot), drugAnnot$NSC)){
		rownames(drugAnnot) <- as.character(drugAnnot$NSC)
	}
	
	drugNames <- vapply(nscSet, 
											function(nsc){
												ifelse((nsc %in% rownames(drugAnnot)),
															 drugAnnot[nsc, "NAME"],
															 NA_character_)
											},
											character(1))
	
	return(drugNames)
}