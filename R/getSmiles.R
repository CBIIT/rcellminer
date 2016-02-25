#' Get the SMILES strings for a set of NSC identifiers.
#' 
#' @param nscSet A character vector of NSC strings 
#' 
#' @return A named character vector indicating the SMILES string for
#' each NSC in nscSet (or NA if no structural information is available).
#' 
#' @examples 
#' nscSet <- c("609699", "94600")
#' getSmiles(nscSet)
#' 
#' @concept rcellminer
#' @export
getSmiles <- function(nscSet){
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	if (!identical(rownames(drugAnnot), drugAnnot$NSC)){
		rownames(drugAnnot) <- as.character(drugAnnot$NSC)
	}
	
	smilesSet <- vapply(nscSet, 
											function(nsc){
												ifelse((nsc %in% rownames(drugAnnot)),
															 drugAnnot[nsc, "SMILES"],
															 NA_character_)
											},
											character(1))
	
	return(smilesSet)
}
