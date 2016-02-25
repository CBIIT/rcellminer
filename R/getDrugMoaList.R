#' Get a list of applicable MOA strings for a drug.
#' 
#' @param nsc An NSC string. 
#' @param moaToCompoundListMap A named list of character vectors, with each name
#' indicating an MOA class, and its corresponding character vector specifying MOA-associated
#' drugs. If unspecified, this is constructed based on MOA information provided by CellMiner.
#' 
#' @return A character vector giving all MOA classes for the drug.
#' 
#' @details LINK TO MOAs? 
#' 
#' @examples 
#' getDrugMoaList("754365")
#' 
#' @concept rcellminer
#' @export
getDrugMoaList <- function(nsc, moaToCompoundListMap=NULL){
	if (is.null(moaToCompoundListMap)){
		moaToCompoundListMap <- getMoaToCompounds()
	}
	nsc <- as.character(nsc)
	drugMoaList <- NULL
	
	moaMatches <- vapply(moaToCompoundListMap, function(x) is.element(nsc, x), logical(1))
	if (any(moaMatches)){
		drugMoaList <- names(moaMatches[which(moaMatches)])
	}
	
	return(drugMoaList)
}
