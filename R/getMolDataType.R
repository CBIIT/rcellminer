#' Get the molecular data type prefixes for a set of features.
#' 
#' @param features A vector of features.
#' @param prefixLen The length of the molecular data type prefix.
#' @return A character vector of molecular data type prefixes.
#' 
#' #' @examples 
#' getMolDataType(c("expTP53", "copMDM2", "mutCHEK2", "mutBRAF"))
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom stringr str_sub
getMolDataType <- function(features, prefixLen = 3){
	molDataTypes <- vapply(features, function(s) str_sub(s, start = 1, end = prefixLen),
												 character(1))
	return(molDataTypes)
}