#' Check if an NSC ID is public 
#' 
#' @param nscs a vector of NSC string IDs
#' @return a vector of boolean values of whether each NSC is public 
#' 
#' @examples 
#' isPublic("-1")
#' isPublic(c("-1", "609699"))
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom Biobase exprs
isPublic <- function(nscs) { 
	pubNscSet <- rownames(exprs(getAct(rcellminerData::drugData)))
  return(nscs %in% pubNscSet)
}

