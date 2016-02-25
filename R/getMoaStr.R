#' Get MOA string 
#' 
#' @param nscStr an NSC string 
#' 
#' @return a comma-delimited string with MOA
#' 
#' @details LINK TO MOAs? 
#' 
#' @examples 
#' getMoaStr("94600")
#' getMoaStr(c("94600", "609699"))
#' 
#' @concept rcellminer
#' @export
getMoaStr <- function(nscStr){
	results <- sapply(nscStr, function(nsc) {
		drugMoaList <- getDrugMoaList(nsc)
		
		if (is.null(drugMoaList)){
			moaStr <- NA_character_
		} else{
			moaStr <- paste(drugMoaList, collapse = ",")
		}
		
		return(moaStr)
	})
	
	return(results)
}