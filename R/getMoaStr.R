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
#' 
#' @concept rcellminer
#' @export
getMoaStr <- function(nscStr){
  drugMoaList <- getDrugMoaList(nscStr)
  if (is.null(drugMoaList)){
  	moaStr <- NA_character_
  } else{
  	moaStr <- paste(drugMoaList, collapse = ",")
  }
  
  return(moaStr)
}