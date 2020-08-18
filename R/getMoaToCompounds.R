#' Get a named list mapping MOA classes to associated compound sets.
#' 
#' @return a named list mapping MOA classes to associated compound sets 
#' (each represented as a character vector).
#' 
#' @examples 
#' moaToCompounds <- getMoaToCompounds()
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom stringr str_split
#' @importFrom Biobase featureData
getMoaToCompounds <- function(){
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	knownMoaDrugAnnot <- drugAnnot[(drugAnnot$MOA != ""), ]
	
	nscToMoaClasses <- str_split(knownMoaDrugAnnot$MOA, pattern = "[|]")
	names(nscToMoaClasses) <- knownMoaDrugAnnot$NSC
	
	moaToCompounds <- list()
	for (nscStr in names(nscToMoaClasses)){
		for (moaClass in nscToMoaClasses[[nscStr]]){
			moaToCompounds[[moaClass]] <- c(moaToCompounds[[moaClass]], nscStr)
		}
	}
	
	return(moaToCompounds)
}