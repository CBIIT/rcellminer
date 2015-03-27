#' Returns a matrix containing activity (-logGI50) data for a set of compounds.
#'
#' @param nscSet A string specifying the NSC identifiers for the compounds.
#' @param onlyCellMinerExps A logical value indicating whether to compute results
#' using only experimental data included in CellMiner (default=TRUE).
#' @return a matrix with NCI-60 average (over experiments) -logGI50 activity
#' data; compound activity profiles are along rows.
#'
#' @examples
#' nscSet <- c("141540", "123127") # Etoposide, Doxorubicin.
#' actData <- getDrugActivityData(nscSet)
#'
#' @concept rcellminer
#' @export
#' 
#' @importFrom Biobase exprs
getDrugActivityData <- function(nscSet, onlyCellMinerExps=TRUE){
	drugAct <- exprs(getAct(rcellminerData::drugData))
	nscSet <- as.character(nscSet)
	actData <- matrix(NA, nrow=length(nscSet), ncol=60)
	rownames(actData) <- nscSet
	colnames(actData) <- colnames(drugAct)
	
	for (nsc in nscSet){
		actExpData <- getDrugActivityRepeatData(nsc, onlyCellMinerExps = onlyCellMinerExps)
		
		if (nrow(actExpData) > 1){
			avgActVec <- colMeans(actExpData, na.rm = TRUE)
			avgActVec[is.nan(avgActVec)] <- NA
		} else{
			avgActVec <- as.numeric(actExpData)
		}
		
		actData[nsc, ] <- avgActVec
	}

  return(actData)
}
