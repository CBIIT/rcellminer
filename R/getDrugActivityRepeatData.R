#' Returns a matrix containing repeat activity experiment data for a compound.
#' 
#' @param nscStr a string specifying the NSC identifier for the compound.
#' @param concFormat a string selected from "NegLogGI50M" or "IC50MicroM".
#' "NegLogGI50M" specifies activities as the negative log of the 50% growth
#' inhibitory concentration (molar). "IC50MicroM" specifies activities as the 50% 
#' growth inhibitory concentration (micromolar). 
#' @param onlyCellMinerExps a logical value indicating whether to only return  
#' experimental data included in CellMiner (default=TRUE).
#' @return a matrix with activity data from each experiment associated with 
#' a compound organized along the rows.
#' 
#' @examples
#' nscStr <- "609699"
#' actData <- getDrugActivityRepeatData(nscStr, concFormat='NegLogGI50M')
#' actData <- getDrugActivityRepeatData(nscStr, concFormat='IC50MicroM')
#'
#' @concept rcellminer
#' @export
#' 
#' @importFrom Biobase exprs featureData
getDrugActivityRepeatData <- function(nscStr, concFormat = 'NegLogGI50M', onlyCellMinerExps=TRUE){
  if (!is.element(concFormat, c("NegLogGI50M", "IC50MicroM"))){
    stop(paste(concFormat, "is not supported. Select from 'NegLogGI50M' or 'IC50MicroM'."))
  }
  drugRepeatAct <- exprs(getRepeatAct(rcellminerData::drugData))
  drugRepeatAnnot <- as(featureData(getRepeatAct(rcellminerData::drugData)), "data.frame")
  
  iNsc <- which(drugRepeatAnnot$nsc == nscStr)
  if (length(iNsc) == 0){
    warning(paste("No -logGI50 activity available for compound ", nscStr, ".", sep = ""))
    return(NULL)
  }
  
  if (onlyCellMinerExps){
    iNsc <- iNsc[ drugRepeatAnnot$used_in_zscore[iNsc] ]
  }
  
  actExpData <- drugRepeatAct[iNsc, , drop=FALSE]
  
  if (concFormat == 'IC50MicroM'){
    # actExpData is -log10GI50 (molar)
    actExpData <- (10^(-actExpData))*(10^6)
  }
  
  rownames(actExpData) <- paste(nscStr, "_", seq(nrow(actExpData)), sep = "")
  
  return(actExpData)
}