#' Returns a table of activity range statistics for a set of compounds.
#' 
#' @param nscSet a character vector specifying NSC identifier(s) for compound(s) of interest.
#' @param concFormat a string selected from "NegLogGI50M" or "IC50MicroM".
#'   "NegLogGI50M" specifies activities as the negative log of the 50% growth
#'   inhibitory concentration (molar). "IC50MicroM" specifies activities as the 50% growth 
#'   inhibitory concentration (micromolar). 
#' @param onlyCellMinerExps a logical value indicating whether to only return  
#' experimental data included in CellMiner (default=TRUE).
#' @return a table of activity range statistics for a set of compounds.
#' 
#' @examples 
#' nscSet <- c("609699", "740")
#' getActivityRangeStats(nscSet)
#' getActivityRangeStats(nscSet, concFormat="IC50MicroM")
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom stats median
getActivityRangeStats <- function(nscSet, concFormat='NegLogGI50M', onlyCellMinerExps=TRUE){
  if (!is.element(concFormat, c("NegLogGI50M", "IC50MicroM"))){
    stop(paste(concFormat, "is not supported. Select from 'NegLogGI50M' or 'IC50MicroM'."))
  }
  
  nscSet <- as.character(unique(nscSet))
  actRangeStats <- data.frame(NSC=nscSet, stringsAsFactors = FALSE)
  actRangeStats$MinActivity <- numeric(nrow(actRangeStats))
  actRangeStats$MedActivity <- numeric(nrow(actRangeStats))
  actRangeStats$MedSenLineActivity <- getMedSenLineActivity(idSet = nscSet, 
  																														onlyCellMinerExps = onlyCellMinerExps)
  actRangeStats$MaxActivity <- numeric(nrow(actRangeStats))
  rownames(actRangeStats) <- actRangeStats$NSC
  
  for (nscStr in nscSet){
    actExpData <- getDrugActivityRepeatData(nscStr, concFormat = concFormat, 
    																				onlyCellMinerExps = onlyCellMinerExps)
    if (nrow(actExpData) > 1){
      meanAct <- colMeans(actExpData, na.rm = TRUE)
    } else{
      meanAct <- as.numeric(actExpData)
    }
    
    actRangeStats[nscStr, "MinActivity"] <- min(meanAct, na.rm = TRUE)
    actRangeStats[nscStr, "MedActivity"] <- median(meanAct, na.rm = TRUE)
    actRangeStats[nscStr, "MaxActivity"] <- max(meanAct, na.rm = TRUE)
  }
  
  if (concFormat == 'IC50MicroM'){
    # Convert from -logGI50(molar) to IC50 (micromolar).
    actRangeStats$MedSenLineActivity <- (10^(-actRangeStats$MedSenLineActivity))*(10^6)
    colnames(actRangeStats) <- c("NSC", "MinIC50MicroM", "MedIC50MicroM", "MedSenLineIC50MicroM", "MaxIC50MicroM")
    actRangeStats <- actRangeStats[, c("NSC", "MinIC50MicroM", "MedSenLineIC50MicroM", "MedIC50MicroM", "MaxIC50MicroM")]
  }
  
  return(actRangeStats)
}