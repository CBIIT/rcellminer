#' Returns a vector of median sensitive cell line activity (-logGI50) values
#' for a set of compounds.
#' 
#' @param idSet a character vector specifying identifier(s) for compound(s) of interest.
#' @param senLineActZThreshold the minimum activity z-score for a sensitive cell line (default=0.5).
#' @param dataSource character string indicating data source (default="NCI60"). 
#' Currently only "NCI60" is supported.
#' @param onlyCellMinerExps a logical value indicating whether to base results strictly on  
#' experimental data included in CellMiner (default=TRUE).
#' @return a numeric vector of median sensitive cell line activity (-logGI50) values indexed 
#' by the identifiers in idSet.
#' 
#' @examples 
#' idSet <- c("609699", "740")
#' getMedSenLineActivity(idSet)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom Biobase exprs
getMedSenLineActivity <- function(idSet, senLineActZThreshold=0.5, onlyCellMinerExps=TRUE, dataSource="NCI60"){
  supportedDataSources <- c("NCI60")
  if (!is.element(dataSource, supportedDataSources)){
    stop(paste(dataSource, " is not among supported data sources.", sep = ""))
  }
  drugAct <- exprs(getAct(rcellminerData::drugData))
  
  idSet <- as.character(idSet)
  medActVals <- numeric(length(idSet))
  names(medActVals) <- idSet
  
  for (idStr in idSet){
    senLines <- colnames(drugAct[idStr, which(drugAct[idStr, ] >= senLineActZThreshold), drop=FALSE])
    if (length(senLines) == 0){
      warning(paste("No sensitive cell lines available for compound ", idStr, ".", sep = ""))
      medActVals[idStr] <- NA_real_
      next
    }
    
    actExpData <- getDrugActivityRepeatData(idStr, onlyCellMinerExps = onlyCellMinerExps)
    
    if (is.null(actExpData)){
      warning(paste("No -logGI50 activity available for compound ", idStr, ".", sep = ""))
      medActVals[idStr] <- NA_real_
      next
    }
      
    if (nrow(actExpData) > 1){
      senLineNeglogGI50 <- colMeans(actExpData[, senLines, drop=FALSE], na.rm = TRUE)
    } else{
      senLineNeglogGI50 <- as.numeric(actExpData[, senLines])
    }
    medActVals[idStr] <- median(senLineNeglogGI50, na.rm = TRUE)
  }
  
  return(medActVals)
}
