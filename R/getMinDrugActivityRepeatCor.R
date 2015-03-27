#' Returns a table indicating, for each compound in a specified set, the least significant 
#' correlation and associated p-value between its replicate experiments.
#' 
#' @param nscSet a character vector specifying NSC identifier(s) for compound(s) of interest.
#' @return a dataframe containing the following columns:
#' NSC, cor, pval
#' 
#' @examples 
#' nscSet <- c("123528", "339316")
#' repExpCorTab <- getMinDrugActivityRepeatCor(nscSet)
#' 
#' @concept rcellminer
#' @export
#'
#' @importFrom Biobase exprs featureData
getMinDrugActivityRepeatCor <- function(nscSet){
	drugRepeatAct <- exprs(getRepeatAct(rcellminerData::drugData))
	drugRepeatAnnot <- as(featureData(getRepeatAct(rcellminerData::drugData)), "data.frame")
    
  nscSet <- as.character(nscSet)
  repExpCorTab <- data.frame(NSC=nscSet, stringsAsFactors=FALSE)
  rownames(repExpCorTab) <- nscSet
  repExpCorTab$cor <- rep(NA, nrow(repExpCorTab))
  repExpCorTab$pval <- rep(NA, nrow(repExpCorTab))
  
  for (nscStr in nscSet){
    iNsc <- which(drugRepeatAnnot$nsc == nscStr)
    iNsc <- iNsc[ drugRepeatAnnot$used_in_zscore[iNsc] ]
    if (length(iNsc) == 0){
    	warning(paste0("Data for NSC ", nscStr, " is not available in CellMiner."))
    	next
    }
    
    repActData <- drugRepeatAct[iNsc, , drop=FALSE]
    rownames(repActData) <- paste(nscStr, "_", (1:nrow(repActData)), sep = "")
    
    if (nrow(repActData) > 1){
      xCorDat <- crossCors(repActData)
      iMax <- which(xCorDat$pval == max(xCorDat$pval, na.rm = TRUE), arr.ind=TRUE)[1, ]
      repExpCorTab[nscStr, "cor"] <- round(xCorDat$cor[iMax["row"], iMax["col"]], digits = 3)
      repExpCorTab[nscStr, "pval"] <- signif(xCorDat$pval[iMax["row"], iMax["col"]], digits = 3)
    }
  }
  
  return(repExpCorTab)
}


