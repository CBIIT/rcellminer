#' Plot NCI-60 drug activity profiles for repeat experiments.
#' 
#' @param nscStr a string specifying the NSC identifier for a compound.
#' @param useZScore a boolean specifying whether to plot z-transformed data 
#'   (as opposed to -logGI50 values).
#' @param maxRepNum an integer specifying the maximum number of repeat experiments to plot.
#' @param pdfFilename name of a PDF output 
#' @param pdfWidth with of the PDF (default: 12)
#' @param pdfHeight with of the PDF (default: 6)
#' 
#' @return NONE
#' 
#' @examples 
#' plotDrugActivityRepeats("609699")
#' plotDrugActivityRepeats("609699", useZScore=TRUE, maxRepNum=3)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom Biobase exprs featureData
plotDrugActivityRepeats <- function(nscStr, useZScore=FALSE, maxRepNum=5, 
                                      pdfFilename=NULL, pdfWidth=12, pdfHeight=6) {
	molDataMats <- getMolDataMatrices()
	drugRepeatAct <- exprs(getRepeatAct(rcellminerData::drugData))
	drugRepeatAnnot <- as(featureData(getRepeatAct(rcellminerData::drugData)), "data.frame")
	
  #load(file.path(.pri, "RData", "lmpdb.Rdata"))
  iNsc <- which(drugRepeatAnnot$nsc == nscStr)
  iNsc <- iNsc[ drugRepeatAnnot$used_in_zscore[iNsc] ]
  stopifnot(length(iNsc) > 0)
  
  repActData <- drugRepeatAct[iNsc, , drop=FALSE]
  rownames(repActData) <- paste(nscStr, "_", (1:nrow(repActData)), sep = "")
  rowsToPlot <- rownames(repActData)[1:min(nrow(repActData), maxRepNum)]
  
  if (useZScore){
    repActData <- t(scale(t(repActData)))
    loXLimit <- -3
    hiXLimit <- 3
    xLab <- "Z-Score"
  } else{
    loXLimit <- floor(min(as.vector(repActData), na.rm = TRUE)) - 0.5
    hiXLimit <- ceiling(max(as.vector(repActData), na.rm = TRUE)) + 0.5
    xLab <- "-logGI50"
  }
  
  if(!is.null(pdfFilename)) {
      pdf(file=pdfFilename, width=pdfWidth, height=pdfHeight)      
  }
  
  plotCellMiner(drugAct = repActData, molData = molDataMats, xLimits = c(loXLimit, hiXLimit), 
                xLabel = xLab, plots=rep("drug", length(rowsToPlot)), nsc = rowsToPlot, 
                verbose = FALSE)
  
  if(!is.null(pdfFilename)) {
      dev.off()
  }    
}
