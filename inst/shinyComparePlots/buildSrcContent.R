library(rcellminer)

config <- jsonlite::fromJSON("config.json")

source("appUtils.R")
source("dataLoadingFunctions.R")

srcContent <- lapply(config, loadSourceContent)
isLoadedSrc <- vapply(srcContent, function(x) { !is.null(x) }, logical(1))
if (any(!isLoadedSrc)){
	srcContent <- srcContent[isLoadedSrc]
}

# For NCI-60, replace default color map to use CellMiner tissue type colors.
nci60ColorTab <- loadNciColorSet(returnDf=TRUE)
nci60ColorTab$OncoTree1 <- srcContent$nci60$sampleData$OncoTree1
srcContent$nci60$tissueColorMap <- c(by(nci60ColorTab, nci60ColorTab$OncoTree1, 
																				FUN = function(x) unique(x$colors)))

saveRDS(srcContent, "srcContent.rds", compress = FALSE)