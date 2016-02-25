## ----knitrSetup, include=FALSE-------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', 
							 fig.align="center", fig.width=8, fig.height=8, tidy=FALSE)

## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("rcellminer")
## biocLite("rcellminerData")

## ----loadLibrary, message=FALSE, warning=FALSE---------------------------
library(rcellminer)
library(rcellminerData)

## ----searchHelp, eval=FALSE, tidy=FALSE----------------------------------
## help.search("rcellminer")

## ----searchForNscs-------------------------------------------------------
searchForNscs("nib$")  

## ----plotCellminer-------------------------------------------------------
# Get Cellminer data
drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()

# One drug
nsc <- "94600"
plots <- c("drug") 
plotCellMiner(drugAct, molData, plots, nsc, NULL)

# One expression
gene <- "TP53"
plots <- c("exp") 
plotCellMiner(drugAct, molData, plots, NULL, gene)

# Two drugs: Irinotecan and topotecan 
nsc <- c("616348", "609699")
plots <- c("drug", "drug") 
plotCellMiner(drugAct, molData, plots, nsc, NULL)

# Two genes 
# NOTE: subscript out of bounds Errors likely mean the gene is not present for that data type
gene <- c("TP53", "MDM2")
plots <- c("exp", "mut", "exp") 
plotCellMiner(drugAct, molData, plots, NULL, gene)

# Gene and drug to plot 
nsc <- "609699"
gene <- "SLFN11"
plots <- c("exp", "cop", "drug") 
plotCellMiner(drugAct, molData, plots, nsc, gene)

## ----plotDrugSets--------------------------------------------------------
# Get CellMiner data
drugAct <- exprs(getAct(rcellminerData::drugData))

drugs <- searchForNscs("camptothecin")
drugAct <- drugAct[drugs,]
mainLabel <- paste("Drug Set: Camptothecin Derivatives, Drugs:", length(drugs), sep=" ")

plotDrugSets(drugAct, drugs, mainLabel) 

## ----plotStructures, fig.width=3, fig.height=3---------------------------
plotStructuresFromNscs("609699")

## ----makeDrugInfoTable, results='hide', message=FALSE--------------------
# getFeatureAnnot() returns a named list of data frames with annotation data 
# for drugs ("drug") and drug activity repeat experiments ("drugRepeat").
drugAnnot <- getFeatureAnnot(rcellminerData::drugData)[["drug"]]

# Get the names and MOA for compounds with MOA entries
knownMoaDrugs <- unique(c(getMoaToCompounds(), recursive = TRUE))
knownMoaDrugInfo <- data.frame(NSC=knownMoaDrugs, stringsAsFactors = FALSE)
knownMoaDrugInfo$Name <- drugAnnot[knownMoaDrugInfo$NSC, "NAME"]

# Process all NSCs 
knownMoaDrugInfo$MOA <- sapply(knownMoaDrugInfo$NSC, getMoaStr)

# Order drugs by mechanism of action.
knownMoaDrugInfo <- knownMoaDrugInfo[order(knownMoaDrugInfo$MOA), ]

# Drug_MOA_Key data frame provides further details on MOA abbrevations.
Drug_MOA_Key[c("A2", "A7"), ]

## ----computeGI50Data, results='hide', message=FALSE----------------------
negLogGI50Data <- getDrugActivityData(nscSet = knownMoaDrugInfo$NSC)

gi50Data <- 10^(-negLogGI50Data)

## ----makeIntegratedTable, results='hide', message=FALSE------------------
knownMoaDrugAct <- as.data.frame(cbind(knownMoaDrugInfo, gi50Data), stringsAsFactors = FALSE)

# This table can be written out to a file
#write.table(knownMoaDrugAct, file="knownMoaDrugAct.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, na="NA")

## ----compareFingerPrints, results='hide', message=FALSE, warning=FALSE----
# Load sqldf
library(sqldf)

# Set up necessary data
# getFeatureAnnot() returns a named list of data frames with annotation data 
# for drugs ("drug") and drug activity repeat experiments ("drugRepeat").
df <- getFeatureAnnot(rcellminerData::drugData)[["drug"]]
## Drug activities 
drugAct <- exprs(getAct(rcellminerData::drugData))
## Molecular profiling data
molData <- getMolDataMatrices() 

# Example filter on particular properties of the compounds
tmpDf <- sqldf("SELECT NSC, SMILES FROM df WHERE SMILES != ''")

# Compare against known MOA compounds demonstration
knownMoaDrugs <- unique(c(getMoaToCompounds(), recursive=TRUE))
ids <- tmpDf[tmpDf$NSC %in% knownMoaDrugs, "NSC"]
smiles <- tmpDf[tmpDf$NSC %in% knownMoaDrugs, "SMILES"]

# All public
#ids <- tmpDf$nsc
#smiles <- tmpDf$smiles

drugOfInterest <- "Camptothecin"
smilesOfInterest <- "CCC1(C2=C(COC1=O)C(=O)N3CC4=CC5=CC=CC=C5N=C4C3=C2)O"

# Make a vector of all the compounds to be pairwise compared 
ids <- c(drugOfInterest, ids)
smiles <- c(smilesOfInterest, smiles)

## ----runComparison, message=FALSE----------------------------------------
# Run fingerprint comparison 
results <- compareFingerprints(ids, smiles, verbose=FALSE)

# Display the top 5 tanimoto similarity coefficient results. The first compound 
# is a self comparison.
results[1:5]

## ----plotSimilarStructures, fig.height=2, fig.width=8--------------------
# Plot top 4 results 
resultsIdx <- sapply(names(results)[2:5], function(x) { which(tmpDf$NSC == x) })
resultsIds <- names(results)[2:5]
resultsSmiles <- tmpDf$SMILES[resultsIdx]

resultsIds <- c(drugOfInterest, resultsIds)
resultsSmiles <- c(smilesOfInterest, resultsSmiles)

# plotStructures() allows plotting user-included SMILES-based structures
plotStructures(resultsIds, resultsSmiles, titleCex=0.5, mainLabel="Fingerprint Results")

## ----plotCellMiner-------------------------------------------------------
nscs <- names(results)[2:5]
plotCellMiner(drugAct=drugAct, molData=molData, plots=rep("drug", length(nscs)), nscs, NULL)

## ----reinhold2015--------------------------------------------------------
# Gene of interest
gene <- "PTEN"

# Get Data
drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()

# Get the cell lines names for cell lines meeting particular thresholds
copKnockdown <- names(which(molData[["cop"]][paste0("cop", gene), ] < -1))
expKnockdown <- names(which(molData[["exp"]][paste0("exp", gene), ] < -1.5))
mutKnockdown <- names(which(molData[["mut"]][paste0("mut", gene), ] == 1))

# Make composite pattern
pattern <- rep(0, length(molData[["cop"]][paste0("cop", gene), ]))
names(pattern) <- names(molData[["cop"]][paste0("cop", gene), ])
tmp <- Reduce(union, list(copKnockdown, expKnockdown, mutKnockdown))
pattern[tmp] <- 1

# Composite plot data
extraPlot <- list()
extraPlot[["title"]] <- "Composite Pattern"
extraPlot[["label"]] <- "Knockdown Composite (Binary)"
extraPlot[["values"]] <- pattern

# Plot data
plotCellMiner(molData=molData, plots=c("cop", "exp", "mut"), gene=gene, extraPlot=extraPlot)

# Significant drug correlations to FDA-approved/Clinical trial drugs
# getFeatureAnnot() returns a named list of data frames with annotation data 
# for drugs ("drug") and drug activity repeat experiments ("drugRepeat").
drugAnnot <- getFeatureAnnot(rcellminerData::drugData)[["drug"]]
tmpDA <- drugAnnot[drugAnnot$FDA_STATUS != "-", c("NSC", "FDA_STATUS")]

tmpDrugAct <- drugAct[rownames(drugAct) %in% tmpDA$NSC,]

results <- patternComparison(pattern, tmpDrugAct)
sigResults <- results[results$PVAL < 0.01, ]

# Get names of drugs that are FDA approved or in clinical trials
getDrugName(rownames(sigResults))


## ----varma2014, fig.height=7.5, fig.width=6------------------------------
molData <- getMolDataMatrices()

# Get Data
copGenes <- removeMolDataType(rownames(molData[["cop"]]))
expGenes <- removeMolDataType(rownames(molData[["exp"]]))
genes <- intersect(copGenes, expGenes)

# Generate the appropriate rownames 
expGeneLabels <- paste0("exp", genes)
copGeneLabels <- paste0("cop", genes)

a <- molData[["exp"]][expGeneLabels,]
b <- molData[["cop"]][copGeneLabels,]
allGenes <- rowCors(a, b)

#selectedOncogenes <- c("FZD1", "JAK2", "ALK", "PIK3CG", "RET", "CDC25A", "PDGFB", "PIK3CB", "WNT3")
selectedOncogenes <- c("ABL1", "ALK", "BRAF", "CCND1", "CCND3", "CCNE1", "CCNE2", 
											 "CDC25A", "EGFR", "ERBB2", "EZH2", "FOS", "FOXL2", "HRAS", 
											 "IDH1", "IDH2", "JAK2", "KIT", "KRAS", "MET", "MOS", "MYC", 
											 "NRAS", "PDGFB", "PDGFRA", "PIK3CA", "PIK3CB", "PIK3CD", 
											 "PIK3CG", "PIM1", "PML", "RAF1", "REL", "RET", "SRC", "STK11", 
											 "TP63", "WNT10B", "WNT4", "WNT2B", "WNT9A", "WNT3", "WNT5A", 
											 "WNT5B", "WNT10A", "WNT11", "WNT2", "WNT1", "WNT7B", "WISP1", 
											 "WNT8B", "WNT7A", "WNT16", "WISP2", "WISP3", "FZD5", "FZD1")

# Generate the appropriate rownames 
expGeneLabels <- paste0("exp", selectedOncogenes)
copGeneLabels <- paste0("cop", selectedOncogenes)

a <- molData[["exp"]][expGeneLabels,]
b <- molData[["cop"]][copGeneLabels,]
selectedOncogenesCor <- rowCors(a, b)

hist(allGenes$cor, main="", xlab="Pearson correlation between expression and copy number", breaks=200, col="lightgray", border="lightgray")

segments(x0=median(allGenes$cor), y0=0, x1=median(allGenes$cor), y1=175, lty=2, lwd=2)
text(median(allGenes$cor)+0.02, y=175, adj=0, labels="Median Correlation\nAll Genes", cex=0.75)

segments(x0=median(selectedOncogenesCor$cor), y0=0, x1=median(selectedOncogenesCor$cor), y1=140, lty=2, lwd=2, col="red")
text(median(selectedOncogenesCor$cor)+0.02, y=140, adj=0, labels="Median Correlation\nOncogenes", cex=0.75)

rug(selectedOncogenesCor$cor, col="red") 

## ----idTopo1InhCors------------------------------------------------------
# Get normalized (Z-score) NCI-60 gene expression and drug activity data.
nci60DrugActZ <- exprs(getAct(rcellminerData::drugData))
nci60GeneExpZ <- getAllFeatureData(rcellminerData::molData)[["exp"]]

antiApoptosisGenes <- c("BAG1", "BAG3", "BAG4", "BCL10", "BCL2", 
												"BCL2A1", "BCL2L1", "BCL2L10", "BCL2L2", 
												"BFAR", "BIRC3", "BIRC6", "BNIP1", "BNIP2", 
												"BNIP3", "BNIP3L", "BRAF", "CASP3", "CD27", 
												"CD40LG", "CFLAR", "CIDEA", "DAPK1", "DFFA",
												"FAS", "IGF1R", "MCL1", "NOL3", "TMBIM1", 
												"TP53", "TP73", "XIAP")
camptothecinNsc <- "94600"

# Compute table of correlations between camptothecin activity and anti-apoptosis gene expression.
pcResults <- patternComparison(nci60DrugActZ[camptothecinNsc, ], nci60GeneExpZ[antiApoptosisGenes, ])

# Adjust p-values for multiple comparisons, sort with respect to q-values.
pcResults$QVAL <- p.adjust(pcResults$PVAL, method = "fdr")
pcResults <- pcResults[order(pcResults$QVAL), ]

# Identify genes with significant negative correlations (FDR < 0.1)
pcResults[((pcResults$COR < 0) & (pcResults$QVAL < 0.10)), ]

## ----plotCptActivityVsGeneExp--------------------------------------------
colorTab <- loadNciColorSet(returnDf = TRUE)
tissueColorTab <- unique(colorTab[, c("tissues", "colors")])

plotData <- as.data.frame(t(rbind(nci60DrugActZ[camptothecinNsc, , drop=FALSE], 
																	nci60GeneExpZ[c("TMBIM1", "BCL2L2"), ])))
colnames(plotData) <- c("Camptothecin", "TMBIM1", "BCL2L2")

plot(x=plotData$TMBIM1, y=plotData$Camptothecin, col=colorTab$colors, pch=16,
		 xlab="TMBIM1 mRNA exp (Z-score)", ylab="Camptothecin Activity (Z-score)",
		 main=paste0("Camptothecin Activity vs. TMBIM Expression (r = ", 
		 						round(pcResults["TMBIM1", "COR"], 2), ")"))
abline(lm(formula("Camptothecin ~ TMBIM1"), plotData), col="red")
legend("bottomleft", legend=tissueColorTab$tissues, col=tissueColorTab$colors, 
			 cex=0.7, pch = 16)

plot(x=plotData$BCL2L2, y=plotData$Camptothecin, col=colorTab$colors, pch=16,
		 xlab="BCL2L2 mRNA exp (Z-score)", ylab="Camptothecin Activity (Z-score)",
		 main=paste0("Camptothecin Activity vs. BCL2L2 Expression (r = ", 
		 						round(pcResults["BCL2L2", "COR"], 2), ")"))
abline(lm(formula("Camptothecin ~ BCL2L2"), plotData), col="red")
legend("bottomleft", legend=tissueColorTab$tissues, col=tissueColorTab$colors, 
			 cex=0.7, pch = 16)

## ----getCellMinerData----------------------------------------------------
# Get Cellminer data
# getFeatureAnnot() returns a named list of data frames with annotation data 
# for drugs ("drug") and drug activity repeat experiments ("drugRepeat").
drugAnnot <- getFeatureAnnot(rcellminerData::drugData)[["drug"]]
fdaDrugAnnot <- drugAnnot[which(drugAnnot$FDA_STATUS == "FDA approved"), ]

nci60FdaDrugAct <- getDrugActivityData(nscSet = fdaDrugAnnot$NSC)
nci60GeneExp <- getAllFeatureData(rcellminerData::molData)[["xai"]]
nci60GeneMut <- getAllFeatureData(rcellminerData::molData)[["mut"]]

## ----vizDnarGeneMutations, fig.width=10, fig.height=7--------------------
dnarGeneSet <- c("APC", "APLF", "ATAD5", "ATM", "CLSPN", "ERCC6", "FANCI", 
                 "FANCM", "GEN1", "HLTF", "MLH1", "POLD1", "POLE", "POLG",
                 "POLQ", "RAD54L", "REV3L", "RMI1", "SLX4", "SMARCA4", "SMC2",
                 "TP53", "WRN")

dnarGeneMut <- nci60GeneMut[dnarGeneSet, ]

# Identify most frequently mutated genes.
numMutLines <- apply(dnarGeneMut, MARGIN = 1, FUN = sum)
sort(numMutLines, decreasing = TRUE)

dnarGeneMutPlotData <- dnarGeneMut[order(numMutLines), ]
heatmap(dnarGeneMutPlotData, scale="none", Rowv = NA, Colv = NA, col = c("grey", "red"), 
				main="Deleterious mutations")

## ----vizDnarGeneExpKnockouts, fig.width=10, fig.height=7-----------------
dnarGeneExp <- nci60GeneExp[dnarGeneSet, ]

hist(dnarGeneExp, xlab="mRNA expression (log2 intensity)",
		 main="Distribution of DNA repair gene expression values")

# Set low expression threshold to 1st percentile of expression values.
lowExpThreshold <- quantile(dnarGeneExp, na.rm = TRUE, probs = c(0.01))
lowExpThreshold

# Construct binary (potential) expression knockout matrix.
dnarGeneExpKo <- matrix(0, nrow=nrow(dnarGeneExp), ncol=ncol(dnarGeneExp))
rownames(dnarGeneExpKo) <- rownames(dnarGeneExp)
colnames(dnarGeneExpKo) <- colnames(dnarGeneExp)
dnarGeneExpKo[which(dnarGeneExp < lowExpThreshold)] <- 1

# Restrict to wild type expression knockouts.
dnarGeneExpKo[which(dnarGeneMut == 1)] <- 0

# Identify genes with most frequent loss of expression.
numExpKoLines <- apply(dnarGeneExpKo, MARGIN = 1, FUN = sum)
sort(numExpKoLines, decreasing = TRUE)

dnarGeneExpKoPlotData <- dnarGeneExpKo[order(numExpKoLines), ]
heatmap(dnarGeneExpKoPlotData, scale="none", Rowv = NA, Colv = NA, 
				col = c("grey", "blue"), main="Potential expression knockouts")

## ----vizDnarGeneAlterations, fig.width=10, fig.height=7------------------
dnarGeneAlt <- matrix(0, nrow=nrow(dnarGeneExp), ncol=ncol(dnarGeneExp))
rownames(dnarGeneAlt) <- rownames(dnarGeneExp)
colnames(dnarGeneAlt) <- colnames(dnarGeneExp)
dnarGeneAlt[which(dnarGeneMut == 1)] <- 1
dnarGeneAlt[which(dnarGeneExpKo == 1)] <- 2

# Identify genes with most frequent alterations (by mutation or lack of expression).
numAltLines <- apply(dnarGeneAlt, MARGIN = 1, FUN = sum)
sort(numAltLines, decreasing = TRUE)

dnarGeneAltPlotData <- dnarGeneAlt[order(numAltLines), ]
heatmap(dnarGeneAltPlotData, scale="none", Rowv = NA, Colv = NA, 
				col = c("grey", "red", "blue"), 
				main="Altered genes by mutation (red) or low expression (blue)")

## ----checkAplfDrugAssocs, fig.width=10, fig.height=7---------------------
geneName <- "SLX4"
altLineIndices <- which(dnarGeneAlt[geneName, ] != 0)
nonAltLineIndices <- which(dnarGeneAlt[geneName, ] == 0)

drugAssocTab <- data.frame(GENE=geneName, DRUG_NAME=fdaDrugAnnot$NAME, 
													 DRUG_NSC=fdaDrugAnnot$NSC, TTEST_PVAL=NA, 
													 ADJ_PVAL=NA, MEAN_ACT_DIFF=NA, 
													 stringsAsFactors = FALSE)
rownames(drugAssocTab) <- drugAssocTab$DRUG_NSC
for (drugNsc in drugAssocTab$DRUG_NSC){
	ttestResult <- t.test(x=nci60FdaDrugAct[drugNsc, altLineIndices], 
												y=nci60FdaDrugAct[drugNsc, nonAltLineIndices])
	drugAssocTab[drugNsc, "TTEST_PVAL"] <- ttestResult$p.value
	drugAssocTab[drugNsc, "MEAN_ACT_DIFF"] <- 
		ttestResult$estimate[1] - ttestResult$estimate[2]
}

drugAssocTab$ADJ_PVAL <- p.adjust(drugAssocTab$TTEST_PVAL, method = "bonferroni")
drugAssocTab <- drugAssocTab[order(drugAssocTab$ADJ_PVAL), ]

drugAssocTab[drugAssocTab$ADJ_PVAL < 0.05, 
						 c("GENE", "DRUG_NAME", "DRUG_NSC", "ADJ_PVAL", "MEAN_ACT_DIFF")]

meanActDiff <- drugAssocTab$MEAN_ACT_DIFF
negLog10Pval <- -log10(drugAssocTab$TTEST_PVAL)
plot(x=meanActDiff, y=negLog10Pval, xlim = c((min(meanActDiff)-1), (max(meanActDiff)+1)),
		 xlab="Difference between -logGI50 means (altered vs. non-altered lines)",
		 ylab="-log10(p-value)", main=(paste0(geneName, " Drug Response Associations")))

ptLabels <- character(length(drugAssocTab$DRUG_NAME))
numLabeledPts <- 7 # label the k most significant drug associations
ptLabels[1:numLabeledPts] <- drugAssocTab$DRUG_NAME[1:numLabeledPts]
text(x=meanActDiff, negLog10Pval, ptLabels, cex=0.9, pos=4, col="red")

## ----comparePlots, eval=FALSE--------------------------------------------
## runShinyApp("shinyComparePlots")

## ----compoundBrowser, eval=FALSE-----------------------------------------
## runShinyApp("shinyReprodPlots")

## ----compareStructures, eval=FALSE---------------------------------------
## runShinyApp("shinyCompareStructures")

## ----sessionInfo---------------------------------------------------------
sessionInfo()

