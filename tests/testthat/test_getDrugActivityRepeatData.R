test_that("getDrugActivityRepeatData returns the correct output", {
	nsc <- "609699"
	drugRepeatAct <- exprs(getRepeatAct(rcellminerData::drugData))
	drugRepeatAnnot <- as(featureData(getRepeatAct(rcellminerData::drugData)), "data.frame")
	
	actData <- getDrugActivityRepeatData(nsc, onlyCellMinerExps = FALSE)
	actDataCm <- getDrugActivityRepeatData(nsc)
	
	iNsc <- which(drugRepeatAnnot$nsc == nsc)
	tmp <- drugRepeatAct[iNsc, ]
	rownames(tmp) <- paste0(nsc, "_", seq(nrow(tmp)))
	expect_identical(tmp, actData)
	
	iNsc <- iNsc[drugRepeatAnnot$used_in_zscore[iNsc]]
	tmp <- drugRepeatAct[iNsc, ]
	rownames(tmp) <- paste0(nsc, "_", seq(nrow(tmp)))
	expect_identical(tmp, actDataCm)
})