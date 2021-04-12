test_that("getNumDrugActivityRepeats returns the correct output", {
	drugAct <- exprs(getAct(rcellminerData::drugData))
	drugRepeatAnnot <- as(featureData(getRepeatAct(rcellminerData::drugData)), "data.frame")
	
	nscSet <- rownames(drugAct)[seq(from = 1, to = nrow(drugAct), by = 200)]
	
	numReps <- integer(length(nscSet))
	names(numReps) <- nscSet
	
	numRepsCm <- integer(length(nscSet))
	names(numRepsCm) <- nscSet
	
	for (nsc in names(numReps)){
		iNsc <- which(as.character(drugRepeatAnnot$NSC) == nsc)
		numReps[nsc] <- length(iNsc)
		
		iNsc <- iNsc[drugRepeatAnnot$USED_IN_ZSCORE[iNsc]]
		numRepsCm[nsc] <- length(iNsc)
	}
	
	expect_identical(numReps, getNumDrugActivityRepeats(nscSet, onlyCellMinerExps = FALSE))
	expect_identical(numRepsCm, getNumDrugActivityRepeats(nscSet))
})