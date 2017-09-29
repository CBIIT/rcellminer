test_that("patternComparison() returns appropriate results", {
	drugAct <- exprs(getAct(rcellminerData::drugData))
	molDataMats <- getMolDataMatrices()
	patternCompMats <- molDataMats[c("exp", "cop", "mut")]
	
	topotecanActVec <- as.numeric(drugAct["609699", ])
	names(topotecanActVec) <- colnames(drugAct)
	
	topotecanActDf <- as.data.frame(drugAct["609699", , drop = FALSE])
	
	vectorPatternResults <- patternComparison(topotecanActVec, patternCompMats)
	matrixPatternResults <- patternComparison(as.matrix(topotecanActDf), patternCompMats)
	dfPatternResults     <- patternComparison(topotecanActDf, patternCompMats)
	
	singlePatternMatrixResults <- patternComparison(topotecanActVec, patternCompMats$exp)
	singlePatternDfResults <- patternComparison(topotecanActVec, as.data.frame(patternCompMats$exp))
	
	expect_identical(vectorPatternResults, matrixPatternResults)
	expect_identical(vectorPatternResults, dfPatternResults)
	expect_identical(singlePatternMatrixResults, singlePatternDfResults)
	
	expect_true(rownames(vectorPatternResults)[1] == "expSLFN11")
	expect_true(rownames(singlePatternMatrixResults)[1] == "expSLFN11")
})