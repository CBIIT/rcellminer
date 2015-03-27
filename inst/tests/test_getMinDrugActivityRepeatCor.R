test_that("getDrugName() returns correct information", {
	testNscSet <- c("609699", "94600")
	
	funcResultsTab <- getMinDrugActivityRepeatCor(testNscSet)
	
	for (nsc in testNscSet){
		repData <- getDrugActivityRepeatData(nsc)
		xCorDat <- crossCors(repData)
		iMax <- which(xCorDat$pval == max(xCorDat$pval, na.rm = TRUE), arr.ind=TRUE)[1, ]
		expect_identical(funcResultsTab[nsc, "cor"], 
										 round(xCorDat$cor[iMax["row"], iMax["col"]], digits = 3))
		expect_identical(funcResultsTab[nsc, "pval"], 
										 signif(xCorDat$pval[iMax["row"], iMax["col"]], digits = 3))
	}
})