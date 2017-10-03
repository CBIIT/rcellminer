test_that("getDrugActivityRange output is correct", {
	
	nscSet <- c("609699", "740", "721")
	logActRanges <- getDrugActivityRange(nscSet)
	logActIQRs <- getDrugActivityRange(nscSet, computeIQR = TRUE)
	
	for (nsc in nscSet){
		avgActVec <- as.numeric(getDrugActivityData(nsc))
		expect_identical((max(avgActVec, na.rm = TRUE) - min(avgActVec, na.rm = TRUE)), 
										 unname(logActRanges[nsc]))
		expect_identical(IQR(avgActVec, na.rm = TRUE), unname(logActIQRs[nsc]))
	}
	
	suppressWarnings(expect_true(is.na(getDrugActivityRange("-1"))))
	suppressWarnings(expect_true(is.na(getDrugActivityRange("-1", computeIQR = TRUE))))
	
})