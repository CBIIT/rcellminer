test_that("getActivityRangeStats output is correct", {
	
	nscSet <- c("609699", "740", "721")
	actRangeStats <- getActivityRangeStats(nscSet)

	expect_identical(unname(getMedSenLineActivity(nscSet)), actRangeStats$MedSenLineActivity)
	
	for (nsc in nscSet){
		repeatActData <- getDrugActivityRepeatData(nsc)
		
		if (nrow(repeatActData) > 1){
			avgAct <- colMeans(repeatActData, na.rm = TRUE)
		} else{
			avgAct <- as.numeric(repeatActData)
		}
		
		expect_identical(c(min(avgAct, na.rm = TRUE), median(avgAct, na.rm = TRUE), max(avgAct, na.rm = TRUE)),
										 as.numeric(actRangeStats[nsc, c("MinActivity", "MedActivity", "MaxActivity")]))
	}
	
})