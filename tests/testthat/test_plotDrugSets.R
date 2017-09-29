test_that("plotDrugSets", {
	drugAct <- exprs(getAct(rcellminerData::drugData))
	drugs <- rownames(drugAct)[1:8]
	
	file <- tempfile() 
	
	pdf(file)
	plotDrugSets(drugAct, drugs, "Test")
	dev.off()
	
	expect_true(file.exists(file))
})
