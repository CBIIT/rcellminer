test_that("plotDrugActivityRepeats", {
	file <- tempfile() 
	
	pdf(file)
	plotDrugActivityRepeats("609699")
	dev.off()
	
	expect_true(file.exists(file))
	
	file <- tempfile() 
	
	pdf(file)
	plotDrugActivityRepeats("609699", useZScore=TRUE, maxRepNum=3)
	dev.off()
	
	expect_true(file.exists(file))
})
