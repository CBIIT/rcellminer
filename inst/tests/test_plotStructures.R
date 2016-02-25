test_that("plotStructures", {
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	file <- tempfile() 
	
	pdf(file)
	plotStructures("94600", drugAnnot["94600","SMILES"])
	dev.off()
	
	expect_true(file.exists(file))
	
	file <- tempfile() 
	
	pdf(file)
	plotStructures(c("609699", "94600"), drugAnnot[c("609699", "94600"),"SMILES"])
	dev.off()
	
	expect_true(file.exists(file))
})

