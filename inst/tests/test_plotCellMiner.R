test_that("plotCellMiner", {
	drugAct <- exprs(getAct(rcellminerData::drugData))
	molDataMats <- getMolDataMatrices()
	plotDataMats <- molDataMats[c("exp", "cop", "mut")]
	
	file <- tempfile() 
	
	options(bitmapType="cairo")
	pdf(file)
	plotCellMiner(drugAct, plotDataMats, plots=c("mut", "drug", "cop"), nsc="94600", gene="TP53")
	dev.off()
	
	expect_true(file.exists(file))
})
