test_that("getDrugActivityData output is computed properly", {
	drugAct <- exprs(getAct(rcellminerData::drugData))
	nscSet <- c("141540", "123127", "609699")
	
	actData <- getDrugActivityData(nscSet)
	
	for (nsc in nscSet){
		negLogGI50Data <- as.numeric(actData[nsc, ])
		zScoreData <- as.numeric(drugAct[nsc, ])
		
		corVal <- cor.test(negLogGI50Data, zScoreData)$estimate
		expect_true(corVal > 0.95)
	}
	
})