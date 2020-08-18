test_that("getDrugName() returns correct information", {
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	expect_identical(rownames(drugAnnot), drugAnnot$NSC)
	drugAnnotNamed <- drugAnnot[!is.na(drugAnnot$NAME), ]
	drugAnnotUnamed <- drugAnnot[is.na(drugAnnot$NAME), ]
	
	nscToDrugName <- getDrugName(rownames(drugAnnotNamed))
	expect_identical(names(nscToDrugName), rownames(drugAnnotNamed))
	expect_identical(unname(nscToDrugName), drugAnnotNamed$NAME)
	expect_identical(unname(nscToDrugName["609699"]), "Topotecan")
	expect_identical(unname(nscToDrugName["94600"]), "CAMPTOTHECIN")
	
	nscToDrugName <- getDrugName(rownames(drugAnnotUnamed))
	expect_identical(names(nscToDrugName), rownames(drugAnnotUnamed))
	expect_true(all(is.na(nscToDrugName)))
	
	expect_true(is.na(getDrugName("NOT_IN_DB")))
})