test_that("getSmiles() returns correct information", {
	expect_identical(unname(getSmiles("7359")), "Nc1ncnc2c1ncn2C3OC(CO)C(O)C3O")
	expect_identical(unname(getSmiles("666307")), "N#Cc1cccc2Nc3ccccc3Sc12")
	expect_identical(unname(getSmiles("700906")), 
		"CCCNC(=O)c1cc(NC(=O)c2cc(NC(=O)CCS(=O)(=O)OC)cn2C)cn1C")
	expect_identical(unname(getSmiles("94600")), 
									 "CC[C@@]1(O)C(=O)OCC2=C1C=C3N(Cc4cc5ccccc5nc34)C2=O")
	expect_identical(unname(getSmiles("778665")), 
									 "COc1cc([nH]c1\\C=C\\2/N=C(C)C(=C2C)C(=O)CCCCCCCCC(=O)O[C@H]3CC[C@H]4C5CCc6cc(O)ccc6C5CC[C@]34C)c7ccc[nH]7")
	expect_true(is.na(getSmiles("NOT_IN_DB")))
	
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	expect_identical(rownames(drugAnnot), as.character(drugAnnot$NSC))
	nscToSmiles <- getSmiles(rownames(drugAnnot[20000:20861, ]))
	expect_identical(names(nscToSmiles), rownames(drugAnnot[20000:20861, ]))
	expect_identical(unname(nscToSmiles), as.character(drugAnnot[20000:20861, "SMILES"]))
})