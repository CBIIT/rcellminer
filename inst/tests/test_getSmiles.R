test_that("getSmiles() returns correct information", {
	expect_identical(unname(getSmiles("7359")), "Nc1ncnc2c1ncn2C3OC(CO)C(O)C3O")
	expect_identical(unname(getSmiles("666307")), "N#Cc1cccc2Nc3ccccc3Sc12")
	expect_identical(unname(getSmiles("700906")), 
		"CCCNC(=O)c1cc(NC(=O)c2cc(NC(=O)CCS(=O)(=O)OC)cn2C)cn1C")
	expect_identical(unname(getSmiles("744570")), 
		"C[C@@H]1CC[C@H]2[C@@H](C)[C@@H](OCCn3cc(nn3)C(=O)N\\N=C\\c4ccccc4O)O[C@@H]5O[C@@]6(C)CC[C@@H]1C25OO6")
	expect_identical(unname(getSmiles("778665")), 
									 "C1(=Cc2c(cc(c3ccc[nH]3)[nH]2)OC)C(=C(C(=N1)C)C(=O)CCCCCCCCC(=O)OC1C2(C(C3CCc4cc(ccc4C3CC2)O)CC1)C)C.Cl")
	expect_true(is.na(getSmiles("NOT_IN_DB")))
	
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	expect_identical(rownames(drugAnnot), drugAnnot$NSC)
	nscToSmiles <- getSmiles(rownames(drugAnnot[20000:20861, ]))
	expect_identical(names(nscToSmiles), rownames(drugAnnot[20000:20861, ]))
	expect_identical(unname(nscToSmiles), as.character(drugAnnot[20000:20861, "SMILES"]))
})