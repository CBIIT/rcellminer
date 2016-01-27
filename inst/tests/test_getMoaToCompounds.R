test_that("getMoaToCompounds() returns correct information", {
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	expect_identical(rownames(drugAnnot), drugAnnot$NSC)
	
	moaToCmpds <- getMoaToCompounds()
	
	knownMoaDrugs <- drugAnnot[!is.na(drugAnnot$MOA), "NSC"]
	testResults <- logical(length(knownMoaDrugs))
	names(testResults) <- knownMoaDrugs
	
	expect_identical(sort(knownMoaDrugs),
									 sort(unique(c(moaToCmpds, recursive = TRUE))))
	
	for (nsc in knownMoaDrugs){
		moaClasses <- strsplit(drugAnnot[nsc, "MOA"], split = "[|]")[[1]]
		
		testResults[nsc] <- all(vapply(moaClasses, 
																	 function(moa) { nsc %in% moaToCmpds[[moa]] }, 
																	 logical(1)))
	}
	
	expect_true(all(testResults))
})