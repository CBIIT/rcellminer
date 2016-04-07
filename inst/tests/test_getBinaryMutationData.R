test_that("getBinaryMutationData() output is correct", {
# Note: CellMiner 2.0 mutation data is no longer binary.
# 	exoDb <- getESetList(rcellminerData::molData)[["exo"]]
# 	exoInfoTab <- as(featureData(exoDb), "data.frame")
# 	exoDataMat <- exprs(exoDb)
# 	binMutMat <- exprs(getESetList(rcellminerData::molData)[["mut"]])
# 	
# 	funcResult <- getBinaryMutationData(exoInfoTab, exoDataMat)
	
	# NOTE: This test is used because the R CHECK command seems to reorder
	# the rows of the matrices.
# 	expect_true(identical(binMutMat[sort(rownames(binMutMat)), ],
# 												funcResult[sort(rownames(funcResult)), ]))
	
# 	binMutMatVals <- as.vector(binMutMat)
# 	if (length(unique(binMutMatVals)) != 2){
# 		warning("binMutMat data is not 0 or 1 as expected. Check if current rcellminerData package is installed.")
# 	} else{
# 		warning("*****OK******")
# 	}
	
# 	expect_identical(colnames(binMutMat), colnames(funcResult))
# 	
# 	expect_identical(nrow(binMutMat), nrow(funcResult))
# 	
# 	diffRownames <- setdiff(rownames(binMutMat), rownames(funcResult))
# 	if (length(diffRownames) > 0){
# 		for (name in diffRownames){
# 			warning(name)
# 		}
# 	} else{
# 		warning("***ROWNAMES OK(1)*****")
# 	}
# 	
# 	diffRownames <- setdiff(rownames(funcResult), rownames(binMutMat))
# 	if (length(diffRownames) > 0){
# 		for (name in diffRownames){
# 			warning(name)
# 		}
# 	} else{
# 		warning("***ROWNAMES OK(2)*****")
# 	}
# 
# 	expect_identical(sort(rownames(binMutMat)), sort(rownames(funcResult)))
# 	expect_identical(rownames(binMutMat), rownames(funcResult))
	

	
# 	expect_equal(setdiff(as.vector(funcResult), as.vector(binMutMat)), integer(0))
# 	expect_true(identical(as.vector(funcResult), as.vector(binMutMat)))
# 	expect_false(identical(c(1, as.vector(funcResult)), as.vector(binMutMat)))
})