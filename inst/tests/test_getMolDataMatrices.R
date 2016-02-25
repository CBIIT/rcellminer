test_that("getMolDataMatrices() returns properly structured results", {
	molDataMats <- getMolDataMatrices()
	expect_true(all(c("cop", "exp", "xai", "exo", "mut", "mir", "pro", "mda") 
									%in% names(molDataMats)))
	
	nci60Names <- loadNciColorSet(returnDf = TRUE)$abbrCellLines
	
	# Make sure all matrices have the same (nci60Names) column names.
	expect_true(all(vapply(molDataMats, 
												 function(X) { identical(colnames(X), nci60Names) },
												 logical(1))))
	
	# Make sure feature names (rownames) of each matrix in molDataMats are
	# prefixed with the appropriate molecular data type (derived from the
	# names of the molDataMats list).
	for (molDataType in names(molDataMats)){
		featurePrefix <- unique(substr(rownames(molDataMats[[molDataType]]), 
														 start = 1, stop = nchar(molDataType)))
		expect_identical(molDataType, featurePrefix)
	}
})