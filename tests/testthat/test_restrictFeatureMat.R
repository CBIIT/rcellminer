test_that("restrictFeatureMat output is correct", {
	
	X <- matrix(1:25, nrow=5)
	rownames(X) <- c("expA", "expB", "copC", "mutC", "expD")
	Xrestricted <- restrictFeatureMat(geneSet = c("B", "C", "NOT_THERE"), X)
	
	expect_identical(X[c("expB", "copC", "mutC"), ], Xrestricted)
	 
})