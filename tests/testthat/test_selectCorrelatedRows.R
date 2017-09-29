test_that("selectCorrelatedRows() returns correct information", {
	X <- getMolDataMatrices()[["exp"]][(1:2000), ]
	corThreshold <- 0.1
	
	# Y is vector, useAbsCor=TRUE
	Y <- as.numeric(getDrugActivityData("609699"))
	corRowDat <- selectCorrelatedRows(Y, X, corThreshold, useAbsCor = TRUE)
	xCorDat <- crossCors(Y, corRowDat)
	expect_true(all(abs(as.numeric(xCorDat$cor)) > corThreshold))
	xCorDat <- crossCors(Y, X[setdiff(rownames(X), rownames(corRowDat)), ])
	expect_true(all(abs(as.numeric(xCorDat$cor)) <= corThreshold))
	
	# Y is vector, useAbsCor=FALSE
	Y <- as.numeric(getDrugActivityData("609699"))
	corRowDat <- selectCorrelatedRows(Y, X, corThreshold, useAbsCor = FALSE)
	xCorDat <- crossCors(Y, corRowDat)
	expect_true(all(as.numeric(xCorDat$cor) > corThreshold))
	xCorDat <- crossCors(Y, X[setdiff(rownames(X), rownames(corRowDat)), ])
	expect_true(all(as.numeric(xCorDat$cor) <= corThreshold))
	
	# Y is matrix, useAbsCor=TRUE
	Y <- getDrugActivityData(c("609699", "760766"))
	corRowDat <- selectCorrelatedRows(Y, X, corThreshold, useAbsCor = TRUE)
	xCorDat <- crossCors(corRowDat, Y)
	maxAbsCors <- apply(abs(xCorDat$cor), MARGIN = 1, max)
	expect_true(all(maxAbsCors > corThreshold))
	xCorDat <- crossCors(X[setdiff(rownames(X), rownames(corRowDat)), ], Y)
	maxAbsCors <- apply(abs(xCorDat$cor), MARGIN = 1, max)
	expect_true(all(maxAbsCors <= corThreshold))
	
	# Y is matrix, useAbsCor=FALSE
	Y <- getDrugActivityData(c("609699", "760766"))
	corRowDat <- selectCorrelatedRows(Y, X, corThreshold, useAbsCor = FALSE)
	xCorDat <- crossCors(corRowDat, Y)
	maxCors <- apply(xCorDat$cor, MARGIN = 1, max)
	expect_true(all(maxCors > corThreshold))
	xCorDat <- crossCors(X[setdiff(rownames(X), rownames(corRowDat)), ], Y)
	maxCors <- apply(xCorDat$cor, MARGIN = 1, max)
	expect_true(all(maxCors <= corThreshold))
})