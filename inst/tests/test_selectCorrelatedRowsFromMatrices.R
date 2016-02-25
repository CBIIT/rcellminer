test_that("selectCorrelatedRowsFromMatrices() returns correct information", {
	XList <- lapply(getMolDataMatrices()[c("exp", "cop")],
									function(x) { x[1:1000, ]} )
	X <- rbind(XList[[1]], XList[[2]])
	corThreshold <- 0.1
	
	# Y is vector, useAbsCor=TRUE
	Y <- as.numeric(getDrugActivityData("609699"))
	names(Y) <- colnames(X)
	expect_identical(
		selectCorrelatedRows(Y, X, corThreshold, useAbsCor = TRUE),
		selectCorrelatedRowsFromMatrices(Y, XList, corThreshold, useAbsCor = TRUE)
	)
	
	# Y is vector, useAbsCor=FALSE
	Y <- as.numeric(getDrugActivityData("609699"))
	names(Y) <- colnames(X)
	expect_identical(
		selectCorrelatedRows(Y, X, corThreshold, useAbsCor = FALSE),
		selectCorrelatedRowsFromMatrices(Y, XList, corThreshold, useAbsCor = FALSE)
	)

	# Y is matrix, useAbsCor=TRUE
	Y <- getDrugActivityData(c("609699", "760766"))
	expect_identical(
		selectCorrelatedRows(Y, X, corThreshold, useAbsCor = TRUE),
		selectCorrelatedRowsFromMatrices(Y, XList, corThreshold, useAbsCor = TRUE)
	)
	
	# Y is matrix, useAbsCor=FALSE
	Y <- getDrugActivityData(c("609699", "760766"))
	expect_identical(
		selectCorrelatedRows(Y, X, corThreshold, useAbsCor = FALSE),
		selectCorrelatedRowsFromMatrices(Y, XList, corThreshold, useAbsCor = FALSE)
	)
})