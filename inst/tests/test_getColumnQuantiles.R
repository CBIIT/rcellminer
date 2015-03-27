test_that("getColumnQuantiles output is correct", {
	
	X <- matrix(as.numeric(1:25), nrow = 5)
	expect_identical(X[1, ], getColumnQuantiles(X, prob = 0))
	expect_identical(X[3, ], getColumnQuantiles(X, prob = 0.5))
	expect_identical(X[5, ], getColumnQuantiles(X, prob = 1))
	
	Y <- rbind(X, matrix(NA, nrow = 1, ncol = 5))
	expect_identical(X[1, ], getColumnQuantiles(Y, prob = 0, naRm = TRUE))
	expect_identical(X[3, ], getColumnQuantiles(Y, prob = 0.5, naRm = TRUE))
	expect_identical(X[5, ], getColumnQuantiles(Y, prob = 1, naRm = TRUE))
	
	Z <- rbind(Y, matrix(0, nrow = 1, ncol = 5))
	expect_identical(X[1, ], getColumnQuantiles(Z, prob = 0, naRm = TRUE, onlyNonzeroVals = TRUE))
	expect_identical(X[3, ], getColumnQuantiles(Z, prob = 0.5, naRm = TRUE, onlyNonzeroVals = TRUE))
	expect_identical(X[5, ], getColumnQuantiles(Z, prob = 1, naRm = TRUE, onlyNonzeroVals = TRUE))

})