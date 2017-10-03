test_that("getRsd", {
	A <- matrix(rnorm(10*60), nrow=10)

	expect_more_than(getRsd(A), 0)
	expect_equal(length(getRsd(A, onlyReturnMedian=FALSE)), 60)
})
