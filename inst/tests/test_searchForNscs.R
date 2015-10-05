test_that("searchForNscs() returns all NSCs for a given pattern", {
	expect_true(length(searchForNscs("nib$")) > 5)
	expect_equal(length(searchForNscs("TEST$")), 0)
})