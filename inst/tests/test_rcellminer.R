# All tests are done on files in package using system.file()

context("RCellmine")

test_that("something", {
  #results <- read.table(system.file("extdata", "10523676-compact.xml", package="paxtoolsr"))
	expect_equal(1, 1)
})

# test_that("something2", {
# 	#results <- read.table(system.file("extdata", "10523676-compact.xml", package="paxtoolsr"))
# 	expect_equal(1, 2)
# })

#DEBUG 
#test_that("FAIL", {    
#   expect_that(FALSE, is_true())
#})
