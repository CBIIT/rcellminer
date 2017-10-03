test_that("getNumMissingLines returns the correct output", {
	drugAct <- exprs(getAct(rcellminerData::drugData))
	expect_identical(sum(is.na(drugAct["1", ])), unname(getNumMissingLines("1")))
	expect_identical(sum(is.na(drugAct["609699", ])), unname(getNumMissingLines("609699")))
	
	nscSet <- rownames(drugAct)[seq(from = 1, to = nrow(drugAct), by = 200)]
	
	numNA <- vapply(nscSet, function(id) sum(is.na(drugAct[id, ])), integer(1))
	
	expect_identical(numNA, getNumMissingLines(nscSet))
})