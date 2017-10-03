test_that("getFingerprintList() output is correct", {
	topotecanNsc <- "609699"
	camptothecinNsc <- "94600"
	ids <- c(topotecanNsc, camptothecinNsc)
	
	results <- getFingerprintList(ids, getSmiles(ids), verbose = FALSE)
	expect_true(distance(results[[topotecanNsc]], results[[camptothecinNsc]]) > 0.9)
})