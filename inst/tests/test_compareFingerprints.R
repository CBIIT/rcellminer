test_that("compareFingerprints() output is correct", {
	topotecanNsc <- "609699"
	top1InhNscs <- c(topotecanNsc, setdiff(getMoaToCompounds()[["T1"]], topotecanNsc))
	results <- compareFingerprints(top1InhNscs, getSmiles(top1InhNscs), verbose = FALSE)
	top3Hits <- names(results)[2:4]
	
  expect_identical(names(results)[1], "609699")
	expect_identical(unname(getDrugName(top3Hits)), 
									 c("Topotecan", "Camptothecin Derivative", "Camptothecin Derivative"))
})