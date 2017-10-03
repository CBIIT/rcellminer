test_that("getMolDataType output is correct", {
	
	expect_identical(unname(getMolDataType(c("expTP53", "copMDM2", "mutCHEK2", "mutBRAF"))),
									 c("exp", "cop", "mut", "mut"))
	
	expect_identical(unname(getMolDataType("expSLFN11")), "exp")
	
})