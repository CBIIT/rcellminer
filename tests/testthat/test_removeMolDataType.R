test_that("removeMolDataType output is correct", {
	
	expect_identical(unname(removeMolDataType(c("expTP53", "copMDM2", "mutCHEK2", "mutBRAF"))),
									 c("TP53", "MDM2", "CHEK2", "BRAF"))
	
	expect_identical(unname(removeMolDataType("expSLFN11")), "SLFN11")
	
})