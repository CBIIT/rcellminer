test_that("getDrugMoaList returns all applicable MOA strings", {
	expect_identical(getDrugMoaList("609699"), "TOP1")
	expect_identical(getDrugMoaList("123127"), "TOP2")
	expect_identical(sort(getDrugMoaList("762")), c("A7", "AlkAg"))
	
	expect_null(getDrugMoaList("0"))
})