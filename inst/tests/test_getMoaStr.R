test_that("getMoaStr() returns all applicable MOA strings", {
	expect_identical(getMoaStr("609699"), "T1")
	expect_identical(getMoaStr("123127"), "T2")
	expect_identical(getMoaStr("762"), "A7,AlkAg")
	
	expect_true(is.na(getMoaStr("0")))
})