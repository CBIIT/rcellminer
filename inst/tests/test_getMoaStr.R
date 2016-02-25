test_that("getMoaStr() returns all applicable MOA strings", {
	expect_identical(getMoaStr("609699"), c("609699"="T1"))
	expect_identical(getMoaStr("123127"), c("123127"="T2"))
	expect_identical(getMoaStr("762"), c("762"="A7,AlkAg"))
	
	expect_true(is.na(getMoaStr("0")))
	expect_true(any(is.na(getMoaStr(c("94600", "609699", "TEST")))))
})
