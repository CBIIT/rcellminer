test_that("passRuleOf5", {
	expect_false(passRuleOf5("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O", verbose=FALSE))
	expect_true(passRuleOf5("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O", acceptableViolations=2, verbose=FALSE))
	expect_true(passRuleOf5("C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)(F)F", verbose=FALSE))
	expect_equivalent(passRuleOf5("BADSMILES", verbose=FALSE), NA)
})
