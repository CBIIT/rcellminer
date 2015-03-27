test_that("getMedSenLineActivity output is computed properly", {
  idSet <- c("609699", "740")
  medSenActVals <- getMedSenLineActivity(idSet)
  
  expectedOutput <- c(7.852208547, 7.495352901)
  names(expectedOutput) <- idSet
  
  expect_equivalent(medSenActVals, expectedOutput)
})