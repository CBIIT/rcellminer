test_that("getMedSenLineActivity output is computed properly", {
  idSet <- c("609699", "740")
  medSenActVals <- getMedSenLineActivity(idSet)
  
  expectedOutput <- c(7.851282, 7.493676)
  names(expectedOutput) <- idSet
  
  expect_equal(medSenActVals, expectedOutput, tolerance = 10^-7)
})