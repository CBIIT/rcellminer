test_that("isPublic", {
  expect_false(isPublic("746998"))
  expect_false(all(isPublic(c("746998", "609699"))))
})