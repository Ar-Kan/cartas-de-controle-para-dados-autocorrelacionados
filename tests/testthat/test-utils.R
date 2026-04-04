testthat::test_that("tryNull retorna valor da expressao quando nao ha erro", {
  out <- tryNull(1 + 1)
  testthat::expect_equal(out, 2)
})

testthat::test_that("tryNull retorna NULL quando ha erro", {
  out <- tryNull(1 + "1")
  testthat::expect_null(out)
})
