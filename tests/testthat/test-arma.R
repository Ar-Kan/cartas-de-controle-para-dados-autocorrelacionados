testthat::test_that("ar_valido identifica parametros validos e invalidos", {
  testthat::expect_true(ar_valido(0.5))
  testthat::expect_false(ar_valido(1.2))
  testthat::expect_false(ar_valido(Inf))
})

testthat::test_that("ma_valido identifica parametros validos e invalidos", {
  testthat::expect_true(ma_valido(0.5))
  testthat::expect_false(ma_valido(1.2))
  testthat::expect_false(ma_valido(NaN))
})

testthat::test_that("cola_series preserva o comprimento da fase I", {
  out <- cola_series(
    serie_fase1 = 1:5,
    serie_fase2 = 101:110,
    numero_de_novas_observacoes = 2
  )

  testthat::expect_equal(out, c(3, 4, 5, 101, 102))
  testthat::expect_length(out, 5)
})

testthat::test_that("cola_series usa apenas fase II quando novas observacoes excedem tamanho inicial", {
  out <- cola_series(
    serie_fase1 = 1:5,
    serie_fase2 = 101:110,
    numero_de_novas_observacoes = 5
  )

  testthat::expect_equal(out, 101:105)
})
