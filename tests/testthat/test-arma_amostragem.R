testthat::test_that("amostrar_coef_validos retorna coeficientes finitos", {
  set.seed(42)

  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)

  out <- amostrar_coef_validos(coef, vc)

  testthat::expect_false(is.null(out))
  testthat::expect_true(all(is.finite(out)))
  testthat::expect_true(all(c("phi_star", "theta_star") %in% names(out)))
})

testthat::test_that("amostrar_coef_validos retorna NULL para matriz invalida", {
  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- matrix(c(1, NA, 0, 1), nrow = 2)

  out <- amostrar_coef_validos(coef, vc)

  testthat::expect_null(out)
})

testthat::test_that("amostrar_coef_validos retorna NULL para coeficientes ausentes ou invalidos", {
  vc <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)

  out1 <- amostrar_coef_validos(c(ar1 = NA, ma1 = 0.5), vc)
  out2 <- amostrar_coef_validos(c(ar1 = 0.2, ma1 = Inf), vc)

  testthat::expect_null(out1)
  testthat::expect_null(out2)
})

testthat::test_that("amostrar_coef_validos retorna valores dentro de (-1, 1)", {
  set.seed(123)

  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)

  out <- amostrar_coef_validos(coef, vc)

  testthat::expect_true(abs(out["phi_star"]) < 1)
  testthat::expect_true(abs(out["theta_star"]) < 1)
})

testthat::test_that("amostrar_coef_validos trunca coeficientes na fronteira antes da transformacao", {
  coef <- c(ar1 = 1, ma1 = -1)
  vc <- diag(c(0.01, 0.01))

  mockery::stub(
    amostrar_coef_validos,
    "chol",
    function(x) diag(2)
  )

  out <- amostrar_coef_validos(coef, vc, eps = 1e-6)

  testthat::expect_false(is.null(out))
  testthat::expect_true(all(is.finite(out)))
  testthat::expect_true(abs(out["phi_star"]) < 1)
  testthat::expect_true(abs(out["theta_star"]) < 1)
})

testthat::test_that("amostrar_coef_validos retorna NULL quando chol falha", {
  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)

  mockery::stub(
    amostrar_coef_validos,
    "chol",
    function(x) stop("matriz nao positiva definida")
  )

  out <- amostrar_coef_validos(coef, vc)

  testthat::expect_null(out)
})

testthat::test_that("amostrar_coef_validos retorna NULL quando matriz_vcov nao e matriz", {
  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- c(0.01, 0, 0, 0.01)

  out <- amostrar_coef_validos(coef, vc)

  testthat::expect_null(out)
})

testthat::test_that("amostrar_coef_validos retorna NULL quando nomes esperados nao existem", {
  vc <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)
  coef <- c(phi = 0.2, theta = 0.5)

  out <- amostrar_coef_validos(coef, vc)

  testthat::expect_null(out)
})

testthat::test_that("amostrar_coef_validos retorna NULL ou falha de forma previsivel para matriz com dimensao invalida", {
  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- diag(3)

  testthat::expect_error(
    amostrar_coef_validos(coef, vc)
  )
})

testthat::test_that("amostrar_coef_validos e reprodutivel com seed fixa", {
  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)

  set.seed(123)
  out1 <- amostrar_coef_validos(coef, vc)

  set.seed(123)
  out2 <- amostrar_coef_validos(coef, vc)

  testthat::expect_equal(out1, out2)
})

testthat::test_that("amostrar_coef_validos constroi Sigma_u corretamente", {
  coef <- c(ar1 = 0.2, ma1 = 0.5)
  vc <- matrix(c(0.01, 0.002, 0.002, 0.03), nrow = 2, byrow = TRUE)

  phi <- 0.2
  theta <- 0.5

  J_esperada <- diag(c(
    1 / (1 - phi^2),
    1 / (1 - theta^2)
  ))

  Sigma_esperada <- J_esperada %*% vc %*% t(J_esperada)
  Sigma_esperada <- as.matrix(Sigma_esperada)
  rownames(Sigma_esperada) <- colnames(Sigma_esperada) <- c("ar1", "ma1")

  mockery::stub(
    amostrar_coef_validos,
    "chol",
    function(x) {
      testthat::expect_equal(x, Sigma_esperada, tolerance = 1e-12)
      diag(2)
    }
  )

  set.seed(1)
  out <- amostrar_coef_validos(coef, vc)

  testthat::expect_false(is.null(out))
})
