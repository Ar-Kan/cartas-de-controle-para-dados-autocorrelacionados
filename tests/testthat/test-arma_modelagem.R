testthat::test_that("simula_arma chama arima.sim e retorna as ultimas n observacoes", {
  mockery::stub(
    simula_arma,
    "arima.sim",
    function(n, model) {
      testthat::expect_equal(n, 250)
      testthat::expect_equal(model, list(ar = 0.2, ma = 0.5))
      1:250
    }
  )

  x <- simula_arma(n = 50, phi = 0.2, theta = 0.5)

  testthat::expect_true(is.numeric(x))
  testthat::expect_length(x, 50)
  testthat::expect_equal(x, 201:250)
})

testthat::test_that("fit_arma retorna estrutura esperada para serie valida", {
  set.seed(42)
  serie <- simula_arma(n = 200, phi = 0.2, theta = 0.5)

  ajuste <- fit_arma(serie)

  testthat::expect_type(ajuste, "list")
  testthat::expect_true(all(c("fit", "warnings", "convergiu", "erro") %in% names(ajuste)))
  testthat::expect_false(ajuste$erro)
  testthat::expect_false(is.null(ajuste$fit))
})

testthat::test_that("fit_arma falha de forma controlada com entrada invalida", {
  ajuste <- fit_arma(c(1, 2))

  testthat::expect_true(ajuste$erro)
  testthat::expect_false(ajuste$convergiu)
  testthat::expect_null(ajuste$fit)
})

# Sucesso sem warnings
testthat::test_that("fit_arma retorna sucesso quando arima ajusta corretamente", {
  fit_fake <- structure(
    list(),
    class = "Arima"
  )

  mockery::stub(
    fit_arma,
    "arima",
    function(serie, order, include.mean, transform.pars, method, init,
             optim.method, optim.control) {
      testthat::expect_equal(order, c(1, 0, 1))
      testthat::expect_false(include.mean)
      testthat::expect_equal(method, "CSS-ML")
      testthat::expect_equal(init, c(ar1 = 0.2, ma1 = 0.5))
      fit_fake
    }
  )

  mockery::stub(
    fit_arma,
    "coef",
    function(object) c(ar1 = 0.2, ma1 = 0.5)
  )

  mockery::stub(
    fit_arma,
    "vcov",
    function(object) {
      matrix(
        c(0.1, 0, 0, 0.1),
        nrow = 2,
        dimnames = list(c("ar1", "ma1"), c("ar1", "ma1"))
      )
    }
  )

  ajuste <- fit_arma(serie = 1:200, phi = 0.2, theta = 0.5)

  testthat::expect_false(ajuste$erro)
  testthat::expect_true(ajuste$convergiu)
  testthat::expect_identical(ajuste$warnings, character())
  testthat::expect_identical(ajuste$fit, fit_fake)
})

testthat::test_that("fit_arma captura warnings e marca convergiu = FALSE", {
  fit_fake <- structure(list(), class = "Arima")

  mockery::stub(
    fit_arma,
    "arima",
    function(...) {
      warning("possible convergence problem: optim gave code = 1")
      fit_fake
    }
  )

  mockery::stub(
    fit_arma,
    "coef",
    function(object) c(ar1 = 0.2, ma1 = 0.5)
  )

  mockery::stub(
    fit_arma,
    "vcov",
    function(object) diag(2)
  )

  ajuste <- fit_arma(serie = 1:200)

  testthat::expect_false(ajuste$erro)
  testthat::expect_false(ajuste$convergiu)
  testthat::expect_length(ajuste$warnings, 1)
  testthat::expect_match(ajuste$warnings[[1]], "possible convergence problem")
  testthat::expect_identical(ajuste$fit, fit_fake)
})

testthat::test_that("fit_arma retorna erro controlado quando arima falha", {
  mockery::stub(
    fit_arma,
    "arima",
    function(...) {
      stop("series is too short")
    }
  )

  ajuste <- fit_arma(serie = c(1, 2))

  testthat::expect_true(ajuste$erro)
  testthat::expect_false(ajuste$convergiu)
  testthat::expect_null(ajuste$fit)
  testthat::expect_identical(ajuste$warnings, character())
})

testthat::test_that("fit_arma retorna erro quando coef(fit) e invalido", {
  fit_fake <- structure(list(), class = "Arima")

  mockery::stub(
    fit_arma,
    "arima",
    function(...) fit_fake
  )

  mockery::stub(
    fit_arma,
    "coef",
    function(object) NULL
  )

  mockery::stub(
    fit_arma,
    "vcov",
    function(object) diag(2)
  )

  ajuste <- fit_arma(serie = 1:200)

  testthat::expect_true(ajuste$erro)
  testthat::expect_false(ajuste$convergiu)
  testthat::expect_null(ajuste$fit)
})

testthat::test_that("fit_arma retorna erro quando vcov(fit) e invalido", {
  fit_fake <- structure(list(), class = "Arima")

  mockery::stub(
    fit_arma,
    "arima",
    function(...) fit_fake
  )

  mockery::stub(
    fit_arma,
    "coef",
    function(object) c(ar1 = 0.2, ma1 = 0.5)
  )

  mockery::stub(
    fit_arma,
    "vcov",
    function(object) NULL
  )

  ajuste <- fit_arma(serie = 1:200)

  testthat::expect_true(ajuste$erro)
  testthat::expect_false(ajuste$convergiu)
  testthat::expect_null(ajuste$fit)
})
