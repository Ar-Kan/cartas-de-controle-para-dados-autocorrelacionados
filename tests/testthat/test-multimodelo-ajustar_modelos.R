testthat::test_that("graus_de_liberdade gera todas as combinações exceto (0,0)", {
  res <- graus_de_liberdade(max_p = 2, max_q = 2)

  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_named(res, c("p", "q"))

  # (2+1) * (2+1) - 1 = 8
  testthat::expect_equal(nrow(res), 8)

  # não deve conter (0,0)
  testthat::expect_false(any(res$p == 0 & res$q == 0))

  # deve conter algumas combinações esperadas
  testthat::expect_true(any(res$p == 1 & res$q == 0))
  testthat::expect_true(any(res$p == 0 & res$q == 1))
  testthat::expect_true(any(res$p == 2 & res$q == 2))
})

testthat::test_that("graus_de_liberdade funciona para max_p = 0 ou max_q = 0", {
  res1 <- graus_de_liberdade(max_p = 0, max_q = 2)
  testthat::expect_equal(nrow(res1), 2)
  testthat::expect_true(all(res1$p == 0))
  testthat::expect_equal(res1$q, c(1, 2))

  res2 <- graus_de_liberdade(max_p = 2, max_q = 0)
  testthat::expect_equal(nrow(res2), 2)
  testthat::expect_true(all(res2$q == 0))
  testthat::expect_equal(res2$p, c(1, 2))
})

testthat::test_that("calcular_pesos funciona com valores normais", {
  ics <- c(100, 102, 105)

  pesos <- calcular_pesos(ics)

  testthat::expect_length(pesos, 3)
  testthat::expect_true(all(pesos >= 0))
  testthat::expect_equal(sum(pesos), 1, tolerance = 1e-8)
})

testthat::test_that("calcular_pesos ignora Inf corretamente", {
  ics <- c(100, Inf, 102)

  pesos <- calcular_pesos(ics)

  testthat::expect_equal(pesos[2], 0)
  testthat::expect_equal(sum(pesos), 1, tolerance = 1e-8)
})

testthat::test_that("calcular_pesos ignora NA corretamente", {
  ics <- c(100, NA, 102)

  pesos <- calcular_pesos(ics)

  testthat::expect_equal(pesos[2], 0)
  testthat::expect_equal(sum(pesos), 1, tolerance = 1e-8)
})

testthat::test_that("calcular_pesos falha se todos inválidos", {
  testthat::expect_error(
    calcular_pesos(c(Inf, Inf)),
    "Nenhum IC válido"
  )
})

testthat::test_that("modelo_ic retorna o critério correto", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 120)
  modelo <- fit_modelo(serie, p = 1, q = 0)

  testthat::expect_equal(modelo_ic(modelo, "aic"), modelo$aic)
  testthat::expect_equal(modelo_ic(modelo, "aicc"), modelo$aicc)
  testthat::expect_equal(modelo_ic(modelo, "bic"), modelo$bic)
})

testthat::test_that("modelo_ic retorna Inf para modelo NULL", {
  testthat::expect_equal(modelo_ic(NULL, "aic"), Inf)
})

testthat::test_that("modelo_ic falha para critério desconhecido", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 120)
  modelo <- fit_modelo(serie, p = 1, q = 0)

  testthat::expect_error(
    modelo_ic(modelo, "foo"),
    "Critério de informação desconhecido"
  )
})

testthat::test_that("ajustar_modelos retorna estrutura esperada sem cutoff", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 200)

  res <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aic")

  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("serie", "modelos", "meta"))

  testthat::expect_equal(res$serie, serie)
  testthat::expect_type(res$meta, "list")
  testthat::expect_equal(res$meta$criterio, "aic")
  testthat::expect_null(res$meta$cutoff)

  testthat::expect_s3_class(res$modelos, "data.frame")
  testthat::expect_named(
    res$modelos,
    c("p", "q", "modelo", "ic", "convergiu", "peso")
  )

  testthat::expect_equal(nrow(res$modelos), 3)
  testthat::expect_true(all(res$modelos$convergiu))
  testthat::expect_true(all(is.finite(res$modelos$ic)))
  testthat::expect_true(all(res$modelos$peso >= 0))
  testthat::expect_equal(sum(res$modelos$peso), 1, tolerance = 1e-8)
})

testthat::test_that("ajustar_modelos calcula pesos corretamente", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.4), n = 150)

  res <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aic")

  pesos_esperados <- calcular_pesos(res$modelos$ic)

  testthat::expect_equal(res$modelos$peso, pesos_esperados, tolerance = 1e-10)
})

testthat::test_that("ajustar_modelos com cutoff filtra modelos e recalcula pesos", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.6), n = 200)

  res_full <- ajustar_modelos(serie, max_p = 2, max_q = 2, criterio = "aic")
  delta_ic <- res_full$modelos$ic - min(res_full$modelos$ic, na.rm = TRUE)
  idx_cut <- which(delta_ic < 2)

  res_cut <- ajustar_modelos(
    serie,
    max_p = 2,
    max_q = 2,
    criterio = "aic",
    cutoff = 2
  )

  testthat::expect_equal(res_cut$meta$cutoff, 2)
  testthat::expect_equal(nrow(res_cut$modelos), length(idx_cut))
  testthat::expect_true(all(res_cut$modelos$peso >= 0))
  testthat::expect_equal(sum(res_cut$modelos$peso), 1, tolerance = 1e-8)
  testthat::expect_equal(res_cut$modelos$ic, res_full$modelos$ic[idx_cut])

  pesos_esperados <- calcular_pesos(res_cut$modelos$ic)
  testthat::expect_equal(res_cut$modelos$peso, pesos_esperados, tolerance = 1e-10)
})

testthat::test_that("ajustar_modelos respeita o critério informado", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 200)

  res_aic <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aic")
  res_aicc <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aicc")
  res_bic <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "bic")

  testthat::expect_equal(res_aic$meta$criterio, "aic")
  testthat::expect_equal(res_aicc$meta$criterio, "aicc")
  testthat::expect_equal(res_bic$meta$criterio, "bic")

  testthat::expect_length(res_aic$modelos$ic, 3)
  testthat::expect_length(res_aicc$modelos$ic, 3)
  testthat::expect_length(res_bic$modelos$ic, 3)

  # Em geral AIC, AICc e BIC não devem ser idênticos para todos os modelos
  testthat::expect_false(isTRUE(all.equal(res_aic$modelos$ic, res_aicc$modelos$ic)))
  testthat::expect_false(isTRUE(all.equal(res_aic$modelos$ic, res_bic$modelos$ic)))
})

testthat::test_that("ajustar_modelos remove do retorno modelos que falham no ajuste", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 100)

  original_fit_modelo <- fit_modelo
  ajustar_modelos_stubbed <- ajustar_modelos

  mockery::stub(
    where = ajustar_modelos_stubbed,
    what = "fit_modelo",
    how = function(serie, p, q) {
      if (p == 1 && q == 1) {
        stop("erro forçado no ajuste")
      }
      original_fit_modelo(serie, p = p, q = q)
    }
  )

  res <- ajustar_modelos_stubbed(
    serie,
    max_p = 1,
    max_q = 1,
    criterio = "aic"
  )

  testthat::expect_equal(nrow(res$modelos), 2)
  testthat::expect_false(any(res$modelos$p == 1 & res$modelos$q == 1))
  testthat::expect_true(all(res$modelos$convergiu))
  testthat::expect_true(all(is.finite(res$modelos$ic)))

  # Os pesos ainda devem somar 1 considerando os modelos válidos
  testthat::expect_equal(sum(res$modelos$peso), 1, tolerance = 1e-8)
})

testthat::test_that("com cutoff muito pequeno pode retornar apenas o melhor modelo", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.7), n = 250)

  res <- ajustar_modelos(
    serie,
    max_p = 2,
    max_q = 2,
    criterio = "aic",
    cutoff = 1e-12
  )

  testthat::expect_gte(nrow(res$modelos), 1)
  testthat::expect_equal(sum(res$modelos$peso), 1, tolerance = 1e-8)

  # Se retornar um único modelo, o peso deve ser 1
  if (nrow(res$modelos) == 1) {
    testthat::expect_equal(res$modelos$peso, 1)
  }
})

testthat::test_that("ajustar_modelos suporta AICc", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 120)

  res <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aicc")

  testthat::expect_equal(res$meta$criterio, "aicc")
  testthat::expect_equal(nrow(res$modelos), 3)
  testthat::expect_true(all(is.finite(res$modelos$ic)))

  # valida se bate com o valor direto do modelo
  modelo_ref <- forecast::Arima(serie, order = c(1, 0, 0), include.mean = FALSE)
  idx <- which(res$modelos$p == 1 & res$modelos$q == 0)

  testthat::expect_length(idx, 1)
  testthat::expect_equal(res$modelos$ic[idx], modelo_ref$aicc, tolerance = 1e-6)
})