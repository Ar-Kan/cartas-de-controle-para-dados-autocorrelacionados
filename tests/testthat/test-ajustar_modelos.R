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

testthat::test_that("ajustar_modelos retorna estrutura esperada sem cutoff", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 200)

  res <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aic")

  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("modelos", "ics", "pesos", "modelos_df"))

  testthat::expect_s3_class(res$modelos_df, "data.frame")
  testthat::expect_equal(nrow(res$modelos_df), 3) # (0,1), (1,0), (1,1)

  testthat::expect_length(res$modelos, 3)
  testthat::expect_length(res$ics, 3)
  testthat::expect_length(res$pesos, 3)

  testthat::expect_true(all(is.finite(res$ics) | is.infinite(res$ics)))
  testthat::expect_true(all(res$pesos >= 0))
  testthat::expect_equal(sum(res$pesos), 1, tolerance = 1e-8)
})

testthat::test_that("ajustar_modelos calcula pesos corretamente", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.4), n = 150)

  res <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aic")

  delta_ic <- res$ics - min(res$ics, na.rm = TRUE)
  pesos_esperados <- exp(-delta_ic / 2) / sum(exp(-delta_ic / 2), na.rm = TRUE)

  testthat::expect_equal(res$pesos, pesos_esperados, tolerance = 1e-10)
})

testthat::test_that("ajustar_modelos com cutoff filtra modelos e recalcula pesos", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.6), n = 200)

  res_full <- ajustar_modelos(serie, max_p = 2, max_q = 2, criterio = "aic")
  delta_ic <- res_full$ics - min(res_full$ics, na.rm = TRUE)
  idx_cut <- which(delta_ic < 2)

  res_cut <- ajustar_modelos(serie, max_p = 2, max_q = 2, criterio = "aic", cutoff = 2)

  testthat::expect_equal(nrow(res_cut$modelos_df), length(idx_cut))
  testthat::expect_equal(length(res_cut$modelos), length(idx_cut))
  testthat::expect_equal(length(res_cut$ics), length(idx_cut))
  testthat::expect_equal(length(res_cut$pesos), length(idx_cut))

  testthat::expect_true(all(res_cut$pesos >= 0))
  testthat::expect_equal(sum(res_cut$pesos), 1, tolerance = 1e-8)

  # Os ICs retornados devem ser exatamente os filtrados
  testthat::expect_equal(res_cut$ics, res_full$ics[idx_cut])
})

testthat::test_that("ajustar_modelos respeita o criterio informado", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 200)

  res_aic <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aic")
  res_bic <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aicc")

  testthat::expect_length(res_aic$ics, 3)
  testthat::expect_length(res_bic$ics, 3)

  # Em geral AIC e AICc não devem ser idênticos para todos os modelos
  testthat::expect_false(isTRUE(all.equal(res_aic$ics, res_bic$ics)))
})

testthat::test_that("ajustar_modelos trata erro de ajuste atribuindo Inf e NULL", {
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

  idx <- which(res$modelos_df$p == 1 & res$modelos_df$q == 1)

  testthat::expect_length(idx, 1)
  testthat::expect_true(is.infinite(res$ics[idx]))
  testthat::expect_equal(length(res$modelos), nrow(res$modelos_df))
  testthat::expect_null(res$modelos[[idx]])

  # Os pesos ainda devem somar 1 considerando os modelos válidos
  testthat::expect_equal(sum(res$pesos), 1, tolerance = 1e-8)
})

testthat::test_that("com cutoff muito pequeno pode retornar apenas o melhor modelo", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.7), n = 250)

  res <- ajustar_modelos(serie, max_p = 2, max_q = 2, criterio = "aic", cutoff = 1e-12)

  testthat::expect_gte(nrow(res$modelos_df), 1)
  testthat::expect_equal(sum(res$pesos), 1, tolerance = 1e-8)

  # Se retornar um único modelo, o peso deve ser 1
  if (length(res$pesos) == 1) {
    testthat::expect_equal(res$pesos, 1)
  }
})

testthat::test_that("ajustar_modelos suporta AICc", {
  set.seed(123)
  serie <- arima.sim(model = list(ar = 0.5), n = 120)

  res <- ajustar_modelos(serie, max_p = 1, max_q = 1, criterio = "aicc")

  testthat::expect_length(res$ics, 3)
  testthat::expect_true(all(is.finite(res$ics)))

  # valida se bate com o valor direto do modelo
  modelo_ref <- forecast::Arima(serie, order = c(1, 0, 0), include.mean = FALSE)

  idx <- which(res$modelos_df$p == 1 & res$modelos_df$q == 0)

  testthat::expect_equal(res$ics[idx], modelo_ref$aicc, tolerance = 1e-6)
})

testthat::test_that("calcular_pesos funciona com valores normais", {
  ics <- c(100, 102, 105)

  pesos <- calcular_pesos(ics)

  testthat::expect_equal(sum(pesos), 1, tolerance = 1e-8)
  testthat::expect_true(all(pesos >= 0))
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
  testthat::expect_error(calcular_pesos(c(Inf, Inf)))
})
