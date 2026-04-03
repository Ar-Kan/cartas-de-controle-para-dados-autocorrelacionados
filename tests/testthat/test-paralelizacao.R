source("R/paralelizacao.R") |> try()

testthat::test_that("faz_uma_execucao_safe retorna resultado com coluna execucao", {

  f <- function(execucao, a) {
    data.table::data.table(valor = execucao + a)
  }

  out <- faz_uma_execucao_safe(
    execucao = 3,
    funcao = f,
    a = 2
  )

  testthat::expect_true(data.table::is.data.table(out))
  testthat::expect_true("execucao" %in% names(out))
  testthat::expect_equal(out$execucao, 3)
  testthat::expect_equal(out$valor, 5)
})


testthat::test_that("faz_uma_execucao_safe preserva execucao se já existir", {

  f <- function(execucao) {
    data.table::data.table(execucao = 999, valor = 1)
  }

  out <- faz_uma_execucao_safe(
    execucao = 3,
    funcao = f
  )

  testthat::expect_equal(out$execucao, 999)
  testthat::expect_equal(out$valor, 1)
})


testthat::test_that("faz_uma_execucao_safe converte erro em data.table padronizada", {

  f <- function(execucao) {
    stop("falha proposital")
  }

  out <- faz_uma_execucao_safe(
    execucao = 2,
    funcao = f
  )

  testthat::expect_true(data.table::is.data.table(out))
  testthat::expect_true("erro_worker" %in% names(out))
  testthat::expect_equal(out$execucao, 2)
  testthat::expect_match(out$erro_worker, "falha proposital")
  testthat::expect_true("pid" %in% names(out))
})


testthat::test_that("faz_uma_execucao_safe retorna data.table vazia quando funcao retorna NULL", {

  f <- function(execucao) {
    NULL
  }

  out <- faz_uma_execucao_safe(
    execucao = 1,
    funcao = f
  )

  testthat::expect_true(data.table::is.data.table(out))
  testthat::expect_equal(nrow(out), 0L)
})


testthat::test_that("consolida_resultados_paralelos junta resultados válidos", {
  res_lista <- list(
    data.table::data.table(execucao = 1, valor = 10),
    data.table::data.table(execucao = 2, valor = 20)
  )

  out <- consolida_resultados_paralelos(res_lista)

  testthat::expect_true(data.table::is.data.table(out))
  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_equal(out$valor, c(10, 20))
})


testthat::test_that("consolida_resultados_paralelos ignora erros e mantém válidos", {
  res_lista <- list(
    data.table::data.table(execucao = 1, valor = 10),
    data.table::data.table(execucao = 2, erro_worker = "erro", pid = 123)
  )

  testthat::expect_warning(
    out <- consolida_resultados_paralelos(res_lista),
    "Ocorreram erros"
  )

  testthat::expect_true(data.table::is.data.table(out))
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$execucao, 1)
  testthat::expect_equal(out$valor, 10)
})

testthat::test_that("consolida_resultados_paralelos emite warning quando há erros", {
  res_lista <- list(
    data.table::data.table(execucao = 1, valor = 10),
    data.table::data.table(execucao = 2, erro_worker = "erro", pid = 123)
  )

  testthat::expect_warning(
    consolida_resultados_paralelos(res_lista),
    "Ocorreram erros"
  )
})

testthat::test_that("consolida_resultados_paralelos falha sem resultados válidos", {
  res_lista <- list(
    data.table::data.table(execucao = 1, erro_worker = "erro", pid = 111),
    data.table::data.table(execucao = 2, erro_worker = "erro", pid = 222)
  )

  testthat::expect_error(
    testthat::expect_warning(
      consolida_resultados_paralelos(res_lista),
      "Ocorreram erros"
    ),
    "Nenhum resultado foi gerado"
  )
})


testthat::test_that("execucao_paralela executa função em paralelo e consolida", {

  f <- function(execucao, multiplicador) {
    data.table::data.table(
      execucao = execucao,
      valor = execucao * multiplicador
    )
  }

  out <- execucao_paralela(
    n_execucoes = 4,
    funcao = f,
    lista_argumentos = list(multiplicador = 2),
    n_processos = 2,
    seed = TRUE
  )

  data.table::setorder(out, execucao)

  testthat::expect_true(data.table::is.data.table(out))
  testthat::expect_equal(nrow(out), 4L)
  testthat::expect_equal(out$execucao, 1:4)
  testthat::expect_equal(out$valor, c(2, 4, 6, 8))
})


testthat::test_that("execucao_paralela mantém resultados válidos mesmo com falhas em alguns workers", {

  f <- function(execucao) {
    if (execucao %% 2 == 0) {
      stop("erro par")
    }

    data.table::data.table(
      execucao = execucao,
      valor = execucao
    )
  }

  testthat::expect_warning(
    out <- execucao_paralela(
      n_execucoes = 5,
      funcao = f,
      n_processos = 2,
      seed = TRUE
    ),
    "Ocorreram erros"
  )

  data.table::setorder(out, execucao)

  testthat::expect_equal(out$execucao, c(1, 3, 5))
  testthat::expect_equal(out$valor, c(1, 3, 5))
})
