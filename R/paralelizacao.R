#' Executa uma função de forma segura dentro de um worker
#'
#' @param execucao Identificador da execução.
#' @param funcao Função a ser executada. Deve aceitar `execucao` como argumento nomeado.
#' @param ... Argumentos adicionais repassados para `funcao`.
#'
#' @return O resultado retornado por `funcao`, ou uma data.table com erro.
faz_uma_execucao_safe <- function(execucao, funcao, ..., p = NULL) {
  stopifnot(length(execucao) == 1L, is.numeric(execucao) || is.integer(execucao))
  stopifnot(is.function(funcao))

  on.exit({
    if (!is.null(p)) {
      p()
    }
  }, add = TRUE)

  tryCatch(
  {
    resultado <- funcao(execucao = execucao, ...)

    if (is.null(resultado)) {
      return(data.table::data.table())
    }

    if (data.table::is.data.table(resultado) || is.data.frame(resultado)) {
      if (!("execucao" %in% names(resultado))) {
        resultado$execucao <- execucao
      }
    }

    resultado
  },
    error = function(e) {
      data.table::data.table(
        execucao = execucao,
        erro_worker = as.character(e$message),
        pid = Sys.getpid()
      )
    }
  )
}


#' Consolida os resultados de uma execução paralela
#'
#' @param res_lista Lista de resultados retornados pelos workers.
#'
#' @return Uma data.table com os resultados válidos consolidados.
consolida_resultados_paralelos <- function(res_lista) {
  stopifnot(is.list(res_lista))

  erros <- Filter(
    function(x) {
      data.table::is.data.table(x) && "erro_worker" %in% names(x)
    },
    res_lista
  )

  if (length(erros) > 0L) {
    erros_dt <- data.table::rbindlist(erros, fill = TRUE)

    warning(sprintf(
      "Ocorreram erros em %d worker(s). Mensagem: %s",
      length(erros),
      paste(unique(erros_dt$erro_worker), collapse = " | ")
    ))
  }

  validos <- Filter(
    function(x) {
      (data.table::is.data.table(x) || is.data.frame(x)) &&
        !("erro_worker" %in% names(x))
    },
    res_lista
  )

  if (length(validos) == 0L) {
    stop("Nenhum resultado foi gerado.")
  }

  dados_out <- data.table::rbindlist(validos, fill = TRUE)

  if (nrow(dados_out) == 0L) {
    stop("Nenhum resultado foi gerado.")
  }

  dados_out
}


#' Executa uma função em paralelo mostrando o progresso de execução
#'
#' @param n_execucoes Número total de execuções.
#' @param funcao Função a ser executada em paralelo retornando um data.table.
#'  Deve aceitar `execucao` (número da execução) como parâmetro.
#' @param lista_argumentos Lista nomeada com argumentos adicionais para `funcao`.
#' @param n_processos Número de workers. Se NULL, a função tentará usar 60% dos núcleos disponíveis, com um mínimo de 2 e máximo de 12.
#' @param seed Controle de semente para future.apply.
#'
#' @return data.table consolidada com os resultados válidos.
#'
#' @examples
#' f <- function(execucao, multiplicador) {
#'   data.table::data.table(
#'     execucao = execucao,
#'     valor = execucao * multiplicador
#'   )
#' }
#' resultado <- execucao_paralela(
#'   n_execucoes = 4,
#'   funcao = f,
#'   lista_argumentos = list(multiplicador = 2),
#' )
#' print(resultado)
execucao_paralela <- function(
  n_execucoes,
  funcao,
  lista_argumentos = list(),
  n_processos = NULL,
  seed = TRUE
) {
  stopifnot(length(n_execucoes) == 1L, n_execucoes >= 1L)
  stopifnot(is.function(funcao))
  stopifnot(is.list(lista_argumentos))

  if (is.null(n_processos)) {
    n_processos <- min(
      12L,
      max(1L, floor(parallelly::availableCores() * 0.6))
    )
  }

  stopifnot(length(n_processos) == 1L, n_processos >= 1L)

  message(sprintf("Workers configurados: %d", n_processos))

  # Salva o plano atual e restaura ao final para evitar efeitos colaterais
  plano_antigo <- future::plan()
  on.exit(future::plan(plano_antigo), add = TRUE)

  future::plan(future::multisession, workers = n_processos)
  if (requireNamespace("progress", quietly = TRUE)) {
    progressr::handlers("progress")
  } else {
    progressr::handlers("txtprogressbar")
  }

  res_lista <- progressr::with_progress({
    p <- progressr::progressor(steps = n_execucoes)

    future.apply::future_lapply(
      X = seq_len(n_execucoes),
      FUN = function(execucao) {
        args <- c(
          list(
            execucao = execucao,
            funcao = funcao
          ),
          lista_argumentos,
          p = p
        )

        do.call(faz_uma_execucao_safe, args)
      },
      future.seed = seed
    )
  })

  consolida_resultados_paralelos(res_lista)
}
