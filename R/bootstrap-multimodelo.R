stationary.arima <- function(n, model, sd = 1, max_iter = 1000) {
  erro <- TRUE
  iter <- 0
  while (erro == TRUE && iter < max_iter) {
    iter <- iter + 1
    tryCatch(
    {
      serie <- arima.sim(n = 200 + n, model = model, sd = sd)

      # Verificar se a série é estacionária e invertível
      arima(serie, include.mean = F, order = c(length(model$ar), 0, length(model$ma)))

      erro <- FALSE
    },
      error = function(e) {
        # Não faz nada, apenas continua o loop
      }
    )
  }
  if (iter >= max_iter) {
    stop(paste("Erro. Valores dos coeficientes = ", model, ". Valor de n = ", n))
  }

  tail(serie, n)
}

coeficientes_bootstrap <- function(modelo) {
  if (length(modelo$coef) == 1) {
    coef_bootstrap <- rnorm(1, modelo$coef, sqrt(modelo$var.coef)) |> tryNull()
  } else {
    # Gerar novos coeficientes com base na matriz de variância-covâriancia da estimativa.
    coef_bootstrap <- MASS::mvrnorm(1, modelo$coef, modelo$var.coef) |> tryNull()
  }
  if (is.null(coef_bootstrap) || any(!is.finite(coef_bootstrap))) return(NULL)
  coef_bootstrap
}

multimodel_Q <- function(serie, resultados_mmi, nlags = 10) {

  modelos_validos <- which(resultados_mmi$pesos > 0)

  pesos <- resultados_mmi$pesos[modelos_validos]
  Qs <- numeric(length(pesos))

  for (m in modelos_validos) {

    tryCatch({
      extract <- arima(
        serie,
        include.mean = FALSE,
        order = forecast::arimaorder(resultados_mmi$modelos[[m]]),
        fixed = c(resultados_mmi$modelos[[m]]$coef)
      )
    }, error = function(e) {
      print(serie)
    })

    box_test_i <- Box.test(extract$residuals, lag = nlags, type = c("Ljung-Box"), fitdf = 1)
    Qi <- box_test_i$statistic

    i <- which(modelos_validos == m)
    Qs[i] <- Qi

  }

  Q_mm <- sum(Qs * pesos)
  return(Q_mm)
}

bootstrap_multimodelo <- function(modelos_candidatos, numero_de_boots, numero_de_novas_observacoes) {
  phi_boot <- rep(NA_real_, numero_de_boots)
  theta_boot <- rep(NA_real_, numero_de_boots)

  serie_original <- modelos_candidatos$serie
  n <- length(serie_original)

  for (i in 1:numero_de_boots) {
    # 1. Selecionar modelo baseado nos pesos de Akaike
    modelo <- amostrar_modelo(modelos_candidatos$modelos)

    # 2. Gerar série bootstrap do modelo selecionado
    iter <- 0
    while (iter < 1000) {
      iter <- iter + 1
      coef_bootstrap <- coeficientes_bootstrap(modelo)
      if (is.null(coef_bootstrap)) next

      # Gerar série temporal com coeficientes bootstrapados
      serie_bootstrap <- stationary.arima(
        n = n,
        model = list(
          ar = coef_bootstrap[grepl("^ar", names(coef_bootstrap))],
          ma = coef_bootstrap[grepl("^ma", names(coef_bootstrap))]
        ),
        sd = sqrt(modelo$sigma2)
      ) |> tryNull()

      if (is.null(serie_bootstrap)) next
    }
    if (iter >= 1000) {
      stop(paste("Erro. Valores dos coeficientes = ", modelo$coef, ". Ordem = ", modelo, modelos_validos, modelo_idx))
    }

    # # Q multi modelo da série bootstrap
    # if (lag < n) {
    #   colagem = c(h0[(lag + 1):n], serie_bootstrap[(n - lag + 1):n])
    #
    #   q_bootstrap[i] <- multimodel_Q(colagem, resultados_mmi)
    # }
    # else {
    #   q_bootstrap[i] <- multimodel_Q(serie_bootstrap, resultados_mmi)
    # }

    serie_colada_bootstrap <- cola_series(
      serie_fase1 = serie_original,
      serie_fase2 = serie_bootstrap,
      numero_de_novas_observacoes = numero_de_novas_observacoes
    )
  }

  return(q_bootstrap)
}
