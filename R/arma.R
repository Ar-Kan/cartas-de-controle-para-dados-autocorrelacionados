# Simula ARMA(1,1) com burn-in de 200
simula_arma <- function(n, phi, theta) {
  obs <- arima.sim(
    n = 200 + n,
    model = list(ar = phi, ma = theta)
  )
  # Retorna as últimas n observações
  tail(obs, n)
}

# Ajusta ARMA(1,1)
#
# Usa uma parametrização alternativa para garantir estacionaridade dos
# termos AR, e trata a invertibilidade dos MA depois da otimização (invertibility enforcement).
# NOTA: Não funciona com o método `CSS` puro
#
# `CSS-ML`: usa CSS para dar um chute inicial e refina com verossimilhança
fit_arma <- function(serie, phi = NULL, theta = NULL, transform.pars = FALSE) {
  avisos <- character()

  fit <- withCallingHandlers(
    arima(
      serie,
      order = c(1, 0, 1),
      include.mean = FALSE,
      transform.pars = transform.pars,
      method = "CSS-ML",
      init = c(ar1 = phi, ma1 = theta),
      optim.method = "BFGS",
      optim.control = list(maxit = 1000, reltol = 1e-10)
    ),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ) |> tryNull()

  if (
    is.null(fit) ||
      coef(fit) |> tryNull() |> is.null() ||
      vcov(fit) |> tryNull() |> is.null()
  ) {
    return(list(
      fit = NULL,
      warnings = avisos,
      convergiu = FALSE,
      erro = TRUE
    ))
  }

  list(
    fit = fit,
    warnings = avisos,
    convergiu = length(avisos) == 0,
    erro = FALSE
  )
}

# Testa estacionaridade (para φ)
ar_valido <- function(phi) {
  is.finite(phi) && all(Mod(polyroot(c(1, -phi))) > 1)
}

# Testa invertibilidade (para θ)
ma_valido <- function(theta) {
  is.finite(theta) && all(Mod(polyroot(c(1, theta))) > 1)
}

# Colagem para série bootstrap: mantém os últimos de Fase I e os primeiros de Fase II*
cola_series <- function(serie_fase1, serie_fase2, numero_de_novas_observacoes) {
  n_inicial <- length(serie_fase1)

  if (numero_de_novas_observacoes >= n_inicial) {
    # Sem colagem
    out <- serie_fase2[seq_len(n_inicial)]
  } else {
    out <- c(
      tail(serie_fase1, n_inicial - numero_de_novas_observacoes),
      head(serie_fase2, numero_de_novas_observacoes)
    )
  }

  stopifnot("Série resultante tem comprimento diferente do esperado" = length(out) == n_inicial)
  out
}


# Gera coeficientes bootstrap válidos para ARMA(1,1) com dist N(coef, matriz_vcov).
# Para séries de ordens maiores devemos usarmor PACF como o `stats::arima()` faz.
#
# A operação é feita em 4 etapas:
# 1. Transformação para o espaço onde os parâmetros são irrestritos (usando atanh para mapear (-1, 1) para ℝ).
# 2. Ajuste da covariância no espaço transformado usando o método delta.
# 3. Sorteio de um novo vetor de parâmetros no espaço transformado.
# 4. Transformação de volta para o espaço original usando tanh.
#
# Método delta para aproximar a covariância no espaço transformado.
# Seja J a jacobiana de g avalidada em coef, temos que se:
#   u = g(coef)
# então:
#   Var(u) ~= J Var(coef) J'
# Onde:
#   g(phi) = atanh(phi) => d/dphi atanh(phi) = 1 / (1 - phi^2)
#   g(theta) = atanh(theta) => d/dtheta atanh(theta) = 1 / (1 - theta^2)
#
# Como a transformação é separada componente a componente, a jacobiana é diagonal.
bootstrap_coef_validos <- function(coef, matriz_vcov, eps = 1e-6) {
  # Validações
  if (!is.matrix(matriz_vcov) || any(!is.finite(matriz_vcov))) {
    return(NULL)
  }

  phi0 <- unname(coef["ar1"])
  theta0 <- unname(coef["ma1"])

  if (!is.finite(phi0) || !is.finite(theta0)) {
    return(NULL)
  }

  # proteção contra valores exatamente na fronteira
  phi0 <- max(min(phi0, 1 - eps), -1 + eps)
  theta0 <- max(min(theta0, 1 - eps), -1 + eps)

  # Média no espaço transformado
  mu_u <- c(
    ar1 = atanh(phi0),
    ma1 = atanh(theta0)
  )

  # Jacobiana da transformação g no ponto coef
  J <- diag(c(
    1 / (1 - phi0^2),
    1 / (1 - theta0^2)
  ))

  Sigma_u <- J %*% matriz_vcov %*% t(J)
  Sigma_u <- as.matrix(Sigma_u)

  rownames(Sigma_u) <- colnames(Sigma_u) <- c("ar1", "ma1")

  # Valida matriz
  if (any(!is.finite(Sigma_u))) {
    return(NULL)
  }

  # Verifica se a matriz é positiva-definida para evitar erros no `mvrnorm`
  if (is.null(chol(Sigma_u) |> tryNull())) {
    return(NULL)
  }

  # Sorteia valores no espaço transformado
  # U* ~ N(mu_u, Sigma_u)
  u_star <- MASS::mvrnorm(1, mu = mu_u, Sigma = Sigma_u)

  if (is.null(u_star) || any(!is.finite(u_star))) {
    return(NULL)
  }

  # Volta ao espaço original
  #   phi* = tanh(u_phi*)
  #   theta* = tanh(u_theta*)
  phi_star <- tanh(unname(u_star["ar1"]))
  theta_star <- tanh(unname(u_star["ma1"]))

  # Proteção final
  if (!is.finite(phi_star) || !is.finite(theta_star)) {
    return(NULL)
  }

  c(phi_star = phi_star, theta_star = theta_star)
}
