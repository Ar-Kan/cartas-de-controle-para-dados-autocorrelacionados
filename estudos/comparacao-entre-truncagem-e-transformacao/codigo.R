# Usa uma parametrização para garantir a validade de φ e θ amostrados durante o bootstrap.
# Também define `transform.pars = TRUE` no ajuste do modelo por meio de `stats::arima()``.
USAR_TRANSFORMACAO_NO_ESPACO_PARAMETRICO <- FALSE


# Ajusta ARMA(1,1)
fit_arma <- function(serie, phi = NULL, theta = NULL) {
  aviso <- NULL

  fit <- tryCatch(
    withCallingHandlers(
      arima(
        serie,
        order = c(1, 0, 1),
        include.mean = FALSE,
        # Usa uma parametrização alternativa para garantir estacionaridade dos
        # termos AR, e trata a invertibilidade dos MA depois da otimização (invertibility enforcement).
        # NOTA: Não funciona com o método `CSS` puro
        transform.pars = USAR_TRANSFORMACAO_NO_ESPACO_PARAMETRICO,
        # `CSS-ML`: usa CSS para dar um chute inicial e refina com verossimilhança
        method = "CSS-ML",
        init = c(ar1 = phi, ma1 = theta),
        optim.method = "BFGS",
        optim.control = list(maxit = 1000, reltol = 1e-10)
      ),
      warning = function(w) {
        aviso <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  list(
    fit = fit,
    warning = aviso,
    convergiu = is.null(aviso)
  )
}

# Gera coeficientes bootstrap válidos para ARMA(1,1) com dist N(coef0$coef, coef0$vcov)
# Para séries de ordens maiores devemos usarmor PACF como o `stats::arima()` faz
bootstrap_coef_validos <- function(coef0, eps = 1e-6) {
  Sigma <- coef0$vcov
  beta0 <- coef0$coef

  # Validações
  if (!is.matrix(Sigma) || any(!is.finite(Sigma))) {
    return(NULL)
  }

  phi0 <- unname(beta0["ar1"])
  theta0 <- unname(beta0["ma1"])

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

  # Método delta para aproximar a covariância no espaço transformado.
  #
  # Seja J a jacobiana de g avalidada em beta0, temos que se:
  #   u = g(beta)
  # então:
  #   Var(u) ≈ J Var(beta) J'
  # Onde:
  #   g(phi) = atanh(phi) => d/dphi atanh(phi) = 1 / (1 - phi^2)
  #   g(theta) = atanh(theta) => d/dtheta atanh(theta) = 1 / (1 - theta^2)
  #
  # Como a transformação é separada componente a componente,
  # a jacobiana é diagonal.
  J <- diag(c(
    1 / (1 - phi0^2),
    1 / (1 - theta0^2)
  ))

  Sigma_u <- J %*% Sigma %*% J
  Sigma_u <- as.matrix(Sigma_u)

  rownames(Sigma_u) <- colnames(Sigma_u) <- c("ar1", "ma1")

  # Valida matriz
  if (any(!is.finite(Sigma_u))) {
    return(NULL)
  }

  if (inherits(try(chol(Sigma_u), silent = TRUE), "try-error")) {
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


# Gera coeficientes e adiciona variabilidade
if (USAR_TRANSFORMACAO_NO_ESPACO_PARAMETRICO) {
  # Gera coeficientes bootstrap válidos para ARMA(1,1) com dist N(coef0$coef, coef0$vcov)
  coef_star <- bootstrap_coef_validos(coef0)
  if (is.null(coef_star)) next

  phi_star <- unname(coef_star["phi_star"])
  theta_star <- unname(coef_star["theta_star"])
} else {
  # Gera coeficientes bootstrap sem restrição, o que pode levar a valores inválidos de φ* e θ*.
  coef_star <- tryCatch(
    MASS::mvrnorm(1, mu = coef0$coef, Sigma = coef0$vcov),
    error = function(e) NULL
  )
  if (is.null(coef_star) || any(!is.finite(coef_star))) next

  phi_star <- unname(coef_star["ar1"])
  theta_star <- unname(coef_star["ma1"])

  # Valida os coeficientes bootstrap gerados
  if (!ar_valido(phi_star) || !ma_valido(theta_star)) next
}
