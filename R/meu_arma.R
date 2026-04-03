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

  fit <- tryCatch(
    withCallingHandlers(
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
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
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

cola_series <- function(serie_fase1, serie_fase2, numero_de_novas_observacoes) {
  n_inicial <- length(serie_fase1)

  if (numero_de_novas_observacoes >= n_inicial) {
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
