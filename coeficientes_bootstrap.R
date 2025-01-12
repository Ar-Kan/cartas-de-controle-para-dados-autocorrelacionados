library(BTSR)

coeficientes <- function(alpha = 0, nu = 20, phi = 0.2) {
  # βARMA: o modelo de Rocha e Cribari-Neto (2009, 2017)
  #        é obtido definindo `coefs$d = 0`
  #        e `d = FALSE` e `error.scale = 1` (escala preditiva)
  #
  # d: parâmetro de longa dependência (long memory), quando
  #    diferente de zero, produz um modelo ARFIMA.
  #
  # ν (nu): parâmetro de precisão, quanto maior o seu valor,
  #         menor a variância condicional (para μ_t fixo).
  #         O pacote BTSR define como padrão ν = 20.
  list(
    alpha = alpha,
    nu = nu,
    phi = phi,
    d = 0
  )
}

amostra_inicial <- BARFIMA.sim(
  n = 100,
  coefs = coeficientes(),
  error.scale = 1,
  complete = FALSE
)

coef_inicial <- BARFIMA.fit(
  yt = amostra_inicial,
  start = list(alpha = 0.1, nu = 20, phi = 0.1),
  p = 1, # Ordem do polinômio AR
  d = 0,
  error.scale = 1,
  report = FALSE
)
coef_inicial
alpha.est <- coef_inicial$coef["alpha"][[1]]
phi.est <- coef_inicial$coef["phi"][[1]]
nu.est <- coef_inicial$coef["nu"][[1]]

bootstrap <- list()
for (i in 1:1000) {
  amostra <- BARFIMA.sim(
    n = 100,
    coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = phi.est),
    error.scale = 1,
    complete = FALSE
  )
  coef <- BARFIMA.fit(
    yt = amostra,
    start = list(alpha = alpha.est, nu = nu.est, phi = phi.est),
    p = 1, # Ordem do polinômio AR
    d = 0,
    error.scale = 1,
    report = FALSE
  )
  bootstrap[[i]] <- list(
    amostra = amostra,
    coef = list(
      alpha = coef$coef["alpha"][[1]],
      phi = coef$coef["phi"][[1]],
      nu = coef$coef["nu"][[1]]
    )
  )
}

# bootstrap

coeficientes_estimados <- list()
for (i in 1:100) {
  boot.amostra <- sample(seq_len(length(bootstrap)), 1)
  amostra <- bootstrap[[boot.amostra]]$amostra
  estimativa <- bootstrap[[boot.amostra]]$coef
  proximas_amostras <- BARFIMA.sim(
    n = i,
    y.start = amostra[length(amostra)],
    coefs = coeficientes(alpha = estimativa$alpha, nu = estimativa$nu, phi = estimativa$phi),
    error.scale = 1,
    complete = FALSE
  )
  if (i == 100) {
    novo_dataset <- proximas_amostras
  } else {
    novo_dataset <- c(
      amostra[(i + 1):length(amostra)],
      proximas_amostras
    )
  }

  coef <- BARFIMA.fit(
    yt = novo_dataset,
    start = list(alpha = alpha.est, nu = nu.est, phi = phi.est),
    p = 1, # Ordem do polinômio AR
    d = 0,
    error.scale = 1,
    report = FALSE
  )
  coeficientes_estimados[[i]] <- list(
    boot.amostra = boot.amostra,
    boot = bootstrap[[n_amostra]],
    amostra = novo_dataset,
    coef = list(
      alpha = coef$coef["alpha"][[1]],
      phi = coef$coef["phi"][[1]],
      nu = coef$coef["nu"][[1]]
    )
  )
}

df <- data.frame(
  alpha = sapply(coeficientes_estimados, function(x) x$coef$alpha),
  phi = sapply(coeficientes_estimados, function(x) x$coef$phi),
  nu = sapply(coeficientes_estimados, function(x) x$coef$nu),
  boot.amostra = sapply(coeficientes_estimados, function(x) x$boot.amostra),
  boot.alpha = sapply(coeficientes_estimados, function(x) x$boot$coef$alpha),
  boot.phi = sapply(coeficientes_estimados, function(x) x$boot$coef$phi),
  boot.nu = sapply(coeficientes_estimados, function(x) x$boot$coef$nu)
)
