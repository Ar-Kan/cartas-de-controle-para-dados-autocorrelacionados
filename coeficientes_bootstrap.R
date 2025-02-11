library(BTSR)
library(ggplot2)
library(dplyr)

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


sim <- function(n, coefs, y.start = NULL) {
  # Simula uma série temporal BAR de tamanho `n`
  BARFIMA.sim(
    n = n,
    coefs = coefs,
    y.start = y.start,
    error.scale = 1,
    complete = FALSE
  )
}

fit <- function(yt, start) {
  # Ajusta um modelo BAR
  BARFIMA.fit(
    yt = yt,
    start = start,
    p = 1, # Ordem do polinômio AR
    d = 0,
    error.scale = 1,
    report = FALSE
  )
}

seq_colagens <- c(1, 10, 25, 40, 50, 75, 99)

coeficientes_estimados <- list()
controle <- list()
for (i in seq_colagens) {
  coeficientes_estimados[[i]] <- list()
  controle[[i]] <- list()
}

TAMANHO <- 500
for (k in 1:TAMANHO) {
  amostra_inicial <- sim(n = 100, coefs = coeficientes())
  
  coef_inicial <- fit(
    yt = amostra_inicial,
    start = list(alpha = 0.1, nu = 20, phi = 0.1)
  )
  coef_inicial
  alpha.est <- coef_inicial$coef["alpha"][[1]]
  phi.est <- coef_inicial$coef["phi"][[1]]
  nu.est <- coef_inicial$coef["nu"][[1]]
  
  print(paste0("Bootstrap (", k, ") - phi: ", phi.est))
  bootstrap <- list()
  for (i in 1:1000) {
    amostra <- sim(
      n = 100,
      coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = phi.est)
    )
    coef <- fit(
      yt = amostra,
      start = list(alpha = alpha.est, nu = nu.est, phi = phi.est)
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
  
  # Quartis do bootstrap
  limite.inf <- quantile(
    unlist(lapply(bootstrap, function(x) x$coef$phi)),
    probs = 0.025
  )
  limite.sup <- quantile(
    unlist(lapply(bootstrap, function(x) x$coef$phi)),
    probs = 0.95
  )
  
  # Príxmas amostras
  proximas_amostras_100 <- sim(
    n = 100,
    coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = phi.est + 0.3)
  )
  for (i in seq_colagens) {
    boot.amostra <- sample(seq_len(length(bootstrap)), 1)
    amostra <- bootstrap[[boot.amostra]]$amostra
    proximas_amostras <- proximas_amostras_100[1:i]
    if (i == 100) {
      novo_dataset <- proximas_amostras
    } else {
      novo_dataset <- c(
        amostra[(i + 1):length(amostra)],
        proximas_amostras
      )
    }
  
    coef <- fit(
      yt = novo_dataset,
      start = list(alpha = alpha.est, nu = nu.est, phi = phi.est)
    )
    phi.nova.amostra <- coef$coef["phi"][[1]]
    coeficientes_estimados[[i]] <- c(
      coeficientes_estimados[[i]],
      phi.nova.amostra
    )
    controle[[i]] <- c(
      controle[[i]],
      phi.nova.amostra < limite.inf | phi.nova.amostra > limite.sup
    )
  }
}

df <- data.frame(
  colagem = rep(seq_colagens, each = TAMANHO),
  controle = unlist(controle)
) %>%
  group_by(colagem) %>%
  summarise(
    porcentagem_fora = sum(controle)/n()
  )


ggplot(df, aes(x = colagem, y = porcentagem_fora)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Colagem",
    y = "Porcentagem de amostras fora do intervalo"
  ) +
  theme_minimal()
