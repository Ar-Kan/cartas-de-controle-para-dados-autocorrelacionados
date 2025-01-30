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

amostra_inicial <- sim(n = 100, coefs = coeficientes())

coef_inicial <- fit(
  yt = amostra_inicial,
  start = list(alpha = 0.1, nu = 20, phi = 0.1)
)
coef_inicial
alpha.est <- coef_inicial$coef["alpha"][[1]]
phi.est <- coef_inicial$coef["phi"][[1]]
nu.est <- coef_inicial$coef["nu"][[1]]

print("Bootstrap")
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


seq_colagens <- c(1, seq(25, 75, 25), 99)

bootstrap_colagem <- list()
for (i in seq_colagens) {
  print(paste0("Colagem com ", i))
  amostras_coladas <- list()
  for (j in seq(1, 1000, 2)) {
    boot.amostra.1 <- j # sample(seq_len(length(bootstrap)), 1)
    boot.amostra.2 <- j + 1 # sample(seq_len(length(bootstrap)), 1)
    amostra.1 <- bootstrap[[boot.amostra.1]]$amostra
    amostra.2 <- bootstrap[[boot.amostra.2]]$amostra
    amostra_colada <- c(amostra.1[1:i], amostra.2[(i + 1):100])
    coef <- fit(
      yt = amostra_colada,
      start = list(alpha = alpha.est, nu = nu.est, phi = phi.est)
    )
    
    amostras_coladas[[j]] <- list(
      amostra = amostra_colada,
      coef = list(
        alpha = coef$coef["alpha"][[1]],
        phi = coef$coef["phi"][[1]],
        nu = coef$coef["nu"][[1]]
      )
    )
  }
  bootstrap_colagem[[i]] <- amostras_coladas
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


coeficientes_estimados <- list()
for (i in seq_colagens) {
  coeficientes_estimados[[i]] <- list()
}

# Monte Carlo
print("Monte Carlo")
for (k in 1:500) {
  proximas_amostras_100 <- sim(
    n = 100,
    coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = phi.est)
  )
  for (i in seq_colagens) {
    boot.amostra <- sample(seq_len(999), 1)
    amostra <- bootstrap_colagem[[i]][[boot.amostra]]$amostra
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
    # coeficientes_estimados[[i]] <- c(
    #   coeficientes_estimados[[i]],
    #   list(
    #     boot.amostra = boot.amostra,
    #     boot = bootstrap[[boot.amostra]],
    #     amostra = novo_dataset,
    #     coef = list(
    #       alpha = coef$coef["alpha"][[1]],
    #       phi = coef$coef["phi"][[1]],
    #       nu = coef$coef["nu"][[1]]
    #     )
    #   )
    # )
    coeficientes_estimados[[i]] <- c(
      coeficientes_estimados[[i]],
      coef$coef["phi"][[1]]
    )
  }
}

df <- data.frame(
  colagem = rep(seq_colagens, each = 50),
  phi = unlist(coeficientes_estimados)
)

library(ggplot2)

ggplot(df, aes(x = colagem, y = phi, group = colagem)) +
  geom_boxplot() +
  geom_hline(yintercept = c(limite.inf, limite.sup), linetype = "dashed") +
  labs(
    title = "Estimativa de phi",
    x = "Número de colagens",
    y = "Estimativa de phi"
  ) +
  facet_wrap(~colagem, scales = "free_y")
