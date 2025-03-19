library(BTSR)
library(ggplot2)
library(ggformula)
library(dplyr)
library(data.table)

BETA_ARMA <- FALSE

coeficientes <- function(alpha = 0, nu = 20, phi = 0.2) {
  if (!BETA_ARMA) {
    # Modelo AR
    return(
      list(ar = c(phi))
    )
  }
  
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
  if (!BETA_ARMA) {
    # Simula uma série temporal AR de tamanho `n`
    return(
      arima.sim(
        n = n,
        model = coefs
      )
    )
  }
  
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
  if (!BETA_ARMA) {
    # Ajusta um modelo AR
    return(
      arima(
        yt,
        order = c(1, 0, 0),
        include.mean = FALSE,
      )
    )
  }
  
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

seq_colagens <- c(1, 24, 49, 75, 99)
phis <- c(-0.1, 0.1, 0.2, 0.3)
tamanho_amostras_iniciais <- c(25, 50, 100, 200)

df <- data.table(
  tamanho_inicial = numeric(),
  tamanho_colagem = numeric(),
  # Phi estimado para cada tamanho de colagem
  phi_estimado = numeric(),
  novo_phi = numeric(),
  # SEQ dos phis para cada tamanho de colagem
  controle = logical()
)

TAMANHO <- 500
for (k in 1:TAMANHO) {
  for (n_inicial in tamanho_amostras_iniciais) {
    amostra_inicial <- sim(n = n_inicial, coefs = coeficientes())
    
    coef_inicial <- fit(
      yt = amostra_inicial,
      start = list(alpha = 0.1, nu = 20, phi = 0.1)
    )
    coef_inicial
    alpha.est <- ifelse(BETA_ARMA, coef_inicial$coef["alpha"][[1]], NA)
    nu.est <- ifelse(BETA_ARMA, coef_inicial$coef["nu"][[1]], NA)
    phi.est <- ifelse(BETA_ARMA, coef_inicial$coef["phi"][[1]], coef_inicial$coef[[1]])
    
    if (BETA_ARMA & phi.est > 0.6) {
      next
    }
    
    print(paste0("Bootstrap (", k, ") - phi: ", phi.est, ", n: ", n_inicial))
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
          alpha = ifelse(BETA_ARMA, coef$coef["alpha"][[1]], NA),
          nu = ifelse(BETA_ARMA, coef$coef["nu"][[1]], NA),
          phi = ifelse(BETA_ARMA, coef$coef["phi"][[1]], coef$coef[[1]])
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
    
    
    for (novo_phi in phis) {
      # Próximas amostras
      proximas_amostras_100 <- sim(
        n = n_inicial,
        coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = novo_phi)
      )
      for (i in seq_colagens) {
        proximas_amostras <- proximas_amostras_100[1:i]
        if (i > n_inicial) {
          next
        } else if (i == n_inicial) {
          novo_dataset <- proximas_amostras
        } else {
          novo_dataset <- c(
            amostra_inicial[(i + 1):length(amostra_inicial)],
            proximas_amostras
          )
        }
        
        ok <- FALSE
        tryCatch({
          coef <- fit(
            yt = novo_dataset,
            start = list(alpha = alpha.est, nu = nu.est, phi = phi.est)
          )
          ok <- TRUE
        }, error = function(e) {
        })
        
        if (!ok) next
        
        phi.nova.amostra <- ifelse(BETA_ARMA, coef$coef["phi"][[1]], coef$coef[[1]])
        
        df <- rbind(
          df,
          data.table(
            tamanho_inicial = n_inicial,
            tamanho_colagem = i,
            phi_estimado = phi.nova.amostra,
            novo_phi = novo_phi,
            controle = phi.nova.amostra < limite.inf | phi.nova.amostra > limite.sup
          )
        )
      }
    }
  }
}

df_resumo <- df %>%
  group_by(tamanho_inicial, tamanho_colagem, novo_phi) %>%
  summarise(
    controle = sum(controle),
    total = n()
  ) %>%
  mutate(
    proporcao = controle / total
  )

p <- ggplot(df_resumo, aes(x = tamanho_colagem, y = proporcao, color = as.factor(novo_phi))) +
  geom_point() +
  geom_line() +
  labs(
    x = "Tamanho da colagem",
    y = "Proporção de controle",
    color = "Novo valor de Phi"
  ) +
  theme_minimal() +
  facet_wrap(~tamanho_inicial)

ggsave("controle_ar.png", p, width = 10, height = 6, units = "in", dpi = 300)

p
