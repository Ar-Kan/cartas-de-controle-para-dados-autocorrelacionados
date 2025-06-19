library(ggplot2)
library(dplyr)
library(data.table)

set.seed(42)

# Parâmetros gerais
PHI_REAL <- 0.2  # valor verdadeiro de Φ
N_INICIAL <- 100  # tamanho da amostra inicial
MC <- 300  # número de iterações Monte Carlo
B <- 800  # número de réplicas bootstrap
CONTROLE_OUT <- matrix(NA, nrow = MC, ncol = 2)  # armazena limites inferior/superior

DESVIOS_PHI <- c(0, 0.2, 0.4)

# Tamanhos das novas observações em H1
passo <- 20 # 20 # Passo para a sequência de colagens
bc <- c(1, 5) # c(1, 3, 5, 7, 9)
SEQ_COLAGENS <- unlist(lapply(seq(0, 100, passo), function(x) bc + x))
SEQ_COLAGENS <- c(SEQ_COLAGENS, 135, 150, 165, 175, 185, 199)

DADOS_OUT <- data.table(
  # Tamanho da amostra inicial
  n0 = numeric(),
  # φ estimado da série inicial
  phi0 = numeric(),
  # Quantidade de novas observações
  n1 = numeric(),
  # φ estimado da série colada
  phi1 = numeric(),
  # Desvio aplicado ao Φ
  desvio = numeric(),
  # Indica se a série colada está fora de controle
  fora_de_controle = logical(),
  limite.inf = numeric(),
  limite.sup = numeric()
)

simula_ar1 <- function(n, phi, sd = NULL, y.start = NULL) {
  # Simula AR(1) com phi e desvio sd (se fornecido), opcionalmente definindo y.start
  args <- list(n = n, model = list(ar = phi))
  if (!is.null(sd))
    args$sd <- sd
  if (!is.null(y.start))
    args$start.innov <- y.start
  do.call(arima.sim, args)
}

fit_ar1 <- function(serie) {
  # Ajusta AR(1) e retorna o modelo
  arima(
    serie,
    order = c(1, 0, 0),
    include.mean = FALSE,
    method = "ML"
  )
}

for (mc in seq_len(MC)) {
  if (mc %% 10 == 0) {
    print(paste0("Iteração Monte Carlo: ", mc, "/", MC, " (", round(mc / MC * 100, 2), "% completo)"))
  }

  ## 1) Simula e ajusta a série inicial
  serie0 <- simula_ar1(N_INICIAL, PHI_REAL)
  fit0 <- fit_ar1(serie0)
  phi0 <- coef(fit0)[1]
  # erro padrão assintótico de phi em AR(1): sqrt(σ² (1-φ²)/n)
  sd_phi0 <- sqrt((fit0$sigma2 * (1 - phi0^2)) / N_INICIAL)
  
  for (desvio in DESVIOS_PHI) {
    # Próximas amostras
    serie1 <- simula_ar1(N_INICIAL, PHI_REAL + desvio)
    
    for (sc in SEQ_COLAGENS) {
      if (sc >= N_INICIAL) {
        # Não realiza colagem após exceder o tamanho inicial
        break
      }

      serie1.controle <- c(tail(serie0, N_INICIAL - sc), head(serie1, sc))
      stopifnot("A série de controle não possui o tamanho esperado" = length(serie1.controle) == N_INICIAL)
      
      # Ajusta a série colada
      phi1 <- coef(fit_ar1(serie1))[1]
  
      ## 2) Bootstrap paramétrico com colagem: reutiliza série0
      phis.boot <- numeric(B)
      for (b in seq_len(B)) {
        # 2.1) gera um φ* via rnorm e força |φ*|<1 para estacionaridade
        repeat {
          phi.star <- rnorm(1, mean = phi0, sd = sd_phi0)
          if (abs(phi.star) < 1)
            break
        }
        # 2.2) simula nova amostra AR(1) com φ*
        serie.b <- simula_ar1(N_INICIAL, phi.star)
        # 2.3) cola: puxa as primeiras k observações de serie.b para frente
        serie.colada <- c(tail(serie0, N_INICIAL - sc), head(serie.b, sc))
        stopifnot("A série bootstrap não possui o tamanho esperado" = length(serie.colada) == N_INICIAL)
        # 2.4) refit
        fit.b <- fit_ar1(serie.colada)
        # 2.5) armazena φ* estimado
        phis.boot[b] <- coef(fit.b)[1]
      }
  
      ## 3) calcula limites de controle 2.5% / 97.5%
      quant.boot <- quantile(phis.boot, probs = c(0.025, 0.975))
      
      DADOS_OUT <- rbind(
        DADOS_OUT,
        data.table(
          n0 = N_INICIAL,
          phi0 = phi0,
          n1 = sc,
          phi1 = phi1,
          desvio = desvio,
          fora_de_controle = (phi1 < quant.boot[1] || phi1 > quant.boot[2]),
          limite.inf = quant.boot[1],
          limite.sup = quant.boot[2]
        )
      )
      
    }
  }
}

# Resumo dos dados: proporção de séries fora de controle
DADOS_RESUMO <- DADOS_OUT %>%
  group_by(
    n0,
    n1,
    desvio,
  ) %>%
  summarise(
    proporcao = mean(fora_de_controle)
  ) %>%
  mutate(
    dp = sqrt((proporcao * (1 - proporcao)) / MC),
    x_label = n1 + n0,
  )

# Gráfico de proporção de séries fora de controle
p1 <- ggplot(DADOS_RESUMO, aes(x = x_label, y = proporcao, color = factor(desvio))) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  # geom_point() +
  geom_errorbar(aes(ymin = proporcao - dp, ymax = proporcao + dp)) +
  labs(
    title = "Proporção de Séries Fora de Controle",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Proporção Fora de Controle",
    color = "Desvio de Φ"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# ggsave("controle_online_ar.png", p1, width = 10, height = 6, units = "in", dpi = 300)

print(p1)

# Gráfico dos limites de controle
p2 <- ggplot(DADOS_OUT, aes(x = n1 + n0, y = phi1, color = fora_de_controle)) +
  geom_errorbar(aes(ymin = limite.inf, ymax = limite.sup), color = "gray", alpha = 0.5) +
  geom_point() +
  geom_hline(yintercept = PHI_REAL, linetype = "dashed", color = "black") +
  labs(
    title = "Limites de Controle para Séries Coladas",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Φ Estimado",
    color = "Fora de Controle"
  ) +
  scale_color_manual(values = c("TRUE" = "#e20f25", "FALSE" = "#3b7fb6")) +
  theme_minimal() +
  facet_wrap(~ desvio) +
  theme(legend.position = "bottom")

# ggsave("controle_limites_ar.png", p2, width = 10, height = 6, units = "in", dpi = 300)

print(p2)
