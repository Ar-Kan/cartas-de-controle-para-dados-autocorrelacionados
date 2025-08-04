library(ggplot2)
library(dplyr)
library(data.table)
library(MASS)

set.seed(42)

# Parâmetros gerais
PHI_REAL <- 0.2 # valor verdadeiro de Φ
THETA_REAL <- 0.5 # valor verdadeiro de Θ
# N_INICIAL <- 100 # tamanho da amostra inicial
tamanhos_amostrais_iniciais <- c(100) # c(50, 100, 200) # tamanho da amostra inicial
MC <- 300 # número de iterações Monte Carlo
B <- 800 # número de réplicas bootstrap

# DESVIOS_PHI <- c(0, 0.2, 0.6)
DESVIOS_PHI <- c(0, 0.2, 0.6)

# Tamanhos das novas observações em H1
passo <- 13 # Passo para a sequência de colagens
bc <- c(1, 3, 7) # c(1, 3, 5, 7, 9)
SEQ_COLAGENS <- unlist(lapply(seq(0, 100, passo), function(x) bc + x))
SEQ_COLAGENS <- c(SEQ_COLAGENS, 135, 150, 165, 175, 185, 199)

DADOS_OUT <- data.table(
  # Tamanho da amostra inicial
  n0 = numeric(),
  # Quantidade de novas observações
  n1 = numeric(),
  # Desvio aplicado ao Φ
  desvio = numeric(),
  # Indica se a série colada está fora de controle
  fora_de_controle = logical(),
  t2.inf = numeric(),
  t2.sup = numeric(),
  t2.controle = numeric()
)

simula_arma <- function(n, phi, theta, sd = NULL) {
  # Simula ARMA(1,1)
  args <- list(n = n, model = list(ar = phi, ma = theta))
  if (!is.null(sd)) {
    args$sd <- sd
  }
  do.call(arima.sim, args)
}

fit_arma <- function(serie) {
  # Ajusta ARMA(1,1) e retorna o modelo
  tryCatch(
    {
      return(
        arima(
          serie,
          order = c(1, 0, 1), # AR(1) e MA(1)
          include.mean = FALSE,
          method = "ML"
        )
      )
    },
    error = \(e) {
      warning("Erro ao ajustar ARMA(1,1): ", e)
      return(NA)
    }
  )
}

# Testa estacionaridade e invertibilidade
arma_valido <- function(phi, theta) {
  stab <- all(Mod(polyroot(c(1, -phi))) > 1)
  inv <- all(Mod(polyroot(c(1, theta))) > 1)
  stab && inv
}

# Extrai coeficientes e matriz de covariância do modelo ARMA
arma_coef <- function(modelo) {
  coefs <- coef(modelo)
  list(
    p = length(grep("^ar", names(coefs))), # número de termos AR
    q = length(grep("^ma", names(coefs))), # número de termos MA
    coef = coef(modelo), # coeficientes do modelo
    vcov = vcov(modelo), # matriz de covariância dos coeficientes
    vcov_inv = solve(vcov(modelo)) # inversa da matriz de covariância
  )
}

print(paste0("Iniciando simulações com ", MC, " iterações Monte Carlo..."))

for (mc in seq_len(MC)) {
  if (mc %% 10 == 0) {
    print(paste0("Iteração Monte Carlo: ", mc, "/", MC, " (", round(mc / MC * 100, 2), "% completo)"))
  }

  for (N_INICIAL in tamanhos_amostrais_iniciais) {
    ## 1) Simula e ajusta a série inicial
    serie0 <- simula_arma(N_INICIAL, PHI_REAL, THETA_REAL)
    fit0 <- fit_arma(serie0)
    coef0 <- arma_coef(fit0)

    for (desvio in DESVIOS_PHI) {
      # Próximas amostras
      serie1 <- simula_arma(N_INICIAL, PHI_REAL + desvio, THETA_REAL)

      for (sc in SEQ_COLAGENS) {
        if (sc >= N_INICIAL) {
          # Não realiza colagem após exceder o tamanho inicial
          break
        }

        # Série de controle: últimas sc observações de serie0 + primeiras sc de serie1
        serie1.controle <- c(tail(serie0, N_INICIAL - sc), head(serie1, sc))
        stopifnot("A série de controle não possui o tamanho esperado" = length(serie1.controle) == N_INICIAL)

        # Ajusta a série colada
        phi1 <- coef(fit_arma(serie1.controle))[1]
        theta1 <- coef(fit_arma(serie1.controle))[2]
        if (is.na(phi1) || is.na(theta1)) {
          print("Ajuste ARMA falhou, pulando iteração.")
          next
        }

        ## 2) Bootstrap paramétrico com colagem, usando a série0
        phis.boot <- numeric(B)
        thetas.boot <- numeric(B)
        for (b in seq_len(B)) {
          # 2.1) Gera um coef* via mvrnorm, garantindo estacionaridade e invertibilidade
          repeat {
            coef.star <- mvrnorm(1, mu = coef0$coef, Sigma = coef0$vcov)
            # TODO: Assim está correto?
            phi.star <- if (coef0$p > 0) coef.star[1:coef0$p] else numeric(0)
            theta.star <- if (coef0$q > 0) coef.star[(coef0$p + 1):(coef0$p + coef0$q)] else numeric(0)
            if (arma_valido(phi.star, theta.star)) break
          }
          # 2.2) Simula nova amostra ARMA com coef*
          serie.b <- simula_arma(N_INICIAL, phi.star, theta.star)
          # 2.3) Une as últimas sc observações de serie0 com as primeiras sc de serie.b
          serie.colada <- c(tail(serie0, N_INICIAL - sc), head(serie.b, sc))
          stopifnot("A série bootstrap não possui o tamanho esperado" = length(serie.colada) == N_INICIAL)

          # 2.4) Refit ARMA na série colada
          fit.b <- fit_arma(serie.colada)
          if (
            !(inherits(fit.b, "Arima") || inherits(fit.b, "ar")) && is.na(fit.b)
          ) {
            print("Ajuste ARMA falhou no bootstrap, pulando iteração.")
            phis.boot[b] <- NA
            thetas.boot[b] <- NA
            next
          }

          # Verifica se o modelo ARMA ajustado está válido
          stopifnot("O modelo ARMA ajustado não é válido" = !is.null(fit.b) &&
            !is.na(coef(fit.b)[1]) &&
            !is.na(coef(fit.b)[2]))

          # 2.5) Armazena φ* e Θ* estimados
          phis.boot[b] <- coef(fit.b)[1]
          thetas.boot[b] <- coef(fit.b)[2]
        }

        phis.boot <- na.omit(phis.boot)
        thetas.boot <- na.omit(thetas.boot)

        coef.vcov <- matrix(
          c(
            var(phis.boot), cov(phis.boot, thetas.boot),
            cov(thetas.boot, phis.boot), var(thetas.boot)
          ),
          nrow = 2, ncol = 2, byrow = TRUE
        )

        coef.vcov.inv <- solve(coef.vcov)
        phi.m <- mean(phis.boot)
        theta.m <- mean(thetas.boot)

        t2 <- numeric(B)
        for (b in seq_len(length(phis.boot))) {
          diff.b <- c(phis.boot[b], thetas.boot[b]) - c(phi.m, theta.m)
          t2[b] <- t(diff.b) %*% coef.vcov.inv %*% diff.b
        }

        ## 3) Calcula limites de controle 2.5% / 97.5%
        quant.boot <- quantile(t2, probs = c(0.025, 0.975))

        # 4) Calcula T² para a série H1
        t2.controle <- (
          t(c(phi1, theta1) - c(phi.m, theta.m)) %*%
            coef.vcov.inv %*%
            (c(phi1, theta1) - c(phi.m, theta.m))
        )[1]

        DADOS_OUT <- rbind(
          DADOS_OUT,
          data.table(
            n0 = N_INICIAL,
            n1 = sc,
            desvio = desvio,
            # 5) Verifica se a série de H1 está fora de controle
            fora_de_controle = (t2.controle < quant.boot[1] || t2.controle > quant.boot[2]),
            t2.inf = quant.boot[1],
            t2.sup = quant.boot[2],
            t2.controle = t2.controle
          )
        )
      }
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
    x_label = n1 + n0,
  )

# Gráfico de proporção de séries fora de controle
p1 <- ggplot(DADOS_RESUMO, aes(x = x_label, y = proporcao, color = factor(desvio))) +
  geom_line(linewidth = 0.8) +
  # geom_smooth(linewidth = 0.8, se = FALSE, method = "gam", formula = y ~ s(x, k = 5)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  annotate(
    "text",
    x = min(DADOS_RESUMO$x_label) - 5,
    y = 0.05,
    label = "5%",
    vjust = -0.2,
    hjust = 0.5,
    size = 4,
    color = "black"
  ) +
  labs(
    title = "Proporção de Séries Fora de Controle",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Proporção Fora de Controle",
    color = "Desvio de Φ"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(
    ~n0,
    labeller = labeller(
      n0 = \(x) paste0("Tamanho Inicial: ", x)
    ),
    scales = "free_x"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ggsave("controle_online_multiplos_ar.png", p1, width = 10, height = 6, units = "in", dpi = 300)

print(p1)
#
# # Gráfico dos limites de controle
# p2 <- ggplot(DADOS_OUT, aes(x = n1 + n0, y = phi1, color = fora_de_controle)) +
#   geom_errorbar(aes(ymin = limite.inf, ymax = limite.sup), color = "gray", alpha = 0.5, width = 5) +
#   geom_point() +
#   geom_hline(yintercept = PHI_REAL, linetype = "dashed", color = "black") +
#   labs(
#     title = "Limites de Controle para Séries Coladas",
#     x = "Tamanho da Série Colada (n0 + n1)",
#     y = "Φ Estimado",
#     color = "Fora de Controle"
#   ) +
#   scale_color_manual(values = c("TRUE" = "#e20f25", "FALSE" = "#3b7fb6")) +
#   theme_minimal() +
#   facet_wrap(
#     ~desvio,
#     labeller = labeller(
#       desvio = \(x) ifelse(
#         x == 0,
#         "Sem Desvio",
#         ifelse(
#           x > 0,
#           paste0("Desvio +", x),
#           paste0("Desvio -", x)
#         )
#       )
#     )
#   ) +
#   theme(legend.position = "bottom")
#
# # ggsave("controle_limites_ar.png", p2, width = 10, height = 6, units = "in", dpi = 300)
#
# print(p2)
