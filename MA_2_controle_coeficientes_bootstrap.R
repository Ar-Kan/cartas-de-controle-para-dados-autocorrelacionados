library(ggplot2)
library(dplyr)
library(data.table)
library(MASS)

set.seed(42)

# ------------------------------
# Parâmetros gerais
# ------------------------------
THETA_REAL <- c(0.5, -0.4) # valores verdadeiros de Θ1 e Θ2
tamanhos_amostrais_iniciais <- c(100)
MC <- 200
B <- 200

DESVIOS_THETA <- c(0, 0.2, 0.6)

passo <- 30
bc <- c(1)
SEQ_COLAGENS <- unlist(lapply(seq(0, 100, passo), function(x) bc + x))
SEQ_COLAGENS <- c(SEQ_COLAGENS, 135, 150, 165, 175, 185, 199)

DADOS_OUT <- data.table(
  n0 = numeric(),
  n1 = numeric(),
  desvio = numeric(),
  fora_de_controle = logical(),
  t2.inf = numeric(),
  t2.sup = numeric(),
  t2.controle = numeric()
)

# ------------------------------
# Funções auxiliares
# ------------------------------

# Simula MA(2)
simula_ma <- function(n, theta, sd = NULL) {
  args <- list(n = 200 + n, model = list(ma = theta))
  if (!is.null(sd)) args$sd <- sd
  obs <- do.call(arima.sim, args)
  return(tail(obs, n))
}

# Ajusta MA(2)
fit_ma <- function(serie) {
  tryCatch(
    {
      return(
        arima(
          serie,
          order = c(0, 0, 2), # MA(2)
          include.mean = FALSE,
          method = "ML"
        )
      )
    },
    error = \(e) {
      warning("Erro ao ajustar MA(2): ", e)
      return(NA)
    }
  )
}

# Testa invertibilidade
ma_valido <- function(theta) {
  inv <- all(Mod(polyroot(c(1, theta))) > 1)
  inv
}

# Extrai coeficientes e matriz de covariância
ma_coef <- function(modelo) {
  coefs <- coef(modelo)
  list(
    q = length(grep("^ma", names(coefs))),
    coef = coefs,
    vcov = vcov(modelo),
    vcov_inv = solve(vcov(modelo))
  )
}

# ------------------------------
# Loop Monte Carlo
# ------------------------------

print(paste0("Iniciando simulações MA(2) com ", MC, " iterações Monte Carlo..."))

for (mc in seq_len(MC)) {
  if (mc %% 10 == 0) {
    print(paste0("Iteração Monte Carlo: ", mc, "/", MC, " (", round(mc / MC * 100, 2), "% completo)"))
  }

  for (N_INICIAL in tamanhos_amostrais_iniciais) {
    ## 1) Simula e ajusta a série inicial
    serie0 <- simula_ma(N_INICIAL, THETA_REAL)
    fit0 <- fit_ma(serie0)
    coef0 <- ma_coef(fit0)

    for (desvio in DESVIOS_THETA) {
      # Próximas amostras (aplica desvio apenas em Θ1)
      serie1 <- simula_ma(N_INICIAL, THETA_REAL + c(desvio, 0))

      for (sc in SEQ_COLAGENS) {
        if (sc >= N_INICIAL) break

        # Série colada (controle)
        serie1.controle <- c(tail(serie0, N_INICIAL - sc), head(serie1, sc))
        stopifnot(length(serie1.controle) == N_INICIAL)

        # Ajusta série colada
        fit_tmp <- fit_ma(serie1.controle)
        theta1 <- coef(fit_tmp)[1]
        theta2 <- coef(fit_tmp)[2]
        if (is.na(theta1) || is.na(theta2)) {
          print("Ajuste MA(2) falhou, pulando iteração.")
          next
        }

        ## 2) Bootstrap paramétrico
        theta1.boot <- numeric(B)
        theta2.boot <- numeric(B)
        for (b in seq_len(B)) {
          repeat {
            coef.star <- mvrnorm(1, mu = coef0$coef, Sigma = coef0$vcov)
            theta.star <- coef.star[1:coef0$q]
            if (ma_valido(theta.star)) break
          }

          serie.b <- simula_ma(N_INICIAL, theta.star)

          serie.colada <- c(tail(serie0, N_INICIAL - sc), head(serie.b, sc))
          stopifnot(length(serie.colada) == N_INICIAL)

          fit.b <- fit_ma(serie.colada)
          if (!(inherits(fit.b, "Arima")) && is.na(fit.b)) {
            print("Ajuste MA(2) falhou no bootstrap, pulando iteração.")
            theta1.boot[b] <- NA
            theta2.boot[b] <- NA
            next
          }

          theta1.boot[b] <- coef(fit.b)[1]
          theta2.boot[b] <- coef(fit.b)[2]
        }

        theta1.boot <- na.omit(theta1.boot)
        theta2.boot <- na.omit(theta2.boot)

        coef.vcov <- matrix(
          c(
            var(theta1.boot), cov(theta1.boot, theta2.boot),
            cov(theta2.boot, theta1.boot), var(theta2.boot)
          ),
          nrow = 2, byrow = TRUE
        )

        coef.vcov.inv <- solve(coef.vcov)
        theta1.m <- mean(theta1.boot)
        theta2.m <- mean(theta2.boot)

        t2 <- numeric(length(theta1.boot))
        for (b in seq_along(theta1.boot)) {
          diff.b <- c(theta1.boot[b], theta2.boot[b]) - c(theta1.m, theta2.m)
          t2[b] <- t(diff.b) %*% coef.vcov.inv %*% diff.b
        }

        ## 3) Limites de controle
        quant.boot <- quantile(t2, probs = c(0.025, 0.975))

        # 4) T² para série de controle
        t2.controle <- (
          t(c(theta1, theta2) - c(theta1.m, theta2.m)) %*%
            coef.vcov.inv %*%
            (c(theta1, theta2) - c(theta1.m, theta2.m))
        )[1]

        DADOS_OUT <- rbind(
          DADOS_OUT,
          data.table(
            n0 = N_INICIAL,
            n1 = sc,
            desvio = desvio,
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

# ------------------------------
# Resumo e gráfico
# ------------------------------

DADOS_RESUMO <- DADOS_OUT %>%
  group_by(n0, n1, desvio) %>%
  summarise(proporcao = mean(fora_de_controle)) %>%
  mutate(x_label = n1 + n0)

p1 <- ggplot(DADOS_RESUMO, aes(x = x_label, y = proporcao, color = factor(desvio))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  annotate("text",
    x = min(DADOS_RESUMO$x_label) - 5, y = 0.05, label = "5%",
    vjust = -0.2, hjust = 0.5, size = 4, color = "black"
  ) +
  labs(
    title = "Proporção de Séries Fora de Controle (MA(2))",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Proporção Fora de Controle",
    color = "Desvio de Θ1"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~n0,
    labeller = labeller(n0 = \(x) paste0("Tamanho Inicial: ", x)),
    scales = "free_x"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("controle_online_multiplos_ma2.png", p1, width = 10, height = 6, units = "in", dpi = 300)

print(p1)
