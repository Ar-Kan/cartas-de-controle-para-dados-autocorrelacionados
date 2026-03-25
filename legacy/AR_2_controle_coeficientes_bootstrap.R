library(ggplot2)
library(dplyr)
library(data.table)
library(MASS)

set.seed(42)

# ------------------------------
# Parâmetros gerais
# ------------------------------
PHI_REAL <- c(0.2, -0.3) # valores verdadeiros de Φ1 e Φ2
tamanhos_amostrais_iniciais <- c(100, 200)
MC <- 300
B <- 200

DESVIOS_PHI <- c(0, 0.2, 0.6)

passo <- 13
bc <- c(1, 3, 7)
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

# Simula AR(2)
simula_ar <- function(n, phi, sd = NULL) {
  args <- list(n = n, model = list(ar = phi))
  if (!is.null(sd)) {
    args$sd <- sd
  }
  do.call(arima.sim, args)
}

# Ajusta AR(2)
fit_ar <- function(serie) {
  tryCatch(
    {
      return(
        arima(
          serie,
          order = c(2, 0, 0), # AR(2)
          include.mean = FALSE,
          method = "ML"
        )
      )
    },
    error = \(e) {
      warning("Erro ao ajustar AR(2): ", e)
      return(NA)
    }
  )
}

# Testa estacionaridade
ar_valido <- function(phi) {
  stab <- all(Mod(polyroot(c(1, -phi))) > 1)
  stab
}

# Extrai coeficientes e matriz de covariância do modelo AR
ar_coef <- function(modelo) {
  coefs <- coef(modelo)
  list(
    p = length(grep("^ar", names(coefs))),
    coef = coefs,
    vcov = vcov(modelo),
    vcov_inv = solve(vcov(modelo))
  )
}

# ------------------------------
# Loop Monte Carlo
# ------------------------------

print(paste0("Iniciando simulações com ", MC, " iterações Monte Carlo..."))

for (mc in seq_len(MC)) {
  if (mc %% 10 == 0) {
    print(paste0("Iteração Monte Carlo: ", mc, "/", MC, " (", round(mc / MC * 100, 2), "% completo)"))
  }

  for (N_INICIAL in tamanhos_amostrais_iniciais) {
    ## 1) Simula e ajusta a série inicial
    serie0 <- simula_ar(N_INICIAL, PHI_REAL)
    fit0 <- fit_ar(serie0)
    coef0 <- ar_coef(fit0)

    for (desvio in DESVIOS_PHI) {
      # Próximas amostras
      serie1 <- simula_ar(N_INICIAL, PHI_REAL + c(desvio, 0))

      for (sc in SEQ_COLAGENS) {
        if (sc >= N_INICIAL) break

        # Série de controle
        serie1.controle <- c(tail(serie0, N_INICIAL - sc), head(serie1, sc))
        stopifnot(length(serie1.controle) == N_INICIAL)

        # Ajusta série colada
        fit_tmp <- fit_ar(serie1.controle)
        phi1 <- coef(fit_tmp)[1]
        phi2 <- coef(fit_tmp)[2]
        if (is.na(phi1) || is.na(phi2)) {
          print("Ajuste AR(2) falhou, pulando iteração.")
          next
        }

        ## 2) Bootstrap paramétrico
        phi1.boot <- numeric(B)
        phi2.boot <- numeric(B)
        for (b in seq_len(B)) {
          repeat {
            coef.star <- mvrnorm(1, mu = coef0$coef, Sigma = coef0$vcov)
            phi.star <- coef.star[1:coef0$p]
            if (ar_valido(phi.star)) break
          }

          serie.b <- simula_ar(N_INICIAL, phi.star)

          serie.colada <- c(tail(serie0, N_INICIAL - sc), head(serie.b, sc))
          stopifnot(length(serie.colada) == N_INICIAL)

          fit.b <- fit_ar(serie.colada)
          if (!(inherits(fit.b, "Arima")) && is.na(fit.b)) {
            print("Ajuste AR(2) falhou no bootstrap, pulando iteração.")
            phi1.boot[b] <- NA
            phi2.boot[b] <- NA
            next
          }

          phi1.boot[b] <- coef(fit.b)[1]
          phi2.boot[b] <- coef(fit.b)[2]
        }

        phi1.boot <- na.omit(phi1.boot)
        phi2.boot <- na.omit(phi2.boot)

        coef.vcov <- matrix(
          c(
            var(phi1.boot), cov(phi1.boot, phi2.boot),
            cov(phi2.boot, phi1.boot), var(phi2.boot)
          ),
          nrow = 2, byrow = TRUE
        )

        coef.vcov.inv <- solve(coef.vcov)
        phi1.m <- mean(phi1.boot)
        phi2.m <- mean(phi2.boot)

        t2 <- numeric(length(phi1.boot))
        for (b in seq_along(phi1.boot)) {
          diff.b <- c(phi1.boot[b], phi2.boot[b]) - c(phi1.m, phi2.m)
          t2[b] <- t(diff.b) %*% coef.vcov.inv %*% diff.b
        }

        ## 3) Limites de controle
        quant.boot <- quantile(t2, probs = c(0.025, 0.975))

        # 4) T² para série de controle
        t2.controle <- (
          t(c(phi1, phi2) - c(phi1.m, phi2.m)) %*%
            coef.vcov.inv %*%
            (c(phi1, phi2) - c(phi1.m, phi2.m))
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
    title = "Proporção de Séries Fora de Controle (AR(2))",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Proporção Fora de Controle",
    color = "Desvio de Φ1"
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

ggsave("controle_online_multiplos_ar2__.png", p1, width = 10, height = 6, units = "in", dpi = 300)

print(p1)
