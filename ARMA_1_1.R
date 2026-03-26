library(ggplot2)
library(dplyr)
library(data.table)
library(MASS)

set.seed(42)

# ------------------------------
# Parâmetros gerais
# ------------------------------
PHI_REAL_LISTA <- c(0.05, 0.3, 0.6) # valor verdadeiro de Φ1
THETA_REAL_LISTA <- c(0.05, 0.3, 0.6) # valor verdadeiro de Θ1

tamanhos_amostrais_iniciais <- c(100)
MC <- 500
B <- 500

DESVIOS_PHI <- c(0) # , 0.2) # , 0.6)
DESVIOS_THETA <- c(0) # , 0.2) # , 0.6)

# passo <- 20
# bc <- c(1, 5)
# SEQ_COLAGENS <- unlist(lapply(seq(0, 100, passo), function(x) bc + x))
# SEQ_COLAGENS <- c(SEQ_COLAGENS, 135, 150, 165, 175, 185, 199)
SEQ_COLAGENS <- c(1, 25, 50, 75, 99, 200, 300, 400, 499)

DADOS_OUT <- data.table(
  phi = numeric(),
  theta = numeric(),
  n0 = numeric(),
  n1 = numeric(),
  desvio_phi = numeric(),
  desvio_theta = numeric(),
  fora_de_controle = logical(),
  t2.inf = numeric(),
  t2.sup = numeric(),
  t2.controle = numeric()
)

# ------------------------------
# Funções auxiliares
# ------------------------------

# Simula ARMA(1,1)
simula_arma <- function(n, phi, theta, sd = NULL) {
  args <- list(n = 200 + n, model = list(ar = phi, ma = theta))
  if (!is.null(sd)) args$sd <- sd
  obs <- do.call(arima.sim, args)
  return(tail(obs, n))
}

# Ajusta ARMA(1,1)
fit_arma <- function(serie) {
  tryCatch(
  {
    arima(
      serie,
      order = c(1, 0, 1), # ARMA(1,1)
      include.mean = FALSE,
      method = "ML"
    )
  },
    error = function(e) {
      warning("Erro ao ajustar ARMA(1,1): ", e)
      return(NA)
    }
  )
}

# Testa estacionaridade (para φ)
ar_valido <- function(phi) {
  all(Mod(polyroot(c(1, -phi))) > 1)
}

# Testa invertibilidade (para θ)
ma_valido <- function(theta) {
  all(Mod(polyroot(c(1, theta))) > 1)
}

# Extrai coeficientes e matriz de covariância
arma_coef <- function(modelo) {
  coefs <- coef(modelo)
  list(
    coef = coefs,
    vcov = vcov(modelo),
    vcov_inv = solve(vcov(modelo))
  )
}

# ------------------------------
# Loop Monte Carlo
# ------------------------------

print(paste0("Iniciando simulações ARMA(1,1) com ", MC, " iterações Monte Carlo..."))

for (PHI_REAL in PHI_REAL_LISTA) {
  for (THETA_REAL in THETA_REAL_LISTA) {
    print(paste0("Simulando para Φ0 = ", PHI_REAL, " e Θ0 = ", THETA_REAL))

    for (mc in seq_len(MC)) {
      if (mc %% 10 == 0) {
        print(paste0("Iteração Monte Carlo: ", mc, "/", MC))
      }

      for (N_INICIAL in tamanhos_amostrais_iniciais) {
        ## 1) Série base
        serie0 <- simula_arma(N_INICIAL, PHI_REAL, THETA_REAL)
        fit0 <- fit_arma(serie0)
        coef0 <- arma_coef(fit0)

        for (desvio_phi in DESVIOS_PHI) {
          for (desvio_theta in DESVIOS_THETA) {
            # Nova série com desvio em φ e θ
            serie1 <- simula_arma(N_INICIAL, PHI_REAL + desvio_phi, THETA_REAL + desvio_theta)

            for (sc in SEQ_COLAGENS) {
              if (sc >= N_INICIAL) break

              # print(paste0(
              #   "MC=", mc,
              #   ", N_inicial=", N_INICIAL,
              #   ", sc=", sc,
              #   ", desvio_phi=", desvio_phi,
              #   ", desvio_theta=", desvio_theta
              # ))
              # Série colada (controle)
              serie1.controle <- c(tail(serie0, N_INICIAL - sc), head(serie1, sc))
              stopifnot(length(serie1.controle) == N_INICIAL)

              fit_tmp <- fit_arma(serie1.controle)
              if (is.na(fit_tmp)[1]) next

              phi_hat <- coef(fit_tmp)["ar1"]
              theta_hat <- coef(fit_tmp)["ma1"]
              if (is.na(phi_hat) || is.na(theta_hat)) next

              ## 2) Bootstrap paramétrico
              phi.boot <- numeric(B)
              theta.boot <- numeric(B)
              for (b in seq_len(B)) {
                tryCatch(
                {
                  repeat {
                    coef.star <- mvrnorm(1, mu = coef0$coef, Sigma = coef0$vcov)
                    phi.star <- coef.star[1]
                    theta.star <- coef.star[2]
                    if (ar_valido(phi.star) && ma_valido(theta.star)) break
                  }
                },
                  error = function(e) {
                    ## 'Sigma' is not positive definite
                    warning("Erro na geração de coeficientes bootstrap: ", e)
                    next
                  }
                )

                serie.b <- simula_arma(N_INICIAL, phi.star, theta.star)
                serie.colada <- c(tail(serie0, N_INICIAL - sc), head(serie.b, sc))
                stopifnot(length(serie.colada) == N_INICIAL)

                fit.b <- fit_arma(serie.colada)
                if (is.na(fit.b)[1]) {
                  phi.boot[b] <- NA
                  theta.boot[b] <- NA
                  next
                }
                if (
                  !(inherits(fit.b, "Arima") || inherits(fit.b, "ar")) && is.na(fit.b)
                ) {
                  print("Ajuste ARMA falhou no bootstrap, pulando iteração.")
                  phi.boot[b] <- NA
                  theta.boot[b] <- NA
                  next
                }

                phi.boot[b] <- coef(fit.b)["ar1"]
                theta.boot[b] <- coef(fit.b)["ma1"]
              }

              phi.boot <- na.omit(phi.boot)
              theta.boot <- na.omit(theta.boot)
              if (length(phi.boot) < 5 || length(theta.boot) < 5) next

              coef.vcov <- matrix(
                c(
                  var(phi.boot), cov(phi.boot, theta.boot),
                  cov(theta.boot, phi.boot), var(theta.boot)
                ),
                nrow = 2, byrow = TRUE
              )

              coef.vcov.inv <- solve(coef.vcov)
              phi.m <- mean(phi.boot)
              theta.m <- mean(theta.boot)

              t2 <- numeric(length(phi.boot))
              for (b in seq_along(phi.boot)) {
                diff.b <- c(phi.boot[b], theta.boot[b]) - c(phi.m, theta.m)
                t2[b] <- t(diff.b) %*% coef.vcov.inv %*% diff.b
              }

              ## 3) Limites de controle
              quant.boot <- quantile(t2, probs = c(0.025, 0.975))

              ## 4) Estatística T²
              t2.controle <- (
                t(c(phi_hat, theta_hat) - c(phi.m, theta.m)) %*%
                  coef.vcov.inv %*%
                  (c(phi_hat, theta_hat) - c(phi.m, theta.m))
              )[1]

              DADOS_OUT <- rbind(
                DADOS_OUT,
                data.table(
                  phi = PHI_REAL + desvio_phi,
                  theta = THETA_REAL + desvio_theta,
                  n0 = N_INICIAL,
                  n1 = sc,
                  desvio_phi = desvio_phi,
                  desvio_theta = desvio_theta,
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
    }
  }
}

# ------------------------------
# Resumo e gráfico
# ------------------------------

DADOS_RESUMO <- DADOS_OUT %>%
  group_by(n0, n1, desvio_phi, desvio_theta) %>%
  summarise(proporcao = mean(fora_de_controle), .groups = "drop") %>%
  mutate(x_label = n1 + n0)

p1 <- ggplot(DADOS_RESUMO, aes(
  x = x_label, y = proporcao,
  color = interaction(desvio_phi, desvio_theta, sep = ";")
)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  labs(
    title = "Proporção de Séries Fora de Controle (ARMA(1,1)).png",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Proporção Fora de Controle",
    color = "Desvios (Φ1; Θ1)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~n0,
             labeller = labeller(n0 = \(x) paste0("Tamanho Inicial: ", x)),
             scales = "free_x"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ggsave("controle_online_arma11_20251115.png", p1, width = 10, height = 6, units = "in", dpi = 300)
print(p1)


DADOS_OUT %>%
  mutate(
    label = paste0(
      "Φ = ", round(phi, 2),
      "; Θ = ", round(theta, 2)
    )
  ) %>%
  group_by(label, n0, n1) %>%
  summarise(proporcao = mean(fora_de_controle), .groups = "drop") %>%
  mutate(
    x_label = n1 + n0,
    se = sqrt(proporcao * (1 - proporcao) / MC),
    ic_inf = pmax(0, proporcao - 1.96 * se),
    ic_sup = pmin(1, proporcao + 1.96 * se)
  ) %>%
  ggplot(aes(x = x_label, y = proporcao, color = label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  # geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 5) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  labs(
    title = "Proporção de Séries Fora de Controle (ARMA(1,1)).png",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Proporção Fora de Controle",
    color = "Parâmetros (Φ; Θ)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(
    palette = "Dark2",
    labels = scales::label_wrap(20),
    guide = guide_legend(nrow = 2)
  ) +
  facet_wrap(~n0,
             labeller = labeller(n0 = \(x) paste0("Tamanho Inicial: ", x)),
             scales = "free_x"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


DADOS_RESUMO__ <- DADOS_OUT %>%
  group_by(n0, n1, desvio_phi, desvio_theta) %>%
  summarise(
    proporcao = mean(fora_de_controle),
    .groups = "drop"
  ) %>%
  mutate(
    x_label = n0 + n1,
    se = sqrt(proporcao * (1 - proporcao) / MC),
    ic_inf = pmax(0, proporcao - 1.96 * se),
    ic_sup = pmin(1, proporcao + 1.96 * se)
  )

p1__ <- ggplot(DADOS_RESUMO__, aes(
  x = x_label, y = proporcao,
  color = interaction(desvio_phi, desvio_theta, sep = ";")
)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 5) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  labs(
    title = "Proporção de Séries Fora de Controle (ARMA(1,1)).png",
    x = "Tamanho da Série Colada (n0 + n1)",
    y = "Proporção Fora de Controle",
    color = "Desvios (Φ1; Θ1)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(
    ~n0,
    labeller = labeller(n0 = \(x) paste0("Tamanho Inicial: ", x)),
    scales = "free_x"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p1__
