rm(list = ls())
options(error = NULL)
Sys.setlocale("LC_CTYPE", "Portuguese_Brazil.utf8")

library(ggplot2)
library(dplyr)

source("R/utils.R")
source("R/paralelizacao.R")
source("R/arma.R")
source("R/bootstrap.R")
source("R/monte_carlo.R")

set.seed(42)

######################################
# PARÂMETROS
PHI_REAL_LISTA <- c(0.2)
THETA_REAL_LISTA <- c(0.5)

TAMANHOS_AMOSTRAIS_INICIAIS <- c(100)

MC <- 500
B <- 500

DESVIOS_PHI <- c(0, 0.2)
DESVIOS_THETA <- c(0)


#####################################################################################
# LOOP PRINCIPAL

message(sprintf(
  "Iniciando simulações ARMA(1,1) com %d iterações Monte Carlo...",
  MC
))


system.time(
  for (phi_real in PHI_REAL_LISTA) {
    for (theta_real in THETA_REAL_LISTA) {
      DADOS_OUT <- execucao_paralela(
        n_execucoes = MC,
        funcao = executa_um_mc,
        lista_argumentos = list(
          phi_real = phi_real,
          theta_real = theta_real,
          tamanhos_amostrais_iniciais = TAMANHOS_AMOSTRAIS_INICIAIS,
          numero_de_boots = B,
          desvios_phi = DESVIOS_PHI,
          desvios_theta = DESVIOS_THETA,
          usar_transformacao = TRUE
        ),
      )
    }
  }
)

###############################################################
# Resultados

message("Concluído com sucesso.")


DADOS_PLOT <- DADOS_OUT %>%
  mutate(
    parametro_real = paste0("(", phi_real, ", ", theta_real, ")"),
    label = paste0("(", phi_alt, ", ", theta_alt, ")")
  ) %>%
  group_by(parametro_real, label, n0, n1) %>%
  summarise(proporcao = mean(fora_de_controle), .groups = "drop") %>%
  mutate(
    se = sqrt(proporcao * (1 - proporcao) / MC),
    ic_inf = pmax(0, proporcao - 1.96 * se),
    ic_sup = pmin(1, proporcao + 1.96 * se)
  )

p <- DADOS_PLOT %>%
  ggplot(aes(x = n1, y = proporcao, color = label, fill = label, group = label)) +
  # geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 5, alpha = 0.6) +
  geom_ribbon(aes(ymin = ic_inf, ymax = ic_sup), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_linerange(aes(ymin = ic_inf, ymax = ic_sup), linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 2) +
  # geom_pointrange(aes(ymin = ic_inf, ymax = ic_sup), linewidth = 0.5, fatten = 1.3) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  labs(
    title = "Monte Carlo para ARMA(1,1)",
    subtitle = sprintf("%d iterações Monte Carlo, %d bootstrap por iteração", MC, B),
    x = "Número de observações da Fase II",
    y = "Proporção fora de controle",
    color = "Parâmetros na Fase II (Φ; Θ)",
    fill = "Parâmetros na Fase II (Φ; Θ)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(
    palette = "Dark2",
    labels = scales::label_wrap(20),
    guide = guide_legend(nrow = 2)
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = scales::label_wrap(20),
    guide = "none"
  ) +
  facet_wrap(
    ~n0 + parametro_real,
    labeller = labeller(
      n0 = \(x) paste0("Tamanho Inicial: ", x),
      parametro_real = \(x) paste0("Parâmetros na Fase I: ", x)
    ),
    scales = "free_x"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16)
  )

print(p)
