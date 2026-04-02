Sys.setlocale("LC_CTYPE", "Portuguese_Brazil.utf8")

library(readr)
library(ggplot2)
library(dplyr)
library(scales)

sem_transformacao <- read_delim(
  "estudos/transformacao-no-espaco-parametrico/controle-sem-transformacao.csv",
  delim = ";"
)

com_transformacao <- read_delim(
  "estudos/transformacao-no-espaco-parametrico/controle-com-transformacao.csv",
  delim = ";"
)

dados <- bind_rows(sem_transformacao, com_transformacao) %>%
  mutate(
    metodo = factor(
      metodo,
      levels = c("Sem transformação", "Espaço transformado")
    )
  )

MC <- 800
B <- 500

dados %>%
  ggplot(aes(x = n1, y = proporcao, color = label, fill = label, group = label)) +
  # geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 5, alpha = 0.6) +
  geom_ribbon(aes(ymin = ic_inf, ymax = ic_sup), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  # geom_linerange(aes(ymin = ic_inf, ymax = ic_sup), linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  labs(
    title = "Monte Carlo para ARMA(1,1)",
    subtitle = sprintf("%d iterações Monte Carlo, %d bootstrap por iteração", MC, B),
    x = "Número de observações da Fase II",
    y = "Proporção fora de controle",
    color = "Parâmetros na Fase II (Φ; Θ)",
    fill = "Parâmetros na Fase II (Φ; Θ)"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_brewer(
    palette = "Dark2",
    labels = label_wrap(20),
    guide = guide_legend(nrow = 2)
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = label_wrap(20),
    guide = "none"
  ) +
  facet_grid(
    n0 + parametro_real ~ metodo,
    labeller = labeller(
      n0 = \(x) paste0("Tamanho Inicial: ", x),
      parametro_real = \(x) paste0("Parâmetros na Fase I: ", x),
      metodo = label_value
    ),
    scales = "fixed"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    text = element_text(size = 16)
  )