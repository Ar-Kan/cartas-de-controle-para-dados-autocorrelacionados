Sys.setenv(lang = "en_US")
Sys.setlocale("LC_ALL", "en_US.UTF-8")
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  out.extra = "keepaspectratio=true",
  fig.align = "center",
  collapse = TRUE
)

required_packages <- c(
  "ggplot2",
  "plotly", # interactive plots
  "knitr",
  "gridExtra", # grid.arrange
  "tidyr",
  "dplyr",
  "purrr",
  "DT", # better html tables
  "qcc",
  "BTSR"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Temas para os gráficos

theme.base <- theme_minimal(base_size = 11) +
  theme(
    axis.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    plot.caption = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(
      colour = adjustcolor("grey90", alpha.f = 0.4),
      linewidth = 0.5
    ),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(colour = "grey"),
    axis.line.y = element_line(colour = "grey"),
  )

theme.no_legend <- theme(legend.position = "none")

plotly.base <- function(p) {
  p %>%
    layout(margin = list(b = 60, t = 80)) %>%
    config(mathjax = "cdn")
}

my.plotly <- function(p) {
  if (knitr::is_latex_output()) {
    return(p)
  }
  ggplotly(p) %>% plotly.base()
}

# Função para salvar e carregar resultados
cache_dados <- function(chave, funcao_geradora) {
  c <- paste("cache", chave, sep = "/")
  c <- paste(c, "rds", sep = ".")
  cache <- file.exists(c)

  if (!cache) {
    resultado <- funcao_geradora()
    saveRDS(resultado, c)
  } else {
    resultado <- readRDS(c)
  }

  resultado
}
