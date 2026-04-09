rm(list = ls())
options(error = NULL)
Sys.setlocale("LC_CTYPE", "Portuguese_Brazil.utf8")

library(ggplot2)
library(scales)
library(dplyr)

# Autocorrelação teórica do ARMA(1,1)
rho <- function(k, phi, theta, tol = 1e-10) {
  stopifnot(k >= 1)

  denominador <- 1 + 2 * phi * theta + theta^2

  if (abs(denominador) < tol) {
    return(NA_real_)
  }

  numerador <- phi^(k - 1) * (phi + theta) * (1 + phi * theta)

  numerador / denominador
}

# Grade apenas na região estacionária/invertível
phi_vals <- seq(-0.99, 0.99, by = 0.02)
theta_vals <- seq(-0.99, 0.99, by = 0.02)
k_vals <- 1:4

pho_long <- expand.grid(
  k = k_vals,
  phi = phi_vals,
  theta = theta_vals
)

pho_long$rho <- mapply(
  rho,
  k = pho_long$k,
  phi = pho_long$phi,
  theta = pho_long$theta
)

# Garante região válida explicitamente
pho_long <- subset(
  pho_long,
  abs(phi) < 1 & abs(theta) < 1
)

pho_long$k <- factor(pho_long$k)


p1 <- ggplot(pho_long, aes(x = theta, y = phi)) +
  geom_tile(aes(fill = rho)) +
  geom_contour(aes(z = rho), color = "white", alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~k, ncol = 2, labeller = label_bquote(k == .(k))) +
  coord_equal(expand = FALSE) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-1, 1),
    oob = squish,
    breaks = seq(-1, 1, by = 0.5),
    name = "ρ"
  ) +
  labs(
    title = "Mapa de calor e curvas de nível de ρ_k(φ, θ)",
    subtitle = "ARMA(1,1) na região |φ| < 1 e |θ| < 1",
    x = "θ",
    y = "φ"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

print(p1)

candidatos <- data.frame(
  phi = c(0.6, 0.7, 0.5, 0.4, 0.3, 0.6, -0.6),
  theta = c(0.2, 0.3, 0.3, 0.5, 0.2, -0.5, 0.2)
) %>%
  mutate(nome = paste0("(Φ=", phi, ", Θ=", theta, ")"))

rho_vec <- function(phi, theta, k_max = 10) {
  sapply(1:k_max, function(k) rho(k, phi, theta))
}

acf_long <- do.call(
  rbind,
  lapply(1:nrow(candidatos), function(i) {
    data.frame(
      lag = 1:10,
      rho = rho_vec(candidatos$phi[i], candidatos$theta[i]),
      nome = candidatos$nome[i]
    )
  })
)

p2 <- ggplot(acf_long, aes(x = lag, y = rho, color = nome)) +
  geom_line() +
  geom_point() +
  theme_minimal()
print(p2)

# ggsave("estudos/comportamento-rho/rho.png", p1, width = 10, height = 8)
# ggsave("estudos/comportamento-rho/candidatos_acf.png", p2, width = 10, height = 6)
