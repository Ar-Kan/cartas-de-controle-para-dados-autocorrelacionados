# Calcula o T² de Hotelling para um par (φ, θ) usando a distribuição bootstrap de T² sob H0
t2_hotelling_boot <- function(phi_boot, theta_boot, phi, theta, vcov_inv) {
  phi_boot_media <- mean(phi_boot)
  theta_boot_media <- mean(theta_boot)

  # Calcula T² para cada par (φ*, θ*) do bootstrap
  # T²(x) = (x - mu_boot)' S^{-1} (x - mu_boot)
  t2_boot <- vapply(
    seq_along(phi_boot),
    function(i) {
      diff_b <- c(phi_boot[i], theta_boot[i]) - c(phi_boot_media, theta_boot_media)
      as.numeric(t(diff_b) %*% vcov_inv %*% diff_b)
    },
    numeric(1)
  )

  # Obtém o limite de controle q95(T²)
  t2_limite <- quantile(t2_boot, probs = 0.95, names = FALSE)
  stopifnot(
    "O quantil de controle T² retornou NA." = !is.na(t2_limite)
  )

  diff_controle <- c(phi, theta) - c(phi_boot_media, theta_boot_media)
  t2_controle <- as.numeric(
    t(diff_controle) %*%
      vcov_inv %*%
      diff_controle
  )

  esta_fora_de_controle <- t2_controle > t2_limite

  data.table::data.table(
    numero_de_bootstrap_validos = length(t2_boot),
    fora_de_controle = esta_fora_de_controle,
    t2_limite = t2_limite,
    t2_controle = t2_controle
  )
}
