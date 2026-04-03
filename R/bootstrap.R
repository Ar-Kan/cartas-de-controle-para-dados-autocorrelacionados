executa_bootstrap <- function(
  serie0,
  coef0,
  n_inicial,
  numero_de_novas_observacoes,
  numero_de_boots,
  usar_transformacao
) {
  phi_boot <- rep(NA_real_, numero_de_boots)
  theta_boot <- rep(NA_real_, numero_de_boots)

  for (b in seq_len(numero_de_boots)) {
    # Gera coeficientes com dist N(coef0$coef, coef0$vcov)
    if (usar_transformacao) {
      # Gera coeficientes bootstrap válidos para ARMA(1,1)
      coef_star <- bootstrap_coef_validos(coef = coef0$coef, matriz_vcov = coef0$vcov)
      if (is.null(coef_star)) next

      phi_star <- unname(coef_star["phi_star"])
      theta_star <- unname(coef_star["theta_star"])
    } else {
      # Gera coeficientes bootstrap sem restrição, o que pode levar a valores inválidos de φ* e θ*.
      coef_star <- MASS::mvrnorm(1, mu = coef0$coef, Sigma = coef0$vcov) |> tryNull()
      if (is.null(coef_star) || any(!is.finite(coef_star))) next

      phi_star <- unname(coef_star["ar1"])
      theta_star <- unname(coef_star["ma1"])

      # Valida os coeficientes bootstrap gerados
      if (!ar_valido(phi_star) || !ma_valido(theta_star)) next
    }

    # Simula série Fase II*
    serie1_b <- simula_arma(n_inicial, phi_star, theta_star) |> tryNull()
    if (is.null(serie1_b)) next

    serie_colada_b <- cola_series(
      serie_fase1 = serie0,
      serie_fase2 = serie1_b,
      numero_de_novas_observacoes = numero_de_novas_observacoes
    )

    ajuste_b <- fit_arma(serie_colada_b, coef0$coef["ar1"], coef0$coef["ma1"], transform.pars = usar_transformacao)
    if (is.null(ajuste_b) || !ajuste_b$convergiu) next
    fit_b <- ajuste_b$fit

    coef_b <- coef(fit_b) |> tryNull()
    if (is.null(coef_b)) next

    phi_boot[b] <- unname(coef_b["ar1"])
    theta_boot[b] <- unname(coef_b["ma1"])
  }

  validos <- is.finite(phi_boot) & is.finite(theta_boot)
  phi_boot <- phi_boot[validos]
  theta_boot <- theta_boot[validos]

  if (length(phi_boot) < numero_de_boots * 0.9) {
    warning(sprintf(
      "Número de bootstrap válidos (%d) é menor que 90%% do total (%d)", length(phi_boot), numero_de_boots
    ))
  }

  list(phi_boot = phi_boot, theta_boot = theta_boot)
}
