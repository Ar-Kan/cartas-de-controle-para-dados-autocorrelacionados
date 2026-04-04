# Gera coeficientes ARMA(1, 1) com dist N(coef_fase1, vcov_fase1)
gera_coeficientes <- function(
  coef_fase1,
  vcov_fase1,
  usar_transformacao
) {
  if (usar_transformacao) {
    # Gera coeficientes bootstrap válidos para ARMA(1,1)
    coef_star <- amostrar_coef_validos(coef = coef_fase1, matriz_vcov = vcov_fase1)
    if (is.null(coef_star)) return(NULL)

    phi_star <- unname(coef_star["phi_star"])
    theta_star <- unname(coef_star["theta_star"])
  } else {
    # Gera coeficientes bootstrap sem restrição, o que pode levar a valores inválidos de φ* e θ*.
    coef_star <- MASS::mvrnorm(1, mu = coef_fase1, Sigma = vcov_fase1) |> tryNull()
    if (is.null(coef_star) || any(!is.finite(coef_star))) return(NULL)

    phi_star <- unname(coef_star["ar1"])
    theta_star <- unname(coef_star["ma1"])

    # Valida os coeficientes bootstrap gerados
    if (!ar_valido(phi_star) || !ma_valido(theta_star)) return(NULL)
  }

  list(phi_star = phi_star, theta_star = theta_star)
}


executa_bootstrap <- function(
  serie_fase1,
  modelo_fase1,
  numero_de_novas_observacoes,
  numero_de_boots,
  usar_transformacao
) {
  phi_boot <- rep(NA_real_, numero_de_boots)
  theta_boot <- rep(NA_real_, numero_de_boots)
  cf1 <- coef(modelo_fase1)
  vf1 <- vcov(modelo_fase1)

  for (b in seq_len(numero_de_boots)) {
    coef_star <- gera_coeficientes(
      coef_fase1 = cf1,
      vcov_fase1 = vf1,
      usar_transformacao = usar_transformacao
    )
    if (is.null(coef_star)) next
    phi_star <- coef_star$phi_star
    theta_star <- coef_star$theta_star

    # Simula série Fase II*
    serie_fase2_b <- simula_arma(length(serie_fase1), phi_star, theta_star) |> tryNull()
    if (is.null(serie_fase2_b)) next

    serie_colada_b <- cola_series(
      serie_fase1 = serie_fase1,
      serie_fase2 = serie_fase2_b,
      numero_de_novas_observacoes = numero_de_novas_observacoes
    )

    ajuste_b <- fit_arma(serie_colada_b, cf1["ar1"], cf1["ma1"], transform.pars = usar_transformacao)
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
