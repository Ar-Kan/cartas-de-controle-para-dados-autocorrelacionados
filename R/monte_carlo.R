seq_novas_observacoes <- function(n_inicial) {
  c(1, round(n_inicial / 2), n_inicial - 1, n_inicial + round(n_inicial / 2))
}

avalia_um_cenario_fase2 <- function(
  serie_fase1,
  serie_fase2,
  modelo_fase1,
  numero_de_novas_observacoes,
  numero_de_boots,
  usar_transformacao
) {
  serie_fase2_controle <- cola_series(
    serie_fase1 = serie_fase1,
    serie_fase2 = serie_fase2,
    numero_de_novas_observacoes = numero_de_novas_observacoes
  )

  # Série de controle: colagem da série de controle da Fase I com as n primeiras observações da Fase II
  cf1 <- coef(modelo_fase1)
  ajuste_fase2_controle <- fit_arma(serie_fase2_controle, cf1["ar1"], cf1["ma1"], transform.pars = usar_transformacao)
  if (is.null(ajuste_fase2_controle) || !ajuste_fase2_controle$convergiu) next
  modelo_fase2_controle <- ajuste_fase2_controle$fit

  coef1_controle <- coef(modelo_fase2_controle) |> tryNull()
  if (is.null(coef1_controle)) next

  coef1_phi <- unname(coef1_controle["ar1"])
  coef1_theta <- unname(coef1_controle["ma1"])
  if (!is.finite(coef1_phi) || !is.finite(coef1_theta)) next

  # Bootstrap para obter distribuição de T² sob H0 (sistema em controle)
  boot <- executa_bootstrap(
    serie_fase1 = serie_fase1,
    modelo_fase1 = modelo_fase1,
    numero_de_novas_observacoes = numero_de_novas_observacoes,
    numero_de_boots = numero_de_boots,
    usar_transformacao = usar_transformacao
  )

  # Matriz da informação de Fisher da série de controle
  amostra_vcov <- modelo_fase2_controle$var.coef

  amostra_vcov_inv <- solve(amostra_vcov) |> tryNull()
  if (is.null(amostra_vcov_inv)) next

  t2_resultado <- t2_hotelling_boot(
    phi_boot = boot$phi_boot,
    theta_boot = boot$theta_boot,
    phi = coef1_phi,
    theta = coef1_theta,
    vcov_inv = amostra_vcov_inv
  )

  t2_resultado
}


executa_um_mc <- function(
  execucao,
  phi_real,
  theta_real,
  tamanhos_amostrais_iniciais,
  numero_de_boots,
  desvios_phi,
  desvios_theta,
  usar_transformacao
) {
  resultados_mc <- list()
  idx_res <- 1L

  for (n_inicial in tamanhos_amostrais_iniciais) {
    # Simula série original sob controle (Fase I)
    # Observações nas quais que se assume que o sistema está em controle
    serie_fase1 <- simula_arma(n_inicial, phi_real, theta_real) |> tryNull()
    if (is.null(serie_fase1)) next

    ajuste_fase1 <- fit_arma(serie_fase1, transform.pars = usar_transformacao)
    if (is.null(ajuste_fase1) || !ajuste_fase1$convergiu) next
    modelo_fase1 <- ajuste_fase1$fit

    for (desvio_phi in desvios_phi) {
      for (desvio_theta in desvios_theta) {
        # Altera os parâmetros reais para simular a série alternativa (Fase II)
        phi_alt <- phi_real + desvio_phi
        theta_alt <- theta_real + desvio_theta

        if (!ar_valido(phi_alt) || !ma_valido(theta_alt)) next

        # Número de novas observações da Fase II
        seq_colagens <- seq_novas_observacoes(n_inicial)

        # Simula série alternativa (Fase II)
        # Observações que foram obtidas opós série de controle da Fase I
        serie_fase2 <- simula_arma(max(seq_colagens), phi_alt, theta_alt) |> tryNull()
        if (is.null(serie_fase2)) next

        for (numero_de_novas_observacoes in seq_colagens) {
          resultado_cenario <- avalia_um_cenario_fase2(
            serie_fase1 = serie_fase1,
            serie_fase2 = serie_fase2,
            modelo_fase1 = modelo_fase1,
            numero_de_novas_observacoes = numero_de_novas_observacoes,
            numero_de_boots = numero_de_boots,
            usar_transformacao = usar_transformacao
          )

          resultados_mc[[idx_res]] <- cbind(
            mc = execucao,
            phi_real = phi_real,
            theta_real = theta_real,
            phi_alt = phi_alt,
            theta_alt = theta_alt,
            n0 = n_inicial,
            n1 = numero_de_novas_observacoes,
            desvio_phi = desvio_phi,
            desvio_theta = desvio_theta,
            resultado_cenario
          )

          idx_res <- idx_res + 1L
        }
      }
    }
  }

  if (length(resultados_mc) == 0L) {
    return(data.table::data.table())
  }

  data.table::rbindlist(resultados_mc, fill = TRUE)
}
