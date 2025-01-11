# Função auxiliar para gerar os dados
# Exemplo 1:
# ```r
# teste1 <- gerador_monte_carlo(
#   parametros = list(
#     lambda = seq(0.1, 0.4, by = 0.1)
#   ),
#   numero_de_execucoes = 2
# )
#
# teste1
# ```
#
# Exemplo 2:
# ```r
# teste2 <- gerador_monte_carlo(
#   parametros = expand.grid(
#     desvio_detectavel = seq(0.6, 1.8, by = 0.2),
#     intervalo_de_decisao = seq(4, 6, by = 1)
#   )
#   numero_de_execucoes = 2
# )
#
# teste2
# ```
#
# Total de execuções é:
# `execucoes = numero_de_execucoes * length(novos_phis) * nrow(parametros)`
#
gerador_monte_carlo <- function(parametros = NULL,
                                numero_de_execucoes = 100,
                                novos_phis = NULL) {
  # Verifica se `parametros` é um data frame, caso contrário converte para um data frame
  if (!is.null(parametros) && !is.data.frame(parametros)) {
    parametros <- as.data.frame(parametros)
  }
  n_par <- ifelse(is.null(parametros), 1, nrow(parametros))

  if (is.null(novos_phis)) {
    novos_phis <- seq(0.2, 0.6, by = 0.1)
  }

  total_execucoes <- numero_de_execucoes * length(novos_phis) * n_par
  print(paste("Total de execuções:", total_execucoes))

  # Expande a grade de parâmetros para cobrir todas as combinações
  result <- expand.grid(
    # Quantidade de execuções para cada combinação de parâmetros
    k = 1:numero_de_execucoes,
    # Parâmetros da Fase II
    f2_phi = novos_phis
  )

  if (!is.null(parametros)) {
    result <- merge(result, parametros, all = TRUE)
  }

  result %>%
    mutate(id = row_number()) %>%
    rowwise() %>%
    mutate(
      # Fase I
      ## Gera amostra de controle
      f1_amostras = list(barma.sim(n = n1, phi = phi_parametro, seed = id)),
      ## Estima o valor de phi da amostra de controle
      f1_phi = barma.phi_estimado(f1_amostras),
      ## Calcula os resíduos da amostra de controle
      f1_residuos = list(barma.residuos(f1_amostras, f1_phi)),
      # Fase IIf
      ## Gera a amostra subsequente
      f2_controle = list(
        barma.sim(
          n = n2,
          phi = f2_phi,
          # Última observação da amostra de controle
          y.start = f1_amostras[n1],
          # Define a semente para garantir a reprodutibilidade
          seed = id + 1337E4
        )
      ),
      ## Calcula os resíduos da amostra subsequente
      f2_residuos = list(barma.residuos(f2_controle, f1_phi))
    )
}
