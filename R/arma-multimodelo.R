# Função para ajustar um modelo ARMA(p, q) específico
# Para faciliar mocks
fit_modelo <- function(serie, p, q) {
  forecast::Arima(serie, order = c(p, 0, q), include.mean = FALSE)
}

# Função para gerar combinações de p e q
# p e q variam de 0 a max_p e max_q, respectivamente, mas excluindo o caso p=0 e q=0
graus_de_liberdade <- function(max_p, max_q) {
  modelos_df <- expand.grid(p = 0:max_p, q = 0:max_q)
  modelos_df[!(modelos_df$p == 0 & modelos_df$q == 0),]
}

# Função para calcular pesos de Akaike a partir de um vetor de critérios de informação (AIC, BIC, etc.)
# A função lida com valores de IC que podem ser Inf ou NA, atribuindo peso zero a esses casos e normalizando os pesos válidos para somarem 1.
# O cálculo é baseado na diferença entre cada IC e o menor IC válido, utilizando a fórmula:
# \deqn{
#   w_i = \frac{
#     \exp\left(-\frac{IC_i - IC_{min}}{2}\right)
#   }{
#     \sum_{j} \exp\left(-\frac{IC_j - IC_{min}}{2}\right)
#   }
# }
# onde \eqn{IC_{min}} é o menor valor de IC entre os modelos válidos.
# A função retorna um vetor de pesos alinhado com o vetor original de ICs, onde os pesos correspondentes a ICs inválidos (Inf ou NA) são zero.
calcular_pesos <- function(ics) {
  # Remove NA explicitamente (são inválidos)
  validos <- is.finite(ics)

  if (!any(validos)) {
    stop("Nenhum IC válido para calcular pesos")
  }

  ics_validos <- ics[validos]

  delta <- ics_validos - min(ics_validos)

  pesos_validos <- exp(-delta / 2)
  pesos_validos <- pesos_validos / sum(pesos_validos)

  # mantém o alinhamento com o vetor original de ICs
  # preenche pesos inválidos com 0
  pesos <- rep(0, length(ics))
  pesos[validos] <- pesos_validos

  return(pesos)
}

# Função para ajustar modelos ARMA(p, q) e calcular o critério de informação
# Retorna uma lista com os modelos ajustados, seus critérios de informação, pesos de Akaike e os graus de liberdade
# Se cutoff for fornecido, filtra os modelos com base na diferença do critério de informação em relação ao melhor modelo e recalcula os pesos de Akaike para os modelos filtrados
ajustar_modelos <- function(serie, max_p, max_q, criterio = "aic", cutoff = NULL) {
  modelos_df <- graus_de_liberdade(max_p, max_q)

  ics <- rep(Inf, nrow(modelos_df))
  resultados <- vector("list", nrow(modelos_df))

  #i=1
  for (i in 1:nrow(modelos_df)) {
    p <- modelos_df$p[i]
    q <- modelos_df$q[i]

    tryCatch({
      modelo <- fit_modelo(serie, p, q)
      resultados[[i]] <- modelo
      ics[i] <- modelo[[criterio]]
    }, error = function(e) {
    })
  }

  # Calcular pesos de Akaike
  pesos <- calcular_pesos(ics)

  if (is.null(cutoff)) {
    list(
      modelos = resultados,
      ics = ics,
      pesos = pesos,
      modelos_df = modelos_df
    )
  }
  else {
    delta_ic <- ics - min(ics, na.rm = TRUE)
    substantial <- which(delta_ic < cutoff)

    ics_filtrados <- ics[substantial]
    pesos_filtrados <- calcular_pesos(ics_filtrados)

    list(
      modelos = resultados[substantial],
      ics = ics_filtrados,
      pesos = pesos_filtrados,
      modelos_df = modelos_df[substantial,]
    )
  }
}
