#
# Função para simulação de amostras de um AR(1)
# 
# @param parametros lista com os parâmetros do modelo AR(1):
# \describe{
#   \item{n}{número de observações}
#   \item{phi}{parâmetro de auto-regressão}
#   \item{media}{média da distribuição do efeito aleatório}
#   \item{sd}{desvio padrão da distribuição do efeito aleatório}
#   \item{x0}{valor inicial da série}
# }
#
# @return vetor com a série simulada
# 
# @examples
# simulacao_ar_1(list(n = 10, phi = 0.5))
#
simulacao_ar_1 <- function(parametros) {
  if (!is.list(parametros)) stop("Os parâmetros devem ser passados em uma lista, por exemplo: list(n = 10, phi = 0.5)")
  if (!all(c("n", "phi") %in% names(parametros))) stop("Os parâmetros 'n' e 'phi' são obrigatórios")
  
  n <- parametros$n + 1
  
  media_pop <- ifelse(!is.numeric(parametros$media), 0, parametros$media)
  sd_pop <- ifelse(!is.numeric(parametros$sd), 1, parametros$sd)
  
  serie <- numeric(n)
  ifelse(is.null(parametros$x1), serie[1] <- rnorm(1), serie[1] <- parametros$x0)
  for (i in 2:n) {
    serie[i] <- parametros$phi * serie[i - 1] + rnorm(1, mean = media_pop, sd = sd_pop)
  }
  
  return(serie[-1])
}

# Função para computar o EWMA
ewma <- function(dados, lambda, x0 = NULL) {
  n <- length(dados)
  final <- ifelse(is.null(x0), n, n + 1)
  ewma_serie <- numeric(final)
  ewma_serie[1] <- ifelse(is.null(x0), dados[1], x0)
  for (i in 2:final) {
    ewma_serie[i] <- lambda * dados[i] + (1 - lambda) * ewma_serie[i - 1]
  }
  return(ewma_serie[ifelse(is.null(x0), 1, 2):final])
}

# Função para computar resíduos do EWMA
ewma_residuos <- function(dados, ewma) {
  n <- length(dados)
  residuos <- numeric(n)
  residuos[1] <- dados[1]
  for (i in 2:n) {
    residuos[i] <- dados[i] - ewma[i - 1]
  }
  return(residuos)
}
