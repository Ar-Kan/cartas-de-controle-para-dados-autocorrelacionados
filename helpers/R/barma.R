#
# Parâmetros
#

n1 <- 100
n2 <- 200
phi_parametro <- 0.2
alpha_teste <- 0.05
quantil_teste <- qnorm(1 - alpha_teste / 2)

#
# Funções auxiliares
#

coeficientes <- function(phi) {
  # βARMA: o modelo de Rocha e Cribari-Neto (2009, 2017)
  #        é obtido definindo `coefs$d = 0`
  #        e `d = FALSE` e `error.scale = 1` (escala preditiva)
  #
  # d: parâmetro de longa dependência (long memory), quando
  #    diferente de zero, produz um modelo ARFIMA.
  #
  # ν (nu): parâmetro de precisão, quanto maior o seu valor,
  #         menor a variância condicional (para μ_t fixo).
  #         O pacote BTSR define como padrão ν = 20.
  list(
    alpha = 0,
    nu = 20,
    phi = phi,
    d = 0
  )
}

barma.sim <- function(n,
                      phi,
                      seed = NULL,
                      y.start = NULL) {
  BARFIMA.sim(
    n = n,
    coefs = coeficientes(phi),
    y.start = y.start,
    error.scale = 1,
    complete = F,
    seed = ifelse(is.null(seed), sample(1:1000, 1), seed)
  )
}

barma.phi_estimado <- function(yt,
                               alpha = 0,
                               nu = 20,
                               phi = 0.1) {
  BARFIMA.fit(
    yt = yt,
    start = list(alpha = alpha, nu = nu, phi = phi),
    p = 1, # Ordem do polinômio AR
    d = 0,
    error.scale = 1,
    report = F
  )$coefficients["phi"][[1]]
}

barma.residuos <- function(yt, phi_estimado) {
  BARFIMA.extract(
    yt = yt,
    coefs = coeficientes(phi_estimado),
    llk = F,
    info = F,
    error.scale = 1
  )$error
}