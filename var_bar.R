library(BTSR)
library(MASS)

parametros <- c(alpha = 0, phi = 0.2, nu = 20)

conf <- c()

for (i in 1:100) {
  amostra <- BARFIMA.sim(
    n = 100,
    coefs = list(
      alpha = parametros["alpha"],
      nu = parametros["nu"],
      phi = parametros["phi"],
      d = 0
    ),
    error.scale = 1,
    complete = FALSE
  )
  
  fit <- BARFIMA.fit(
    yt = amostra,
    start = list(alpha = 0.1, nu = 20, phi = 0.1),
    p = 1, # Ordem do polinômio AR
    d = 0,
    error.scale = 1,
    report = FALSE,
    info = TRUE
  )
  
  estimativas <- c(
    alpha = fit$coefficients["alpha"][[1]],
    phi = fit$coefficients["phi"][[1]],
    nu = fit$coefficients["nu"][[1]]
  )
  
  info.inversa <- solve(fit$info.Matrix)
  
  # estimativas.novas <- mvrnorm(1, estimativas, info.inversa, empirical = FALSE)
  novo_phi <- rnorm(1, estimativas["phi"], sqrt(info.inversa["phi", "phi"]))
  
  ampl <- 1.96 * sqrt(info.inversa["phi", "phi"])
  sup <- parametros["phi"] + ampl
  inf <- parametros["phi"] - ampl
  # conf <- c(conf, estimativas.novas["phi"] > inf & estimativas.novas["phi"] < sup)
  conf <- c(conf, novo_phi > inf & novo_phi < sup)
  
  # Primeira estimativa de Φ
  # estimativas["phi"]
  
  # Nova estimativa de Φ
  # estimativas.novas["phi"]
  
  # Variância de Φ
  # info.inversa["phi", "phi"]
}

# Verifica a proporção de intervalos de confiança que contém o valor verdadeiro
proporcao <- sum(conf) / length(conf)
print(paste0(
  "A proporção de intervalos de confiança que contém o valor verdadeiro é: ",
  round(proporcao, 4)
))
