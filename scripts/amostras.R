source("scripts/arima_funcs.R")

PARAMETROS <- data.frame(
  phi = 0.2,
  media = 0,
  sd = 1,
  n_h0 = 100,
  n_h1 = 200
)

gerar_amostra_inicial <- function(){
  amostra_inicial <- ts(simulacao_ar_1(list(n = PARAMETROS$n_h0, phi = PARAMETROS$phi)))
  modelo <- arima(amostra_inicial, order = c(1, 0, 0), include.mean = FALSE)
  residuos <- modelo$residuals
  
  list(
    serie = amostra_inicial,
    x0 = amostra_inicial[1],
    xn = amostra_inicial[PARAMETROS$n_h0],
    modelo = modelo,
    residuos = residuos,
    phi = coef(modelo)[[1]],
    sd_res = sd(residuos),
    media_res = mean(residuos)
  )
}

amostra_inicial <- gerar_amostra_inicial()

phi_0_para_amostras_seguintes <- seq(0, 0.8, by = 0.1)
n_simulacoes <- length(phi_0_para_amostras_seguintes)

amostras_subsequentes <- lapply(
  phi_0_para_amostras_seguintes,
  function(phi_amostra) ts(simulacao_ar_1(list(n = PARAMETROS$n_h1, phi = phi_amostra, x0 = amostra_inicial$xn)))
)

amostras_df <- data.frame()

for (i in 1:n_simulacoes) {
  amostras_df <- rbind(
    amostras_df,
    data.frame(
      phi = phi_0_para_amostras_seguintes[i],
      serie = as.numeric(amostras_subsequentes[[i]]),
      residuos = as.numeric(Arima(amostras_subsequentes[[i]], model = amostra_inicial$modelo, include.mean = F)$residuals)
    )
  )
}
