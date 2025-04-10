library(BTSR)
library(ggplot2)
library(ggformula)
library(dplyr)
library(data.table)
library(writexl)

BETA_ARMA <- FALSE
SALVAR_DADOS <- TRUE
SALVAR_DF <- TRUE

PHI_REAL <- 0.2

coeficientes <- function(alpha = 0, nu = 20, phi = PHI_REAL) {
  if (!BETA_ARMA) {
    # Modelo AR
    return(
      list(ar = c(phi))
    )
  }
  
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
    alpha = alpha,
    nu = nu,
    phi = phi,
    d = 0
  )
}


sim <- function(n, coefs, y.start = NULL) {
  if (!BETA_ARMA) {
    # Simula uma série temporal AR de tamanho `n`
    return(
      arima.sim(
        n = n,
        model = coefs
      )
    )
  }
  
  # Simula uma série temporal BAR de tamanho `n`
  BARFIMA.sim(
    n = n,
    coefs = coefs,
    y.start = y.start,
    error.scale = 1,
    complete = FALSE
  )
}

fit <- function(yt, start) {
  if (!BETA_ARMA) {
    # Ajusta um modelo AR
    return(
      arima(
        yt,
        order = c(1, 0, 0),
        include.mean = FALSE,
        method = "ML"
      )
    )
  }
  
  # Ajusta um modelo BAR
  BARFIMA.fit(
    yt = yt,
    start = start,
    p = 1, # Ordem do polinômio AR
    d = 0,
    error.scale = 1,
    report = FALSE
  )
}

seq_colagens <- c(1, 24, 49, 75, 99, 199)
desvios_do_phi <- c(0, 0.1, -0.3)
tamanho_amostras_iniciais <- c(25, 50, 100, 200)

df <- data.table(
  tamanho_inicial = numeric(),
  tamanho_colagem = numeric(),
  # Phi estimado para cada tamanho de colagem
  phi_estimado = numeric(),
  novo_phi_estimado = numeric(),
  desvio = numeric(),
  # SEQ dos phis para cada tamanho de colagem
  fora_de_controle = logical(),
  limite.inf = numeric(),
  limite.sup = numeric(),
  # SEQ dos phis para cada tamanho de colagem, no ARMA
  fora_de_controle.teorico = logical(),
  limite.inf.teorico = numeric(),
  limite.sup.teorico = numeric()
)

tempo_inicial <- Sys.time()
print(paste0("Começando simulações para ", ifelse(BETA_ARMA, "BARMA", "AR"), "... ", tempo_inicial))

TAMANHO <- 300
for (k in 1:TAMANHO) {
  for (n_inicial in tamanho_amostras_iniciais) {
    amostra_inicial <- sim(n = n_inicial, coefs = coeficientes())
    
    coef_inicial <- fit(
      yt = amostra_inicial,
      start = list(alpha = 0.1, nu = 20, phi = 0.1)
    )
    alpha.est <- ifelse(BETA_ARMA, coef_inicial$coef["alpha"][[1]], NA)
    nu.est <- ifelse(BETA_ARMA, coef_inicial$coef["nu"][[1]], NA)
    phi.est <- ifelse(BETA_ARMA, coef_inicial$coef["phi"][[1]], coef_inicial$coef[[1]])
    
    if (BETA_ARMA & phi.est > 0.6) {
      print("Phi maior que 0.6, pulando...")
      next
    }
    
    print(paste0("Bootstrap (", k, ") - phi: ", phi.est, ", n: ", n_inicial))
    bootstrap <- list()
    for (boots in 1:1000) {
      amostra <- sim(
        n = n_inicial,
        coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = phi.est)
      )
      coef <- fit(
        yt = amostra,
        start = list(alpha = alpha.est, nu = nu.est, phi = phi.est)
      )
      bootstrap[[boots]] <- list(
        amostra = amostra,
        coef = list(
          alpha = ifelse(BETA_ARMA, coef$coef["alpha"][[1]], NA),
          nu = ifelse(BETA_ARMA, coef$coef["nu"][[1]], NA),
          phi = ifelse(BETA_ARMA, coef$coef["phi"][[1]], coef$coef[[1]])
        )
      )
    }
    
    # Quartis do bootstrap paramétrico
    lista.phis <- unlist(lapply(bootstrap, function(x) x$coef$phi))
    limite.inf <- quantile(lista.phis, probs = 0.025)
    limite.sup <- quantile(lista.phis, probs = 0.975)
    
    if (!BETA_ARMA) {
      # Quantis teóricos para o Phi estimado do modelo AR(1)
      sd.teorico <- sqrt((sd(coef_inicial$residuals)^2 * (1 - phi.est^2)) / n_inicial)
      limite.inf.teorico <- qnorm(0.025, mean = phi.est, sd = sd.teorico)
      limite.sup.teorico <- qnorm(0.975, mean = phi.est, sd = sd.teorico)
    } else {
      limite.inf.teorico <- NA
      limite.sup.teorico <- NA
    }
    
    
    for (desvio in desvios_do_phi) {
      # Próximas amostras
      novo_phi <- PHI_REAL + desvio
      proximas_amostras_100 <- sim(
        n = 200,
        coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = novo_phi)
      )
      for (sc in seq_colagens) {
        proximas_amostras <- proximas_amostras_100[1:sc]
        if (sc >= n_inicial) {
          # sem colagem, apenas amostras novas
          novo_dataset <- proximas_amostras[max(c(1, sc - n_inicial)):sc]
        } else {
          novo_dataset <- c(
            amostra_inicial[(sc + 1):length(amostra_inicial)],
            proximas_amostras
          )
        }
        
        ok <- FALSE
        tryCatch({
          coef <- fit(
            yt = novo_dataset,
            start = list(alpha = alpha.est, nu = nu.est, phi = phi.est)
          )
          ok <- TRUE
        }, error = function(e) {
        })
        
        if (!ok) {
          print("Erro no ajuste do modelo, pulando...")
          next
        }
        
        phi.nova.amostra <- ifelse(BETA_ARMA, coef$coef["phi"][[1]], coef$coef[[1]])
        phi.nova.amostra.round <- round(phi.nova.amostra, 4)
        
        df <- rbind(
          df,
          data.table(
            tamanho_inicial = n_inicial,
            tamanho_colagem = sc,
            phi_estimado = phi.est,
            novo_phi_estimado = phi.nova.amostra,
            desvio = desvio,
            fora_de_controle = (phi.nova.amostra < limite.inf) | (phi.nova.amostra > limite.sup),
            limite.inf = limite.inf,
            limite.sup = limite.sup,
            fora_de_controle.teorico = (phi.nova.amostra < limite.inf.teorico) | (phi.nova.amostra > limite.sup.teorico),
            limite.inf.teorico = limite.inf.teorico,
            limite.sup.teorico = limite.sup.teorico
          )
        )
      }
    }
  }
}

tempo_final <- Sys.time()
print(paste0("Fim das simulações para ", ifelse(BETA_ARMA, "BARMA", "AR"), "... ", tempo_final))
print(paste0("Tempo total (horas): ", round(difftime(tempo_final, tempo_inicial, units = "hours"), 2)))


if (BETA_ARMA) {
  df_resumo <- df %>%
    group_by(
      tamanho_inicial,
      tamanho_colagem,
      desvio
    ) %>%
    summarise(
      controle = sum(fora_de_controle),
      dp = abs(sd(phi_estimado)),
      total = n()
    ) %>%
    mutate(
      proporcao = controle / total
    )
  
  p <- ggplot(df_resumo, aes(x = tamanho_colagem, y = proporcao, color = as.factor(desvio))) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_line() +
    geom_errorbar(
      aes(ymin = proporcao - dp, ymax = proporcao + dp),
      width = 0.2,
      position = position_dodge(0.5)
    ) +
    labs(
      x = "Quantidade de novas observações",
      y = "Proporção fora de controle",
      color = "Desvio do Phi",
      title = "Proporção de controle fora dos limites - βAR(1)"
    ) +
    theme_minimal() +
    facet_wrap(~tamanho_inicial)

  if (SALVAR_DADOS) {
    ggsave("controle_bar.png", p, width = 10, height = 6, units = "in", dpi = 300)
    
    if (SALVAR_DF) {
      write_xlsx(df, "dados_simulacao_bar.xlsx")
    }
  }
}

if (!BETA_ARMA) {
  df_resumo <- df %>%
    group_by(
      tamanho_inicial,
      tamanho_colagem,
      desvio
    ) %>%
    summarise(
      controle = sum(fora_de_controle),
      controle.teorico = sum(fora_de_controle.teorico),
      dp = abs(sd(phi_estimado)),
      total = n()
    ) %>%
    mutate(
      proporcao = controle / total,
      proporcao.teorico = controle.teorico / total
    )
  
  p <- ggplot(df_resumo, aes(x = tamanho_colagem, y = proporcao, color = as.factor(desvio))) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_line() +
    geom_line(aes(y = proporcao.teorico), linetype = "dashed") +
    labs(
      x = "Quantidade de novas observações",
      y = "Proporção fora de controle",
      color = "Desvio do Phi",
      title = "Proporção de controle fora dos limites - AR(1)"
    ) +
    theme_minimal() +
    facet_wrap(~tamanho_inicial) +
    scale_color_manual(values = c("red", "blue", "green")) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 200, 25)) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0.5, 'cm')
    )
  
  if (SALVAR_DADOS) {
    ggsave("controle_ar.png", p, width = 10, height = 6, units = "in", dpi = 300)
    
    if (SALVAR_DF) {
      write_xlsx(df, "dados_simulacao_ar.xlsx")
    }
  }
}

p
