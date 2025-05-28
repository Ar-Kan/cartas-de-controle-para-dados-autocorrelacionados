library(BTSR)
library(ggplot2)
library(ggformula)
library(dplyr)
library(data.table)
library(writexl)
library(patchwork)
library(cowplot)

DEBUG <- FALSE

USAR_BETA_ARMA <- TRUE # se TRUE, simula o modelo βARMA
USAR_PHI_RNORM <- TRUE # se TRUE, simula o phi com rnorm(1, mean = phi.est, sd = sd.phi)
CONTINUAR_COLAGEM_APOS_EXCEDER_N_INICIAL <- FALSE # se TRUE, continua colando após exceder o tamanho inicial (H1 com apenas novas amostras)

SALVAR_GRAFICO <- FALSE
SALVAR_DF <- FALSE

PHI_REAL <- 0.2

coeficientes <- function(alpha = 0, nu = 20, phi = PHI_REAL) {
  if (!USAR_BETA_ARMA) {
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


sim <- function(n, coefs, y.start = NULL, sd = NULL) {
  if (!USAR_BETA_ARMA) {
    # Simula uma série temporal AR de tamanho `n`
    if (is.null(sd)) {
      return(
        arima.sim(
          n = n,
          model = coefs
        )
      )
    }
    return(
      arima.sim(
        n = n,
        model = coefs,
        sd = sd
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
  if (!USAR_BETA_ARMA) {
    # Ajusta um modelo AR
    # NOTA: O método CSS-ML começa com uma estimação por mínimos quadrados condicionais (CSS)
    #       e depois refina com ML, sendo mais seguro e recomendável para simulações em massa
    return(
      arima(
        yt,
        order = c(1, 0, 0),
        include.mean = FALSE,
        method = "ML" #"CSS-ML"
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
    report = FALSE,
    info = TRUE
  )
}

bc <- c(1, 3, 5, 7, 9)
seq_colagens <- unlist(lapply(seq(0, 100, 10), function(x) bc + x))
seq_colagens <- c(seq_colagens, 135, 150, 165, 175, 185, 199)

if (USAR_BETA_ARMA) {
  desvios_do_phi <- c(0, 0.1, -0.3, -0.6)
} else {
  desvios_do_phi <- c(0, 0.1, -0.3, -1)
}
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
  z.estatistica = numeric(),
  fora_de_controle.z = logical()
)

tempo_inicial <- Sys.time()
print(paste0("Começando simulações para ", ifelse(USAR_BETA_ARMA, "βARMA", "AR"), "... ", tempo_inicial))

TAMANHO_MONTE_CARLO <- 300
TAMANHO_BOOTSTRAP <- 1000

for (k in 1:TAMANHO_MONTE_CARLO) {
  if (k %% 10 == 0) {
    porcento_Completo <- round(k / TAMANHO_MONTE_CARLO * 100, 2)
    print(paste0("Iteração Monte Carlo: ", k, " (", porcento_Completo, "% completo)"))
  }
  
  for (n_inicial in tamanho_amostras_iniciais) {
    amostra_inicial <- sim(n = n_inicial, coefs = coeficientes())
    
    # NOTA: Estimativa dos parâmetros da amostra inicial
    coef_inicial <- fit(
      yt = amostra_inicial,
      start = list(alpha = 0.1, nu = 20, phi = 0.1)
    )
    alpha.est <- ifelse(USAR_BETA_ARMA, coef_inicial$coef["alpha"][[1]], NA)
    nu.est <- ifelse(USAR_BETA_ARMA, coef_inicial$coef["nu"][[1]], NA)
    phi.est <- ifelse(USAR_BETA_ARMA, coef_inicial$coef["phi"][[1]], coef_inicial$coef[[1]])
    
    if (USAR_BETA_ARMA & phi.est > 0.6) {
      print("Phi maior que 0.6, pulando...")
      next
    }
    
    if (!USAR_BETA_ARMA) {
      # Variância teórica para o Phi estimado do modelo AR(1)
      var.phi <- (coef_inicial$sigma2 * (1 - phi.est^2)) / n_inicial
      sd.phi <- sqrt(var.phi)
      residuos.modelo <- residuals(coef_inicial)
    } else {
      # Variância estimada pela máxima verossimilhança
      info.inversa <- solve(coef_inicial$info.Matrix)
      var.phi <- info.inversa["phi", "phi"]
      sd.phi <- sqrt(var.phi)
      residuos.modelo <- coef_inicial$residuals
    }
    
    if (DEBUG)
      print(paste0("Monte Carlo (", k, ") - phi: ", phi.est, ", n: ", n_inicial))
    
    # Obs.: Limite para manter a estabilidade do processo
    # Def1: o modelo AR(1) é estacionário, sse, |Φ| < 1
    # NOTA 1: Como N(μ, σ²) ⊂ ℝ, precisamos garantir Def1
    # NOTA 2: com |Φ| > 0.9, o modelo AR(1) é instável, especialmente com amostras pequenas.
    # NOTA 3: βAR(1) é instável com |Φ| > 0.6.
    if (USAR_BETA_ARMA) {
      phi.maximo <- 0.6
    }
    else {
      phi.maximo <- ifelse(n_inicial < 50, 0.8, 1)
    }
    
    # Primeiro nível do Bootstrap Paramétrico
    # NOTA: reflete o comportamento do estimador de phi, assumindo que o processo não mudou
    phis.bootstrap <- replicate(TAMANHO_BOOTSTRAP,  {
      if (USAR_PHI_RNORM) {
        repeat {
          # Simula um novo phi com desvio teórico
          # NOTA 1: essa abordagem assume que o phi.est da série inicial é um ponto fixo confiável
          # NOTA 2: o rnorm(...) serve como proxy para o erro amostral do phi
          # NOTA 3: estamos dizendo: "Não conheço o phi verdadeiro, mas sei que ele pode ser qualquer
          #         valor em torno de phi.est, com desvio padrão sd.phi".
          phi.boot <- rnorm(1, mean = phi.est, sd = sd.phi)
          if (abs(phi.boot) < phi.maximo) break
        }
      } else {
        phi.boot <- phi.est
      }
      
      repeat {
        # Simula uma nova amostra com os parâmetros estimados, a partir da amostra inicial
        amostra <- sim(
          n = n_inicial,
          coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = phi.boot),
          sd = sqrt(coef_inicial$sigma2)
        )
        
        # Ajusta o modelo na nova amostra
        coef <- fit(
          yt = amostra,
          start = list(alpha = alpha.est, nu = nu.est, phi = phi.boot)
        )
        
        phi.novo <- ifelse(USAR_BETA_ARMA, coef$coef["phi"][[1]], coef$coef[[1]])
        if (abs(phi.boot) < phi.maximo) break
      }
      
      phi.novo
    })

    limite.inf <- quantile(phis.bootstrap, probs = 0.025)
    limite.sup <- quantile(phis.bootstrap, probs = 0.975)
    
    if (DEBUG)
      print(paste0("Limites inferiores e superiores: ", limite.inf, ", ", limite.sup))
    
    for (desvio in desvios_do_phi) {
      # Próximas amostras
      novo_phi <- PHI_REAL + desvio
      amostras_futuras <- sim(
        n = 200, # 200 é o maior tamanho de amostra neste estudo
        coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = novo_phi)
      )
      for (sc in seq_colagens) {
        proximas_amostras <- amostras_futuras[1:sc]
        if (sc >= n_inicial) {
          if (!CONTINUAR_COLAGEM_APOS_EXCEDER_N_INICIAL) {
            # Não realiza colagem após exceder o tamanho inicial
            
            break
          }
          # NOTA: Sem colagem, apenas amostras novas
          novo_dataset <- proximas_amostras[(sc - n_inicial + 1):sc]
        } else {
          # NOTA: Colagem de amostras antigas (amostra inicial) e novas
          novo_dataset <- c(
            amostra_inicial[(sc + 1):n_inicial],
            proximas_amostras
          )
        }
        
        stopifnot("A nova amostra não possui o tamanho esperado" = length(novo_dataset) == n_inicial)
        
        ok <- FALSE
        tryCatch({
          novo_coef <- fit(
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
        
        phi.nova.amostra <- ifelse(USAR_BETA_ARMA, novo_coef$coef["phi"][[1]], novo_coef$coef[[1]])
        
        if (!USAR_BETA_ARMA) {
          # Variância teórica para o Phi estimado do modelo AR(1)
          var.phi.nova.amostra <- (novo_coef$sigma2 * (1 - phi.nova.amostra^2)) / n_inicial
        } else {
          # Variância estimada pela máxima verossimilhança
          info.inversa <- solve(novo_coef$info.Matrix)
          var.phi.nova.amostra <- info.inversa["phi", "phi"]
        }
        
        # Z-score
        z.estatistica <- (phi.nova.amostra - phi.est) / sqrt(var.phi.nova.amostra + var.phi)
        limites.z.inf <- qt(0.025, df = n_inicial - 1)
        limites.z.sup <- qt(0.975, df = n_inicial - 1)
        
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
            z.estatistica = z.estatistica,
            fora_de_controle.z = (z.estatistica < limites.z.inf) | (z.estatistica > limites.z.sup)
          )
        )
      }
    }
  }
}

tempo_final <- Sys.time()
print(paste0("Fim das simulações para ", ifelse(USAR_BETA_ARMA, "BARMA", "AR"), "... ", tempo_final))
print(paste0("Tempo total (horas): ", round(difftime(tempo_final, tempo_inicial, units = "hours"), 2)))


df_resumo <- df %>%
  group_by(
    tamanho_inicial,
    tamanho_colagem,
    desvio
  ) %>%
  summarise(
    controle = sum(fora_de_controle),
    dp = abs(sd(phi_estimado)),
    controle.z = sum(fora_de_controle.z),
    total = n()
  ) %>%
  mutate(
    proporcao = controle / total,
    proporcao.z = controle.z / total,
    x_label = tamanho_colagem + tamanho_inicial,
  )

# Função para gerar o gráfico para cada tamanho_inicial
gerar_grafico <- function(tam_inicial_dado, break_by = 10) {
  data_breaks <- seq(tam_inicial_dado, tam_inicial_dado + 200, by = break_by)
  data_breaks <- c(
    tam_inicial_dado + 1,
    data_breaks[data_breaks > tam_inicial_dado + 1]
  )
  
  p <- ggplot(
      df_resumo %>% filter(tamanho_inicial == tam_inicial_dado),
      aes(x = x_label, y = proporcao, color = as.factor(desvio))
    ) +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    geom_line(linewidth = 1, alpha = 0.6) +
    geom_line(aes(y = proporcao.z),linewidth = 1, alpha = 0.6) +
    geom_line(aes(y = proporcao.z), linetype = "dashed", linewidth = 1, alpha = 0.6) +
    labs(
      title = paste("Tamanho inicial:", tam_inicial_dado),
      x = "Quantidade total de observações",
      y = "Proporção fora de controle",
      color = "Desvio do Φ"
    ) +
    scale_color_manual(values = c("red", "blue", "darkgreen", "orange", "purple", "brown")) +
    scale_x_continuous(breaks = data_breaks) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(p)
}

# Gerar gráficos individuais
g1 <- gerar_grafico(25, 5)
g2 <- gerar_grafico(50, 10)
g3 <- gerar_grafico(100, 25)
g4 <- gerar_grafico(200, 50)

# Montar layout com patchwork
layout <- (g1 | g2) / (g3 | g4)

# Adicionar título e legenda única
final <- layout +
  plot_annotation(
    title = paste0(
      "Proporção de controle fora dos limites - ",
      ifelse(USAR_BETA_ARMA, "β", ""),
      "AR(1)"
    ),
    subtitle = "Linha contínua: controle com bootstrap | Linha tracejada: controle com z-score",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  ) & theme(legend.position = "bottom")


if (!USAR_BETA_ARMA) {
  if (SALVAR_GRAFICO) {
    ggsave("controle_ar.png", final, width = 10, height = 6, units = "in", dpi = 300)
    print("Salvo gráfico AR(1)...")
    
    if (SALVAR_DF) {
      write_xlsx(df, "dados_simulacao_ar.xlsx")
      print("Salvos dados AR(1)...")
    }
  }
} else {
  if (SALVAR_GRAFICO) {
    ggsave("controle_bar.png", final, width = 10, height = 6, units = "in", dpi = 300)
    print("Salvo gráfico BARMA...")
    
    if (SALVAR_DF) {
      write_xlsx(df, "dados_simulacao_bar.xlsx")
      print("Salvos dados BARMA...")
    }
  }
}

print(final)
