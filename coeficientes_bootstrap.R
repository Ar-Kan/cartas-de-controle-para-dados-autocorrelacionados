library(BTSR)
library(ggplot2)
library(ggformula)
library(dplyr)
library(data.table)
library(writexl)

DEBUG <- FALSE

USAR_BETA_ARMA <- TRUE
# USAR_PHI_RNORM: se TRUE, simula o phi com rnorm(1, mean = phi.est, sd = sd.phi)
USAR_PHI_RNORM <- TRUE

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

seq_colagens <- c(1, 24, 49, 75, 99, 199)
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
    
    # Primeiro nível do Bootstrap
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
        
        # residuo.aleatoriezado <- sample(residuos.modelo, size = n_inicial, replace = TRUE)
        # amostra <- numeric(n_inicial)
        # amostra[1] <- amostra_inicial[1]
        # 
        # for (t in 2:n_inicial) {
        #   amostra[t] <- phi.boot * amostra[t - 1] + residuo.aleatoriezado[t]
        # }
        
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
    
    # # Segundo nível (Bootstrap interno ou Monte Carlo futuro)
    # # NOTA: reflete como o phi se comportaria em novas amostras
    # phis.bootstrap.segundo <- numeric(TAMANHO_BOOTSTRAP)
    # for (i in 1:TAMANHO_BOOTSTRAP) {
    #   # Simula uma nova amostra com os parâmetros estimados
    #   amostra <- sim(
    #     n = n_inicial * 4,
    #     coefs = coeficientes(alpha = alpha.est, nu = nu.est, phi = phis.bootstrap[i]),
    #     sd = sqrt(coef_inicial$sigma2)
    #   )
    #   # Ajusta o modelo na nova amostra
    #   coef <- fit(
    #     yt = amostra,
    #     start = list(alpha = alpha.est, nu = nu.est, phi = phis.bootstrap[i])
    #   )
    # 
    #   phis.bootstrap.segundo[i] <- ifelse(USAR_BETA_ARMA, coef$coef["phi"][[1]], coef$coef[[1]])
    # }
    # 
    # # Quantis do bootstrap paramétrico
    # limite.inf <- quantile(phis.bootstrap.segundo, probs = 0.025)
    # limite.sup <- quantile(phis.bootstrap.segundo, probs = 0.975)
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


if (USAR_BETA_ARMA) {
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
    )
  
  p <- ggplot(df_resumo, aes(x = tamanho_colagem, y = proporcao, color = as.factor(desvio))) +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    geom_line(size = 1, alpha = 0.6) +
    geom_line(aes(y = proporcao.z), linetype = "dashed", linewidth = 1, alpha = 0.6) +
    labs(
      x = "Quantidade de novas observações",
      y = "Proporção fora de controle",
      color = "Desvio do Phi",
      title = "Proporção de controle fora dos limites - βAR(1)",
      subtitle = "Linha contínua: controle com bootstrap\nLinha tracejada: controle com z-score"
    ) +
    theme_minimal() +
    facet_wrap(~tamanho_inicial) +
    scale_color_manual(values = c("red", "blue", "darkgreen", "orange", "purple", "brown")) +
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

  if (SALVAR_GRAFICO) {
    ggsave("controle_bar.png", p, width = 10, height = 6, units = "in", dpi = 300)
    print("Salvando gráfico BARMA...")
    
    if (SALVAR_DF) {
      write_xlsx(df, "dados_simulacao_bar.xlsx")
      print("Salvando dados BARMA...")
    }
  }
}

if (!USAR_BETA_ARMA) {
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
    )
  
  p <- ggplot(df_resumo, aes(x = tamanho_colagem, y = proporcao, color = as.factor(desvio))) +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    geom_line(size = 1, alpha = 0.6) +
    geom_line(aes(y = proporcao.z), linetype = "dashed", linewidth = 1, alpha = 0.6) +
    labs(
      x = "Quantidade de novas observações",
      y = "Proporção fora de controle",
      color = "Desvio do Phi",
      title = "Proporção de controle fora dos limites - AR(1)",
      subtitle = "Linha contínua: controle com bootstrap\nLinha tracejada: controle com z-score"
    ) +
    theme_minimal() +
    facet_wrap(~tamanho_inicial) +
    scale_color_manual(values = c("red", "blue", "darkgreen", "orange", "purple", "brown")) +
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
  
  if (SALVAR_GRAFICO) {
    ggsave("controle_ar.png", p, width = 10, height = 6, units = "in", dpi = 300)
    print("Salvando gráfico AR(1)...")
    
    if (SALVAR_DF) {
      write_xlsx(df, "dados_simulacao_ar.xlsx")
      print("Salvando dados AR(1)...")
    }
  }
}

print(p)
