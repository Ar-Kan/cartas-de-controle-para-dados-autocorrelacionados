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
CONTINUAR_COLAGEM_APOS_EXCEDER_N_INICIAL <- FALSE # se TRUE, continua colando após exceder o tamanho inicial (H1 com apenas novas amostras)

SALVAR_GRAFICO <- FALSE
SALVAR_DF <- FALSE

PHI_REAL <- 0.2

TAMANHO_MONTE_CARLO <- 300
TAMANHO_BOOTSTRAP <- 600

coeficientes <- function(alpha = 0, nu = 20, phi = PHI_REAL) {
  if (!USAR_BETA_ARMA) {
    # Modelo AR
    return(
      list(ar = phi)
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
    ar.parametros <- list(
      n = n,
      model = coefs
    )
    if (!is.null(sd)) {
      ar.parametros$sd <- sd
    }
    return(do.call(arima.sim, ar.parametros))
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
    report = F,
    info = TRUE
  )
}

estimativas.modelo <- function(modelo.ar) {
  alpha.est <- ifelse(USAR_BETA_ARMA, modelo.ar$coef["alpha"][[1]], NA)
  nu.est <- ifelse(USAR_BETA_ARMA, modelo.ar$coef["nu"][[1]], NA)
  phi.est <- ifelse(USAR_BETA_ARMA, modelo.ar$coef["phi"][[1]], modelo.ar$coef[[1]])

  if (!USAR_BETA_ARMA) {
    # Variância teórica para o Phi estimado do modelo AR(1)
    var.phi <- (modelo.ar$sigma2 * (1 - phi.est^2)) / n_inicial
    sd.phi <- sqrt(var.phi)
    residuos.modelo <- residuals(modelo.ar)
  } else {
    # Variância estimada pela máxima verossimilhança
    info.inversa <- solve(modelo.ar$info.Matrix)
    var.phi <- info.inversa["phi", "phi"]
    sd.phi <- sqrt(var.phi)
    residuos.modelo <- modelo.ar$residuals
  }

  list(
    alpha = alpha.est,
    nu = nu.est,
    phi = phi.est,
    var.phi = var.phi,
    sd.phi = sd.phi,
    residuos.modelo = residuos.modelo
  )
}

colar.series <- function(serie.h0, amostras_futuras, n.novos) {
  n_inicial <- length(serie.h0)
  proximas_amostras <- amostras_futuras[1:n.novos]
  if (n.novos >= n_inicial) {
    # NOTA: Sem colagem, apenas amostras novas
    novo_dataset <- proximas_amostras[(n.novos - n_inicial + 1):n.novos]
  } else {
    # NOTA: Colagem de amostras antigas (amostra inicial) e novas
    novo_dataset <- c(
      amostra_inicial[(n.novos + 1):n_inicial],
      proximas_amostras
    )
  }

  stopifnot("A nova amostra não possui o tamanho esperado" = length(novo_dataset) == n_inicial)

  novo_dataset
}

executar.bootstrap <- function(n_inicial, phi.est, sd.phi, sigma2 = NULL, usar.phi.rnorm = FALSE) {
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
  phis.bootstrap <- replicate(TAMANHO_BOOTSTRAP, {
    if (usar.phi.rnorm) {
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
        coefs = coeficientes(phi = phi.boot),
        sd = sqrt(sigma2)
      )

      # Ajusta o modelo na nova amostra
      coef <- fit(yt = amostra, start = coeficientes(phi = phi.boot))

      phi.novo <- ifelse(USAR_BETA_ARMA, coef$coef["phi"][[1]], coef$coef[[1]])
      if (abs(phi.boot) < phi.maximo) break
    }

    phi.novo
  })

  list(
    phis = phis.bootstrap,
    limite.inf = quantile(phis.bootstrap, probs = 0.025),
    limite.sup = quantile(phis.bootstrap, probs = 0.975)
  )
}

realizar.controle <- function(n.inicial, n.novos, phi, y.start = NULL, sd = NULL, serie.h0 = NULL, usar.phi.rnorm = FALSE) {
  # Simula a série de controle
  serie.controle <- sim(n.novos, coefs = coeficientes(phi = phi), y.start = y.start, sd = sd)

  if (is.null(serie.h0)) {
    # Se não foi fornecida uma série inicial, usa a série de controle como base
    serie.h0 <- serie.controle
  } else {
    # Caso contrário, cola a série de controle na série inicial
    serie.h0 <- colar.series(serie.h0, serie.controle, n.novos)
  }

  # Ajusta o modelo AR ou BAR na série de controle
  modelo.controle <- fit(yt = serie.h0, start = coeficientes(phi = phi))

  # Estima os parâmetros do modelo
  estimativas.controle <- estimativas.modelo(modelo.controle)
  # Realiza o bootstrap para obter os limites de controle
  limites.bootstrap <- executar.bootstrap(
    n_inicial = n.inicial,
    phi.est = estimativas.controle$phi,
    sd.phi = estimativas.controle$sd.phi,
    usar.phi.rnorm = usar.phi.rnorm
  )

  # Retorna as estimativas e limites de controle
  list(
    limite.inf = limites.bootstrap$limite.inf,
    limite.sup = limites.bootstrap$limite.sup,
    modelo = modelo.controle,
    estimativas = estimativas.controle
  )
}

# Simulação de controle para modelos AR(1) e βAR(1)

passo <- 20 # Passo para a sequência de colagens
bc <- c(1, 5) # c(1, 3, 5, 7, 9)
seq_colagens <- unlist(lapply(seq(0, 100, passo), function(x) bc + x))
seq_colagens <- c(seq_colagens, 135, 150, 165, 175, 185, 199)

if (USAR_BETA_ARMA) {
  desvios_do_phi <- c(0, 0.1, -0.3) #, -0.6)
} else {
  desvios_do_phi <- c(0, 0.1, -0.3, -1)
}
tamanho_amostras_iniciais <- c(100) # c(25, 50, 100, 200)

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
  fora_de_controle.z = logical(),
  fora_de_controle.online = logical(),
  limite.inf.online = numeric(),
  limite.sup.online = numeric()
)

tempo_inicial <- Sys.time()
print(paste0("Começando simulações para ", ifelse(USAR_BETA_ARMA, "βARMA", "AR"), "... ", tempo_inicial))

for (k in 1:TAMANHO_MONTE_CARLO) {
  if (k %% 10 == 0) {
    porcento_Completo <- round(k / TAMANHO_MONTE_CARLO * 100, 2)
    print(paste0("Iteração Monte Carlo: ", k, " (", porcento_Completo, "% completo)"))
  }

  for (n_inicial in tamanho_amostras_iniciais) {
    controle <- realizar.controle(
      n.inicial = n_inicial,
      n.novos = n_inicial,
      phi = PHI_REAL,
      usar.phi.rnorm = TRUE, # Usar rnorm para simular o erro amostral do phi
    )
    amostra_inicial <- controle$modelo$series
    # NOTA: Estimativa dos parâmetros da amostra inicial
    estimativas_iniciais <- controle$estimativas

    if (USAR_BETA_ARMA & estimativas_iniciais$phi > 0.6) {
      print("Phi maior que 0.6, pulando...")
      next
    }

    if (DEBUG) {
      print(paste0("Monte Carlo (", k, ") - phi: ", estimativas_iniciais$phi, ", n: ", n_inicial))
      print(paste0("Limites inferiores e superiores: ", controle$limite.inf, ", ", controle$limite.sup))
    }

    for (desvio in desvios_do_phi) {
      # Próximas amostras
      novo_phi <- PHI_REAL + desvio
      amostras_futuras <- sim(
        n = max(tamanho_amostras_iniciais),
        coefs = coeficientes(phi = novo_phi),
        y.start = amostra_inicial[n_inicial],
      )
      for (sc in seq_colagens) {
        if (sc >= n_inicial && !CONTINUAR_COLAGEM_APOS_EXCEDER_N_INICIAL) {
          # Não realiza colagem após exceder o tamanho inicial
          break
        }
        novo_dataset <- colar.series(
          serie.h0 = amostra_inicial,
          amostras_futuras = amostras_futuras,
          n.novos = sc
        )

        ok <- FALSE
        tryCatch({
          novo_modelo <- fit(
            yt = novo_dataset,
            start = coeficientes(phi = estimativas_iniciais$phi)
          )
          ok <- TRUE
        }, error = \() { })

        if (!ok) {
          print("Erro no ajuste do modelo, pulando...")
          next
        }

        estimativas_novas <- estimativas.modelo(novo_modelo)

        controle_novo <- realizar.controle(
          n.inicial = n_inicial,
          n.novos = sc,
          phi = estimativas_novas$phi,
          y.start = amostra_inicial[n_inicial],
          sd = estimativas_novas$sd.phi,
          serie.h0 = amostra_inicial
        )

        # Z-score
        z.estatistica <- (
          (estimativas_novas$phi - estimativas_iniciais$phi)
            / sqrt(estimativas_novas$var.phi + estimativas_iniciais$var.phi)
        )
        limites.z.inf <- qt(0.025, df = n_inicial - 1)
        limites.z.sup <- qt(0.975, df = n_inicial - 1)

        df <- rbind(
          df,
          data.table(
            tamanho_inicial = n_inicial,
            tamanho_colagem = sc,
            phi_estimado = estimativas_iniciais$phi,
            novo_phi_estimado = estimativas_novas$phi,
            desvio = desvio,
            fora_de_controle = (estimativas_novas$phi < controle$limite.inf) | (estimativas_novas$phi > controle$limite.sup),
            limite.inf = controle$limite.inf,
            limite.sup = controle$limite.sup,
            z.estatistica = z.estatistica,
            fora_de_controle.z = (z.estatistica < limites.z.inf) | (z.estatistica > limites.z.sup),
            fora_de_controle.online = (estimativas_novas$phi < controle_novo$limite.inf) | (estimativas_novas$phi > controle_novo$limite.sup),
            limite.inf.online = controle_novo$limite.inf,
            limite.sup.online = controle_novo$limite.sup
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
    desvio,
  ) %>%
  summarise(
    controle = sum(fora_de_controle),
    controle.z = sum(fora_de_controle.z),
    controle.online = sum(fora_de_controle.online),
    total = n()
  ) %>%
  mutate(
    proporcao = controle / total,
    proporcao.z = controle.z / total,
    proporcao.online = controle.online / total,
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
    # geom_line(linewidth = 1, alpha = 0.6, linetype = "solid") +
    geom_line(aes(y = proporcao.z), linetype = "dashed", alpha = 0.6) +
    geom_line(aes(y = proporcao.online), linetype = "solid", alpha = 0.6) +
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
# final <- layout +
final <- gerar_grafico(100, 25)+
  plot_annotation(
    title = paste0(
      "Proporção de controle fora dos limites - ",
      ifelse(USAR_BETA_ARMA, "β", ""),
      "AR(1)"
    ),
    # subtitle = "Linha tracejada: controle com bootstrap | Linha pontilhada: controle z | Linha sólida: controle online",
    subtitle = "Linha sólida: controle online | Linha tracejada: controle z",
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
    ggsave("controle_bar_online.png", final, width = 10, height = 6, units = "in", dpi = 300)
    print("Salvo gráfico BARMA...")

    if (SALVAR_DF) {
      write_xlsx(df, "dados_simulacao_bar.xlsx")
      print("Salvos dados BARMA...")
    }
  }
}

print(final)
