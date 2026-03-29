rm(list = ls())
graphics.off()
options(error = NULL)

library(ggplot2)
library(dplyr)
library(data.table)
library(MASS)
library(future.apply)
library(parallelly)
library(scales)
library(progressr)

set.seed(42)

######################################
# PARÂMETROS
PHI_REAL_LISTA <- c(0.2)
THETA_REAL_LISTA <- c(0.5)

TAMANHOS_AMOSTRAIS_INICIAIS <- c(100)

MC <- 500
B <- 500

DESVIOS_PHI <- c(0)
DESVIOS_THETA <- c(0)

SEQ_COLAGENS <- c(1, 50, 99, 150)

N_WORKERS <- max(1L, parallelly::availableCores() - 1L)

message(sprintf("Workers configurados: %d", N_WORKERS))

# Windows: multisession = processos separados
future::plan(future::multisession, workers = N_WORKERS)

# habilita barra de progresso
progressr::handlers(global = TRUE)
progressr::handlers("progress") # txtprogressbar

########################################
# FUNÇÕES AUXILIARES

# Simula ARMA(1,1) com burn-in de 200
simula_arma <- function(n, phi, theta) {
  obs <- arima.sim(
    n = 200 + n,
    model = list(ar = phi, ma = theta)
  )
  tail(obs, n)
}

# Ajusta ARMA(1,1)
fit_arma <- function(serie, phi = NULL, theta = NULL) {
  aviso <- NULL

  fit <- tryCatch(
    withCallingHandlers(
      arima(
        serie,
        order = c(1, 0, 1),
        include.mean = FALSE,
        # Usa uma parametrização alternativa para garantir estacionaridade dos
        # termos AR, e trata a invertibilidade dos MA depois da otimização (invertibility enforcement).
        # NOTA: Não funciona com o método `CSS` puro
        transform.pars = TRUE,
        # `CSS-ML`: usa CSS para dar um chute inicial e refina com verossimilhança
        method = "CSS-ML",
        init = c(ar1 = phi, ma1 = theta),
        optim.method = "BFGS",
        optim.control = list(maxit = 1000, reltol = 1e-10)
      ),
      warning = function(w) {
        aviso <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  list(
    fit = fit,
    warning = aviso,
    convergiu = is.null(aviso)
  )
}

# Testa estacionaridade (para φ)
ar_valido <- function(phi) {
  is.finite(phi) && all(Mod(polyroot(c(1, -phi))) > 1)
}

# Testa invertibilidade (para θ)
ma_valido <- function(theta) {
  is.finite(theta) && all(Mod(polyroot(c(1, theta))) > 1)
}

# Extrai coeficientes e matriz de covariância
arma_coef <- function(modelo) {
  coefs <- coef(modelo)
  vc <- vcov(modelo)
  list(
    coef = coefs,
    vcov = vc,
    vcov_inv = solve(vc)
  )
}

# Gera coeficientes ~ N(coef0$coef, coef0$vcov) e verifica se são válidos
bootstrap_coef_validos <- function(coef0, max_tentativas = 100L) {
  Sigma <- coef0$vcov

  # A matriz é válida?
  if (!is.matrix(Sigma) ||
    any(!is.finite(Sigma))) {
    return(NULL)
  }

  # É positiva definida?
  # As vezes o modelo pode não convergir ou a matriz de covariância pode ser mal comportada
  if (inherits(try(chol(Sigma), silent = TRUE), "try-error")) {
    return(NULL)
  }

  for (i in seq_len(max_tentativas)) {
    coef_star <- MASS::mvrnorm(1, mu = coef0$coef, Sigma = Sigma)

    if (is.null(coef_star)) {
      return(NULL)
    }

    phi_star <- unname(coef_star["ar1"])
    theta_star <- unname(coef_star["ma1"])

    if (ar_valido(phi_star) && ma_valido(theta_star)) {
      return(c(phi_star = phi_star, theta_star = theta_star))
    }
  }

  NULL
}

executa_um_mc <- function(
  mc,
  PHI_REAL,
  THETA_REAL) {

  resultados_mc <- list()
  idx_res <- 1L

  for (N_INICIAL in TAMANHOS_AMOSTRAIS_INICIAIS) {
    # Simula série original sob controle (Fase I)
    # Observações nas quais que se assume que o sistema está em controle
    serie0 <- tryCatch(
      simula_arma(N_INICIAL, PHI_REAL, THETA_REAL),
      error = function(e) NULL
    )
    if (is.null(serie0)) next

    ajuste0 <- fit_arma(serie0)
    if (is.null(ajuste0) || !ajuste0$convergiu) next
    fit0 <- ajuste0$fit

    coef0 <- tryCatch(arma_coef(fit0), error = function(e) NULL)
    if (is.null(coef0)) next

    for (desvio_phi in DESVIOS_PHI) {
      for (desvio_theta in DESVIOS_THETA) {
        # Altera os parâmetros reais para simular a série alternativa (Fase II)
        phi_alt <- PHI_REAL + desvio_phi
        theta_alt <- THETA_REAL + desvio_theta

        if (!ar_valido(phi_alt) || !ma_valido(theta_alt)) next

        # Simula série alternativa (Fase II)
        # Observações que foram obtidas opós série de controle da Fase I
        serie1 <- tryCatch(
          simula_arma(max(SEQ_COLAGENS), phi_alt, theta_alt),
          error = function(e) NULL
        )
        if (is.null(serie1)) next

        for (sc in SEQ_COLAGENS) {
          ## Verifica se o número de novas observações é maior que a série da Fase I
          # if (sc >= N_INICIAL) break

          # Colagem: mantém os últimos da Fase I e os primeiros da Fase II
          serie1_controle <- c(
            tail(serie0, N_INICIAL - sc),
            head(serie1, sc)
          )


          ajuste1 <- fit_arma(serie1_controle, coef0$coef["ar1"], coef0$coef["ma1"])
          if (is.null(ajuste1) || !ajuste1$convergiu) next
          fit1_controle <- ajuste1$fit

          coef1_controle <- tryCatch(coef(fit1_controle), error = function(e) NULL)
          if (is.null(coef1_controle)) next

          phi_hat <- unname(coef1_controle["ar1"])
          theta_hat <- unname(coef1_controle["ma1"])
          if (!is.finite(phi_hat) || !is.finite(theta_hat)) next

          # Bootstrap para obter distribuição de T² sob H0 (sistema em controle)
          # COMEÇO DO BOOTSTRAP

          phi_boot <- rep(NA_real_, B)
          theta_boot <- rep(NA_real_, B)

          for (b in seq_len(B)) {
            # Gera coeficientes e adiciona variabilidade
            coef_star <- bootstrap_coef_validos(coef0)
            if (is.null(coef_star)) next

            phi_star <- unname(coef_star["phi_star"])
            theta_star <- unname(coef_star["theta_star"])

            # Simula série Fase II*
            serie1_b <- tryCatch(
              simula_arma(N_INICIAL, phi_star, theta_star),
              error = function(e) NULL
            )
            if (is.null(serie1_b)) next

            # Colagem para série bootstrap: mantém os últimos de Fase I e os primeiros de Fase II*
            serie_colada_b <- c(
              tail(serie0, N_INICIAL - sc),
              head(serie1_b, sc)
            )

            ajuste_b <- fit_arma(serie_colada_b, coef0$coef["ar1"], coef0$coef["ma1"])
            if (is.null(ajuste_b) || !ajuste_b$convergiu) next
            fit_b <- ajuste_b$fit

            coef_b <- tryCatch(coef(fit_b), error = function(e) NULL)
            if (is.null(coef_b)) next

            phi_boot[b] <- unname(coef_b["ar1"])
            theta_boot[b] <- unname(coef_b["ma1"])
          }

          # FIM DO BOOTSTRAP

          validos <- is.finite(phi_boot) & is.finite(theta_boot)
          phi_boot <- phi_boot[validos]
          theta_boot <- theta_boot[validos]

          if (length(phi_boot) < 5L) next

          # Matriz da informação de Fisher da série de controle
          coef_vcov <- fit1_controle$var.coef

          coef_vcov_inv <- tryCatch(solve(coef_vcov), error = function(e) NULL)
          if (is.null(coef_vcov_inv)) next

          phi_m <- mean(phi_boot)
          theta_m <- mean(theta_boot)

          # Calcula T² para cada par (φ*, θ*) do bootstrap
          t2_boot <- vapply(seq_along(phi_boot), function(i) {
            diff_b <- c(phi_boot[i], theta_boot[i]) - c(phi_m, theta_m)
            as.numeric(t(diff_b) %*% coef_vcov_inv %*% diff_b)
          }, numeric(1))

          # Obtém quantis de T² para o intervalo de confiança
          limite_sup <- quantile(t2_boot, probs = 0.95, na.rm = TRUE, names = FALSE)

          diff_controle <- c(phi_hat, theta_hat) - c(phi_m, theta_m)
          t2_controle <- as.numeric(
            t(diff_controle) %*%
              coef_vcov_inv %*%
              diff_controle
          )

          esta_fora_de_controle <- t2_controle > limite_sup

          resultados_mc[[idx_res]] <- data.table(
            mc = mc,
            phi_real = PHI_REAL,
            theta_real = THETA_REAL,
            phi_alt = phi_alt,
            theta_alt = theta_alt,
            n0 = N_INICIAL,
            n1 = sc,
            desvio_phi = desvio_phi,
            desvio_theta = desvio_theta,
            fora_de_controle = esta_fora_de_controle,
            t2_sup = limite_sup,
            t2_controle = t2_controle,
            numero_de_bootstrap_validos = length(t2_boot)
          )

          idx_res <- idx_res + 1L
        }
      }
    }
  }

  if (length(resultados_mc) == 0L) {
    return(data.table())
  }

  rbindlist(resultados_mc, fill = TRUE)
}

executa_um_mc_safe <- function(mc, ..., p = NULL) {
  on.exit({
    if (!is.null(p)) p(sprintf("MC %d/%d", mc, MC))
  }, add = TRUE)

  tryCatch(
    executa_um_mc(mc = mc, ...),
    error = function(e) {
      data.table(
        mc = mc,
        erro_worker = as.character(e$message),
        pid = Sys.getpid()
      )
    }
  )
}

#####################################################################################
# LOOP PRINCIPAL

message(sprintf(
  "Iniciando simulações ARMA(1,1) com %d iterações Monte Carlo...",
  MC
))

DADOS_OUT_LIST <- list()
idx_global <- 1L

for (PHI_REAL in PHI_REAL_LISTA) {
  for (THETA_REAL in THETA_REAL_LISTA) {
    message(sprintf("Simulando para Φ = %.2f e Θ = %.2f", PHI_REAL, THETA_REAL))

    res_lista <- progressr::with_progress({
      p <- progressr::progressor(steps = MC)

      future_lapply(
        X = seq_len(MC),
        FUN = executa_um_mc_safe,
        future.seed = TRUE,
        PHI_REAL = PHI_REAL,
        THETA_REAL = THETA_REAL,
        p = p
      )
    })

    erros <- Filter(function(x) is.data.table(x) && "erro_worker" %in% names(x), res_lista)
    if (length(erros) > 0L) {
      message("Ocorreram erros em alguns workers:")
      print(rbindlist(erros, fill = TRUE))
    }

    validos <- Filter(function(x) is.data.table(x) && !("erro_worker" %in% names(x)), res_lista)

    if (length(validos) > 0L) {
      bloco <- rbindlist(validos, fill = TRUE)
      if (nrow(bloco) > 0L) {
        DADOS_OUT_LIST[[idx_global]] <- bloco
        idx_global <- idx_global + 1L
      }
    }
  }
}


###############################################################
# Resultados

if (length(DADOS_OUT_LIST) == 0L) {
  stop("Nenhum resultado válido foi gerado.")
}

DADOS_OUT <- rbindlist(DADOS_OUT_LIST, fill = TRUE)

message("Concluído com sucesso.")
print(head(DADOS_OUT))


DADOS_OUT %>%
  mutate(
    label = paste0(
      "(", phi_alt, ", ", theta_alt, ")"
    )
  ) %>%
  group_by(label, n0, n1) %>%
  summarise(proporcao = mean(fora_de_controle), .groups = "drop") %>%
  mutate(
    x_label = n1 + n0,
    se = sqrt(proporcao * (1 - proporcao) / MC),
    ic_inf = pmax(0, proporcao - 1.96 * se),
    ic_sup = pmin(1, proporcao + 1.96 * se)
  ) %>%
  ggplot(aes(x = x_label, y = proporcao, color = label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 5, alpha = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  labs(
    title = "Monte Carlo para ARMA(1,1)",
    x = "Número de novas observações",
    y = "Proporção Fora de Controle",
    color = "Parâmetros (Φ; Θ)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_brewer(
    palette = "Dark2",
    labels = scales::label_wrap(20),
    guide = guide_legend(nrow = 2)
  ) +
  facet_wrap(~n0,
    labeller = labeller(n0 = \(x) paste0("Tamanho Inicial: ", x)),
    scales = "free_x"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
