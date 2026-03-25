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

set.seed(42)

# ==============================
# PARÂMETROS
# ==============================
PHI_REAL_LISTA <- c(0.3)
THETA_REAL_LISTA <- c(0.6)

tamanhos_amostrais_iniciais <- c(100)

MC <- 500
B <- 500

DESVIOS_PHI <- c(0)
DESVIOS_THETA <- c(0)

SEQ_COLAGENS <- c(1, 25, 50, 75, 99, 200, 300, 400, 499)

N_WORKERS <- max(1L, parallelly::availableCores() - 1L)

message(sprintf("Workers configurados: %d", N_WORKERS))

# Windows: multisession = processos separados
future::plan(future::multisession, workers = N_WORKERS)

# ==============================
# FUNÇÕES AUXILIARES
# ==============================

# Simula ARMA(1,1)
simula_arma <- function(n, phi, theta, sd = NULL) {
  args <- list(
    n = 200 + n,
    model = list(ar = phi, ma = theta)
  )
  if (!is.null(sd)) args$sd <- sd
  obs <- do.call(arima.sim, args)
  tail(obs, n)
}

# Ajusta ARMA(1,1)
fit_arma <- function(serie) {
  tryCatch(
    arima(
      serie,
      order = c(1, 0, 1),
      include.mean = FALSE,
      method = "ML"
    ),
    error = function(e) NA
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

    phi_star <- unname(coef_star[1])
    theta_star <- unname(coef_star[2])

    if (ar_valido(phi_star) && ma_valido(theta_star)) {
      return(c(phi_star = phi_star, theta_star = theta_star))
    }
  }

  NULL
}

executa_um_mc <- function(
    mc,
    PHI_REAL,
    THETA_REAL,
    tamanhos_amostrais_iniciais,
    DESVIOS_PHI,
    DESVIOS_THETA,
    SEQ_COLAGENS,
    B,
    MC_TOTAL) {
  # cat(sprintf("mc=%d | pid=%d\n", mc, Sys.getpid()))

  resultados_mc <- list()
  idx_res <- 1L

  for (N_INICIAL in tamanhos_amostrais_iniciais) {
    # Simula série original sob controle (H0)
    serie0 <- tryCatch(
      simula_arma(N_INICIAL, PHI_REAL, THETA_REAL),
      error = function(e) NULL
    )
    if (is.null(serie0)) next

    fit0 <- fit_arma(serie0)
    if (length(fit0) == 1 && is.na(fit0)) next

    coef0 <- tryCatch(arma_coef(fit0), error = function(e) NULL)
    if (is.null(coef0)) next

    for (desvio_phi in DESVIOS_PHI) {
      for (desvio_theta in DESVIOS_THETA) {
        # Altera os parâmetros reais para simular a série alternativa (H1)
        phi_alt <- PHI_REAL + desvio_phi
        theta_alt <- THETA_REAL + desvio_theta

        if (!ar_valido(phi_alt) || !ma_valido(theta_alt)) next

        # Simula série alternativa (H1)
        serie1 <- tryCatch(
          simula_arma(N_INICIAL, phi_alt, theta_alt),
          error = function(e) NULL
        )
        if (is.null(serie1)) next

        for (sc in SEQ_COLAGENS) {
          if (sc >= N_INICIAL) break

          # Colagem: mantém os últimos de H0 e os primeiros de H1
          serie1_controle <- c(
            tail(serie0, N_INICIAL - sc),
            head(serie1, sc)
          )

          fit_tmp <- fit_arma(serie1_controle)
          if (length(fit_tmp) == 1 && is.na(fit_tmp)) next

          coef_tmp <- tryCatch(coef(fit_tmp), error = function(e) NULL)
          if (is.null(coef_tmp)) next

          phi_hat <- unname(coef_tmp["ar1"])
          theta_hat <- unname(coef_tmp["ma1"])

          if (!is.finite(phi_hat) || !is.finite(theta_hat)) next

          # Bootstrap para obter distribuição de T² sob H0
          # COMEÇO DO BOOTSTRAP

          phi_boot <- rep(NA_real_, B)
          theta_boot <- rep(NA_real_, B)

          for (b in seq_len(B)) {
            # Gera coeficientes e adiciona variabilidade
            coef_star <- bootstrap_coef_validos(coef0)
            if (is.null(coef_star)) next

            phi_star <- unname(coef_star["phi_star"])
            theta_star <- unname(coef_star["theta_star"])

            # Simula série bootstrap com os coeficientes gerados (H1*)
            serie_b <- tryCatch(
              simula_arma(N_INICIAL, phi_star, theta_star),
              error = function(e) NULL
            )
            if (is.null(serie_b)) next

            # Colagem para série bootstrap: mantém os últimos de H0 e os primeiros de H1*
            serie_colada <- c(
              tail(serie0, N_INICIAL - sc),
              head(serie_b, sc)
            )

            fit_b <- fit_arma(serie_colada)
            if (length(fit_b) == 1 && is.na(fit_b)) next

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

          # Calcula matriz de covariância dos coeficientes bootstrap
          coef_vcov <- matrix(
            c(
              var(phi_boot), cov(phi_boot, theta_boot),
              cov(theta_boot, phi_boot), var(theta_boot)
            ),
            nrow = 2,
            byrow = TRUE
          )

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
          quant_boot <- quantile(
            t2_boot,
            probs = c(0.025, 0.975),
            na.rm = TRUE,
            names = FALSE
          )

          # Calcula T² para o par (φ^, θ^) da série colada (H1)
          diff_ctrl <- c(phi_hat, theta_hat) - c(phi_m, theta_m)
          t2_controle <- as.numeric(
            t(diff_ctrl) %*% coef_vcov_inv %*% diff_ctrl
          )

          esta_fora_de_controle <- t2_controle < quant_boot[1] || t2_controle > quant_boot[2]

          resultados_mc[[idx_res]] <- data.table(
            mc = mc,
            phi = phi_alt,
            theta = theta_alt,
            n0 = N_INICIAL,
            n1 = sc,
            desvio_phi = desvio_phi,
            desvio_theta = desvio_theta,
            fora_de_controle = esta_fora_de_controle,
            t2_inf = quant_boot[1],
            t2_sup = quant_boot[2],
            t2_controle = t2_controle,
            numero_de_bootstrap_validos = length(t2_boot)
          )

          idx_res <- idx_res + 1L
        }
      }
    }
  }

  # cat(sprintf(
  #   "Monte Carlo: %d/%d (%.2f%%) | Φ0=%.2f | Θ0=%.2f | pid=%d\n",
  #   mc, MC_TOTAL, 100 * mc / MC_TOTAL, PHI_REAL, THETA_REAL, Sys.getpid()
  # ))

  if (length(resultados_mc) == 0L) {
    return(data.table())
  }

  rbindlist(resultados_mc, fill = TRUE)
}

executa_um_mc_safe <- function(mc, ...) {
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

# ==============================
# LOOP PRINCIPAL
# ==============================

message(sprintf(
  "Iniciando simulações ARMA(1,1) com %d iterações Monte Carlo...",
  MC
))

DADOS_OUT_LIST <- list()
idx_global <- 1L

for (PHI_REAL in PHI_REAL_LISTA) {
  for (THETA_REAL in THETA_REAL_LISTA) {
    message(sprintf("Simulando para Φ0 = %.2f e Θ0 = %.2f", PHI_REAL, THETA_REAL))

    res_lista <- future_lapply(
      X = seq_len(MC),
      FUN = executa_um_mc_safe,
      future.seed = TRUE,
      PHI_REAL = PHI_REAL,
      THETA_REAL = THETA_REAL,
      tamanhos_amostrais_iniciais = tamanhos_amostrais_iniciais,
      DESVIOS_PHI = DESVIOS_PHI,
      DESVIOS_THETA = DESVIOS_THETA,
      SEQ_COLAGENS = SEQ_COLAGENS,
      B = B,
      MC_TOTAL = MC
    )

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

if (length(DADOS_OUT_LIST) == 0L) {
  stop("Nenhum resultado válido foi gerado.")
}

DADOS_OUT <- rbindlist(DADOS_OUT_LIST, fill = TRUE)

message("Concluído com sucesso.")
print(head(DADOS_OUT))
