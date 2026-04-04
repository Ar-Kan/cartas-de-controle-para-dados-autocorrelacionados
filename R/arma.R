#' Simula uma série ARMA(1,1) com burn-in
#'
#' Gera uma série ARMA(1,1) via [stats::arima.sim()] e descarta um período
#' inicial de burn-in para reduzir o efeito das condições iniciais.
#'
#' @param n Número de observações desejadas na série retornada.
#' @param phi Coeficiente autorregressivo AR(1).
#' @param theta Coeficiente de médias móveis MA(1).
#'
#' @return Vetor numérico de comprimento `n`.
#'
#' @details
#' A função simula `200 + n` observações e retorna apenas as `n` últimas.
#'
#' @examples
#' set.seed(42)
#' x <- simula_arma(n = 100, phi = 0.2, theta = 0.5)
#' length(x)
simula_arma <- function(n, phi, theta) {
  obs <- arima.sim(
    n = 200 + n,
    model = list(ar = phi, ma = theta)
  )
  # Retorna as últimas n observações
  tail(obs, n)
}

#' Ajusta um modelo ARMA(1,1) sem média
#'
#' Ajusta um modelo ARMA(1,1) usando [stats::arima()] com método `"CSS-ML"`.
#' Warnings do ajuste são capturados e retornados, em vez de serem emitidos
#' diretamente na console.
#'
#' @param serie Vetor numérico ou série temporal.
#' @param phi Valor inicial opcional para o coeficiente AR(1).
#' @param theta Valor inicial opcional para o coeficiente MA(1).
#' @param transform.pars Lógico; repassado para o argumento homônimo de
#'   [stats::arima()].
#'
#' @return Lista com os elementos:
#' \describe{
#'   \item{fit}{Objeto retornado por [stats::arima()], ou `NULL` em caso de falha.}
#'   \item{warnings}{Vetor de caracteres com warnings capturados.}
#'   \item{convergiu}{`TRUE` quando não houve warnings capturados.}
#'   \item{erro}{`TRUE` quando houve falha no ajuste ou no pós-processamento.}
#' }
#'
#' @details
#' A função:
#' \itemize{
#'   \item ajusta `order = c(1, 0, 1)`
#'   \item usa `include.mean = FALSE`
#'   \item usa `method = "CSS-ML"`
#'   \item usa `optim.method = "BFGS"`
#' }
#'
#' Se o ajuste falhar, ou se `coef()`/`vcov()` falharem, o retorno marca
#' `erro = TRUE` e `convergiu = FALSE`.
#'
#' @examples
#' set.seed(42)
#' serie <- simula_arma(200, 0.2, 0.5)
#' ajuste <- fit_arma(serie)
#' ajuste$erro
fit_arma <- function(serie, phi = NULL, theta = NULL, transform.pars = FALSE) {
  avisos <- character()

  fit <- withCallingHandlers(
    arima(
      serie,
      order = c(1, 0, 1),
      include.mean = FALSE,
      transform.pars = transform.pars,
      method = "CSS-ML",
      init = c(ar1 = phi, ma1 = theta),
      optim.method = "BFGS",
      optim.control = list(maxit = 1000, reltol = 1e-10)
    ),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ) |> tryNull()

  if (
    is.null(fit) ||
      coef(fit) |> tryNull() |> is.null() ||
      vcov(fit) |> tryNull() |> is.null()
  ) {
    return(list(
      fit = NULL,
      warnings = avisos,
      convergiu = FALSE,
      erro = TRUE
    ))
  }

  list(
    fit = fit,
    warnings = avisos,
    convergiu = length(avisos) == 0,
    erro = FALSE
  )
}


#' Verifica se um coeficiente AR(1) é estacionário
#'
#' @param phi Valor candidato para o coeficiente AR(1).
#'
#' @return `TRUE` se o parâmetro for finito e corresponder a um modelo
#'   estacionário; `FALSE` caso contrário.
#'
#' @details
#' Para AR(1), a condição de estacionaridade equivale a exigir que as raízes
#' do polinômio autorregressivo estejam fora do círculo unitário.
#'
#' @examples
#' ar_valido(0.5)
#' ar_valido(1.2)
ar_valido <- function(phi) {
  is.finite(phi) && all(Mod(polyroot(c(1, -phi))) > 1)
}


#' Verifica se um coeficiente MA(1) é invertível
#'
#' @param theta Valor candidato para o coeficiente MA(1).
#'
#' @return `TRUE` se o parâmetro for finito e corresponder a um modelo
#'   invertível; `FALSE` caso contrário.
#'
#' @details
#' Para MA(1), a invertibilidade equivale a exigir que as raízes do polinômio
#' de médias móveis estejam fora do círculo unitário.
#'
#' @examples
#' ma_valido(0.5)
#' ma_valido(1.2)
ma_valido <- function(theta) {
  is.finite(theta) && all(Mod(polyroot(c(1, theta))) > 1)
}

#' Combina uma série de Fase I com novas observações da Fase II
#'
#' Constrói uma nova série mantendo as observações mais recentes da Fase I e
#' substituindo as últimas posições por observações iniciais da Fase II.
#'
#' @param serie_fase1 Vetor numérico da Fase I.
#' @param serie_fase2 Vetor numérico da Fase II.
#' @param numero_de_novas_observacoes Número de observações da Fase II a serem
#'   incorporadas.
#'
#' @return Vetor numérico com o mesmo comprimento de `serie_fase1`.
#'
#' @details
#' Se `numero_de_novas_observacoes >= length(serie_fase1)`, a função retorna
#' apenas as primeiras `length(serie_fase1)` observações de `serie_fase2`.
#'
#' @examples
#' cola_series(
#'   serie_fase1 = 1:5,
#'   serie_fase2 = 101:110,
#'   numero_de_novas_observacoes = 2
#' )
cola_series <- function(serie_fase1, serie_fase2, numero_de_novas_observacoes) {
  n_inicial <- length(serie_fase1)

  if (numero_de_novas_observacoes >= n_inicial) {
    # Sem colagem
    out <- serie_fase2[seq_len(n_inicial)]
  } else {
    out <- c(
      tail(serie_fase1, n_inicial - numero_de_novas_observacoes),
      head(serie_fase2, numero_de_novas_observacoes)
    )
  }

  stopifnot("Série resultante tem comprimento diferente do esperado" = length(out) == n_inicial)
  out
}

#' Gera coeficientes bootstrap válidos para ARMA(1,1)
#'
#' Gera um novo par de coeficientes AR(1) e MA(1) a partir de uma aproximação
#' normal multivariada, utilizando uma transformação que garante que os valores
#' finais permaneçam no intervalo admissível `(-1, 1)`.
#'
#' @param coef Vetor nomeado contendo ao menos `ar1` e `ma1`.
#' @param matriz_vcov Matriz de covariância associada aos coeficientes.
#' @param eps Pequeno valor positivo utilizado para evitar avaliação exatamente
#'   na fronteira de `-1` e `1`.
#'
#' @return
#' Vetor nomeado com os elementos `phi_star` e `theta_star`, ou `NULL` quando:
#' \itemize{
#'   \item a matriz de covariância for inválida
#'   \item os coeficientes forem não finitos
#'   \item a matriz transformada não for positiva-definida
#' }
#'
#' @details
#' O método consiste em transformar os coeficientes para um espaço irrestrito,
#' realizar o sorteio nesse espaço, e então retornar ao espaço original.
#'
#' Seja o vetor de parâmetros:
#' \deqn{
#'   \boldsymbol{\beta} = (\phi, \theta)
#' }
#'
#' A transformação é aplicada componente a componente:
#' \deqn{
#'   u_{\phi} = \operatorname{atanh}(\phi), \quad
#'   u_{\theta} = \operatorname{atanh}(\theta)
#' }
#'
#' Em forma vetorial:
#' \deqn{
#'   \mathbf{u} = g(\boldsymbol{\beta})
#' }
#'
#' onde \eqn{g} atua separadamente em cada componente.
#'
#' A função \eqn{\operatorname{atanh}} mapeia o intervalo \eqn{(-1, 1)} em
#' \eqn{\mathbb{R}}, permitindo modelar os parâmetros com distribuição normal.
#'
#' \strong{Aproximação da covariância (método delta)}
#'
#' Seja \eqn{J} a jacobiana da transformação \eqn{g} avaliada em
#' \eqn{\boldsymbol{\beta}}. Então:
#' \deqn{
#'   \operatorname{Var}(\mathbf{u}) \approx
#'   J \, \operatorname{Var}(\boldsymbol{\beta}) \, J'
#' }
#'
#' Como a transformação é separável por componente, a jacobiana é diagonal:
#' \deqn{
#'   J =
#'   \begin{pmatrix}
#'     \frac{1}{1 - \phi^2} & 0 \\
#'     0 & \frac{1}{1 - \theta^2}
#'   \end{pmatrix}
#' }
#'
#' \strong{Geração bootstrap}
#'
#' No espaço transformado:
#' \deqn{
#'   \mathbf{u}^* \sim \mathcal{N}(\boldsymbol{\mu}_u, \Sigma_u)
#' }
#'
#' onde:
#' \itemize{
#'   \item \eqn{\boldsymbol{\mu}_u = g(\boldsymbol{\beta})}
#'   \item \eqn{\Sigma_u = J \, \operatorname{Var}(\boldsymbol{\beta}) \, J'}
#' }
#'
#' \strong{Transformação inversa}
#'
#' Os coeficientes são retornados ao espaço original via:
#' \deqn{
#'   \phi^* = \tanh(u_{\phi}^*), \quad
#'   \theta^* = \tanh(u_{\theta}^*)
#' }
#'
#' garantindo que ambos permaneçam no intervalo `(-1, 1)`.
#'
#' Antes da transformação e depois de aplicar a função inversa, há uma etapa de truncamento
#' para garantir que os valores permaneçam estritamente dentro do intervalo `(-1, 1)`,
#' evitando problemas numéricos associados a valores exatamente na fronteira.
#'
#' @examples
#' set.seed(42)
#'
#' coef <- c(ar1 = 0.2, ma1 = 0.5)
#' vc <- matrix(c(0.01, 0, 0, 0.01), 2, 2)
#'
#' bootstrap_coef_validos(coef, vc)
amostrar_coef_validos <- function(coef, matriz_vcov, eps = 1e-6) {
  # Validações
  if (!is.matrix(matriz_vcov) || any(!is.finite(matriz_vcov))) {
    return(NULL)
  }

  phi <- unname(coef["ar1"])
  theta <- unname(coef["ma1"])

  if (!is.finite(phi) || !is.finite(theta)) {
    return(NULL)
  }

  # proteção contra valores exatamente na fronteira
  phi <- max(min(phi, 1 - eps), -1 + eps)
  theta <- max(min(theta, 1 - eps), -1 + eps)

  # Média no espaço transformado
  mu_u <- c(
    ar1 = atanh(phi),
    ma1 = atanh(theta)
  )

  # Jacobiana da transformação g no ponto \beta
  J <- diag(c(
    1 / (1 - phi^2),
    1 / (1 - theta^2)
  ))

  Sigma_u <- J %*% matriz_vcov %*% t(J)
  Sigma_u <- as.matrix(Sigma_u)

  rownames(Sigma_u) <- colnames(Sigma_u) <- c("ar1", "ma1")

  # Valida matriz
  if (any(!is.finite(Sigma_u))) {
    return(NULL)
  }

  # Verifica se a matriz é positiva-definida para evitar erros no `mvrnorm`
  if (is.null(chol(Sigma_u) |> tryNull())) {
    return(NULL)
  }

  # Sorteia valores no espaço transformado
  # U* ~ N(mu_u, Sigma_u)
  u_star <- MASS::mvrnorm(1, mu = mu_u, Sigma = Sigma_u)

  if (is.null(u_star) || any(!is.finite(u_star))) {
    return(NULL)
  }

  # Volta ao espaço original
  #   phi* = tanh(u_phi*)
  #   theta* = tanh(u_theta*)
  phi_star <- tanh(unname(u_star["ar1"]))
  theta_star <- tanh(unname(u_star["ma1"]))

  # Proteção final
  if (!is.finite(phi_star) || !is.finite(theta_star)) {
    return(NULL)
  }

  # Proteção final contra saturação numérica em +/-1
  phi_star <- max(min(phi_star, 1 - eps), -1 + eps)
  theta_star <- max(min(theta_star, 1 - eps), -1 + eps)

  c(phi_star = phi_star, theta_star = theta_star)
}
