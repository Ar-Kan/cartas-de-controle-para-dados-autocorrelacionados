Sys.setlocale("LC_CTYPE", "Portuguese_Brazil.UTF-8")

# Precisa ter a seguinte estrutura:
# Exemplo:
# r"(
# #' Título da função
# minha_funcao <- function(param1, param2) {}
# )"
txt <- enc2utf8(
  r"(
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
#' Antes da transformação, os coeficientes são truncados usando `eps` para evitar
#' valores exatamente na fronteira, onde a derivada diverge.
#'
#' @examples
#' set.seed(42)
#'
#' coef <- c(ar1 = 0.2, ma1 = 0.5)
#' vc <- matrix(c(0.01, 0, 0, 0.01), 2, 2)
#'
#' bootstrap_coef_validos(coef, vc)
bootstrap_coef_validos <- function(coef, matriz_vcov, eps = 1e-6) {}
)"
)

topic <- roxygen2::roc_proc_text(roxygen2::rd_roclet(), txt)[[1]]

rd_text <- paste(format(topic), collapse = "\n")
writeLines(rd_text, ".doc-preview.Rd", useBytes = TRUE)

tools::Rd2HTML(".doc-preview.Rd", out = ".doc-preview.html")
# browseURL("preview.html")
