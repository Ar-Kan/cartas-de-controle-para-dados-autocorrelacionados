# mixture_lm_em.R
# EM algorithm for a mixture of two linear regressions

# --- 1. Generate synthetic data ----
set.seed(123)
n <- 300
p <- 2
X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n)) # intercept + 1 covariate
# True parameters
beta1_true <- c(2, 1.5)
beta2_true <- c(-1, 2.5)
sigma1_true <- 0.5
sigma2_true <- 1.0
pi_true <- 0.6

# Simulate latent component labels
z <- rbinom(n, size = 1, prob = pi_true) + 1
y <- numeric(n)
for (i in 1:n) {
  if (z[i] == 1) {
    y[i] <- X[i, ] %*% beta1_true + rnorm(1, 0, sigma1_true)
  } else {
    y[i] <- X[i, ] %*% beta2_true + rnorm(1, 0, sigma2_true)
  }
}

# --- 2. EM initialization ----
# Random starts
pi_k <- 0.5
beta1 <- lm(y ~ X[, -1])$coef + rnorm(p, 0, 0.1)
beta2 <- lm(y ~ X[, -1])$coef + rnorm(p, 0, 0.1)
sigma1 <- sd(y)
sigma2 <- sd(y)

max_iter <- 1000
tol <- 1e-6
loglik_prev <- -Inf

# --- 3. EM loop ----
for (iter in 1:max_iter) {
  # E‑step: compute responsibilities gamma[i,1]
  dens1 <- dnorm(y, mean = X %*% beta1, sd = sigma1)
  dens2 <- dnorm(y, mean = X %*% beta2, sd = sigma2)
  w1 <- pi_k * dens1
  w2 <- (1 - pi_k) * dens2
  gamma1 <- w1 / (w1 + w2)
  gamma2 <- 1 - gamma1

  # M‑step: update parameters
  # mixing proportion
  pi_k <- mean(gamma1)
  # weighted least squares for beta1 and beta2
  W1 <- diag(gamma1)
  W2 <- diag(gamma2)
  beta1 <- solve(t(X) %*% W1 %*% X, t(X) %*% W1 %*% y)
  beta2 <- solve(t(X) %*% W2 %*% X, t(X) %*% W2 %*% y)
  # variances
  sigma1 <- sqrt(sum(gamma1 * (y - X %*% beta1)^2) / sum(gamma1))
  sigma2 <- sqrt(sum(gamma2 * (y - X %*% beta2)^2) / sum(gamma2))

  # log‑likelihood for convergence check
  loglik <- sum(log(pi_k * dens1 + (1 - pi_k) * dens2))
  if (abs(loglik - loglik_prev) < tol) {
    message("Converged at iter = ", iter)
    break
  }
  loglik_prev <- loglik
}

# --- 4. Results ----
cat("Estimated mixing proportion π:\n", round(pi_k, 3), "\n\n")
cat("Estimated beta1:\n")
print(round(beta1, 3))
cat("\n")
cat("Estimated sigma1:\n", round(sigma1, 3), "\n\n")
cat("Estimated beta2:\n")
print(round(beta2, 3))
cat("\n")
cat("Estimated sigma2:\n", round(sigma2, 3), "\n")


# after your EM has converged, you have:
#   pi_k, beta1, beta2, sigma1, sigma2
# and your data X, y, n

library(numDeriv)

# 1. Pack parameters into a single vector
theta_hat <- c(
  log(pi_k / (1 - pi_k)), # reparametrize pi on the real line via logit
  beta1, # length p
  beta2, # length p
  log(sigma1), # log-scale
  log(sigma2)
)

# 2. Define the log‑likelihood of the *observed* data
loglik_obs <- function(theta) {
  # unpack
  eta <- theta[1]
  pi <- exp(eta) / (1 + exp(eta))
  p <- length(beta1)
  b1 <- theta[2:(1 + p)]
  b2 <- theta[(1 + p + 1):(1 + 2 * p)]
  s1 <- exp(theta[2 * p + 2])
  s2 <- exp(theta[2 * p + 3])
  # densitie
  mu1 <- X %*% b1
  mu2 <- X %*% b2
  d1 <- dnorm(y, mean = mu1, sd = s1, log = FALSE)
  d2 <- dnorm(y, mean = mu2, sd = s2, log = FALSE)
  ll_vec <- log(pi * d1 + (1 - pi) * d2)
  sum(ll_vec)
}

# 3. Compute the Hessian at the MLE
H <- hessian(loglik_obs, theta_hat)

# 4. Invert the *observed* information (–Hessian) to get vcov
vcov_theta <- solve(-H)

# 5. Back‐transform to original scales
#    variances of logit(pi) map to pi by delta method
#    likewise for log(sigma)
#    betas are on identity scale, so vcov entries correspond directly.

# Example: se(pi) ≈ sqrt( (dπ/dη)^2 * Var(η̂) ), with dπ/dη = π(1–π)
se_eta <- sqrt(vcov_theta[1, 1])
se_pi <- abs((pi_k * (1 - pi_k)) * se_eta)

# collect results
list(
  vcov_full = vcov_theta,
  se = c(
    pi = se_pi,
    beta1 = sqrt(diag(vcov_theta[2:(1 + p), 2:(1 + p)])),
    beta2 = sqrt(diag(vcov_theta[(2 + p):(1 + 2 * p), (2 + p):(1 + 2 * p)])),
    sigma1 = abs(sigma1) * sqrt(vcov_theta[2 + 2 * p, 2 + 2 * p]),
    sigma2 = abs(sigma2) * sqrt(vcov_theta[3 + 2 * p, 3 + 2 * p])
  )
)
