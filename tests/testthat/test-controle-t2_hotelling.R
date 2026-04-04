testthat::test_that("t2_hotelling_boot retorna data.table com colunas esperadas", {
  phi_boot <- c(0.1, 0.2, 0.3)
  theta_boot <- c(0.4, 0.5, 0.6)
  vcov_inv <- diag(2)

  out <- t2_hotelling_boot(
    phi_boot = phi_boot,
    theta_boot = theta_boot,
    phi = 0.2,
    theta = 0.5,
    vcov_inv = vcov_inv
  )

  testthat::expect_s3_class(out, "data.table")
  testthat::expect_equal(nrow(out), 1)
  testthat::expect_true(all(c(
    "numero_de_bootstrap_validos",
    "fora_de_controle",
    "t2_limite",
    "t2_controle"
  ) %in% names(out)))
})

testthat::test_that("t2_hotelling_boot calcula corretamente em caso simples", {
  # Sejam:
  #
  # phi_boot = c(0, 1)
  # theta_boot = c(0, 1)
  # média mu_boot = (0.5, 0.5)
  # vcov_inv = I
  #
  # Então, seja d a diferença entre o controle e a média, d = (phi - 0.5, theta - 0.5)
  #
  # T²_1: d = (-0.5, -0.5) => d' I d => T² = 0.5
  # T²_2: d = (0.5, 0.5) => d' I d => T² = 0.5
  #
  # Logo:
  #
  # t2_boot = c(0.5, 0.5)
  # t2_limite = 0.5
  #
  # Se phi = 0.5, theta = 0.5, então:
  #
  # diff_controle = (0, 0)
  # t2_controle = 0
  # fora_de_controle = FALSE

  phi_boot <- c(0, 1)
  theta_boot <- c(0, 1)
  vcov_inv <- diag(2)

  out <- t2_hotelling_boot(
    phi_boot = phi_boot,
    theta_boot = theta_boot,
    phi = 0.5,
    theta = 0.5,
    vcov_inv = vcov_inv
  )

  testthat::expect_equal(out$numero_de_bootstrap_validos, 2)
  testthat::expect_equal(out$t2_limite, 0.5, tolerance = 1e-12)
  testthat::expect_equal(out$t2_controle, 0, tolerance = 1e-12)
  testthat::expect_false(out$fora_de_controle)
})

testthat::test_that("t2_hotelling_boot marca fora_de_controle quando t2_controle excede limite", {
  # Sejam phi = 2, theta = 2 e phi_boot = theta_boot = c(0, 1), então:
  #
  # média mu_boot = (0.5, 0.5)
  # diferença = (1.5, 1.5)
  # T² = 1.5^2 + 1.5^2 = 4.5
  #
  # Como 4.5 > 0.5, deve indicar fora de controle.

  phi_boot <- c(0, 1)
  theta_boot <- c(0, 1)
  vcov_inv <- diag(2)

  out <- t2_hotelling_boot(
    phi_boot = phi_boot,
    theta_boot = theta_boot,
    phi = 2,
    theta = 2,
    vcov_inv = vcov_inv
  )

  testthat::expect_equal(out$t2_limite, 0.5, tolerance = 1e-12)
  testthat::expect_equal(out$t2_controle, 4.5, tolerance = 1e-12)
  testthat::expect_true(out$fora_de_controle)
})

testthat::test_that("t2_hotelling_boot nao marca fora_de_controle quando t2_controle e igual ao limite", {
  # A regra é:
  #
  # Se t2_controle > t2_limite então fora_de_controle = TRUE

  phi_boot <- c(0, 1)
  theta_boot <- c(0, 1)
  vcov_inv <- diag(2)

  out <- t2_hotelling_boot(
    phi_boot = phi_boot,
    theta_boot = theta_boot,
    phi = 1,
    theta = 1,
    vcov_inv = vcov_inv
  )

  testthat::expect_equal(out$t2_limite, 0.5, tolerance = 1e-12)
  testthat::expect_equal(out$t2_controle, 0.5, tolerance = 1e-12)
  testthat::expect_false(out$fora_de_controle)
})

testthat::test_that("t2_hotelling_boot usa corretamente a matriz vcov_inv", {
  # matriz inversa não identidade
  # Se vcov_inv = diag(c(2, 3)) e phi_boot = theta_boot = c(0, 2)
  #
  # Teremos média mu_boot = (1, 1)
  #
  # E que, para cada par (phi_boot[i], theta_boot[i]):
  #
  # diff = (-1, -1) ou (1, 1)
  # T² = 2*(1^2) + 3*(1^2) = 5
  #
  # Então o limite será 5.
  #
  # Se o controle for (2, 1):
  #
  # diff = (1, 0) => T² = 2*(1^2) + 3*(0^2) = 2

  phi_boot <- c(0, 2)
  theta_boot <- c(0, 2)
  vcov_inv <- diag(c(2, 3))

  out <- t2_hotelling_boot(
    phi_boot = phi_boot,
    theta_boot = theta_boot,
    phi = 2,
    theta = 1,
    vcov_inv = vcov_inv
  )

  testthat::expect_equal(out$t2_limite, 5, tolerance = 1e-12)
  testthat::expect_equal(out$t2_controle, 2, tolerance = 1e-12)
  testthat::expect_false(out$fora_de_controle)
})
