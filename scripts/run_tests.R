options(testthat.default_reporter = "progress")
options(encoding = "UTF-8")

# Evita paralelismo dentro dos testes
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

testthat::test_dir("tests/testthat", reporter = "progress")
