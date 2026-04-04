options(
  testthat.default_reporter = "progress",
  encoding = "UTF-8"
)

# Evita paralelismo dentro dos testes
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

devtools::load_all()
devtools::test()
