# Simples função para tentar avaliar uma expressão e retornar NULL em caso de erro
# Uso: `tryNull(expr)` é equivalente a `tryCatch(expr, error = function(e) NULL)`
# `(1 + "1") |> tryNull()` retorna NULL em vez de lançar um erro
tryNull <- \(expr) tryCatch(expr, error = \(e) NULL)
