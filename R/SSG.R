# Gompertz equation
G.4.fun <- function(X, b, c, d, e) {
  c + (d - c) * (exp (- exp ( - b * ( X - e))))
}

G.3.fun <- function(X, b, d, e) {
  d * (exp (- exp ( - b * ( X - e))))
}

G.2.fun <- function(X, b, e) {
  exp (- exp ( - b * ( X - e)))
}
