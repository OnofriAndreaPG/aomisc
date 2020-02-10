E.4.fun <- function(X, b, c, d, e) {
  c + (d - c) * (1 - exp( - exp ( b * ( X - e))))
}

E.3.fun <- function(X, b, d, e) {
  d * (1 - exp( - exp ( b * ( X - e))))
}

E.2.fun <- function(X, b, e) {
  1 - exp( - exp ( b * ( X - e)))
}
