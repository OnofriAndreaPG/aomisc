# Bragg's equation
bragg.3.fun <- function(X, b, d, e){
  d * exp(- b * (X - e)^2)
}

bragg.4.fun <- function(X, b, c, d, e){
  c + (d - c) * exp(- b * (X - e)^2)
}
