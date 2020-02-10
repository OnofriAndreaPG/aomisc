# Brain-Cousens model for hormesis
BC.5.fun <- function(predictor, b, c, d, e, f) {
  x <- predictor
  c + (d - c + f*x)/(1 + exp(-b * (log(x + 0.00001) - log(e))))
}

BC.4.fun <- function(predictor, b, d, e, f) {
  x <- predictor
  (d + f*x)/(1 + exp(-b * (log(x + 0.00001) - log(e))))
}

BC.3.fun <- function(predictor, b, e, f) {
  x <- predictor
  (f * x)/(1 + exp(-b * (log(x + 0.00001) - log(e))))
}
