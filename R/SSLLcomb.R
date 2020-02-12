LL.Comb.fun <- function(predictor, b1, e1, b2, e2) {
  x <- predictor
  Y1 <- 1/(1+exp(-b1*(log(x+0.000001)-log(e1))))
  Y2 <- 1/(1+exp(-b2*(log(x+0.000001)-log(e2))))
  Y1 + Y2
}
