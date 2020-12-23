angularTransform <- function(percentage){
  # angular transformation for proportions
  result <- asin(sqrt(percentage/100))*(180/pi)
  return(result)
}

ma <- function(x, n = 5, sides = 2){
  # moving average for a vector
  res <- stats::filter(x, rep(1 / n, n), method = "convolution", sides = sides)
  as.numeric(res)
  }


