ma <- function(x, n = 5, sides = 2){
  # moving average for a vector with the convolution method
  res <- stats::filter(x, rep(1 / n, n), 
                       method = "convolution", sides = sides)
  as.numeric(res)
  }


