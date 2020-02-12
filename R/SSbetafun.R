beta.fun <- function(X, a, b, c, d){
  ifelse(X > b & X < d, 
         (((X - b)/(c - b) * (d - X)/(d - c))^((d - c)/(c - b)))^a,
         0)
}
