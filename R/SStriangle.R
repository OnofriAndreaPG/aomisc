triangle.fun <- function(X, k, Xb, Xo, Xc){
  Y <- ifelse(X < Xb, 0,  
              ifelse(X < Xo, k * (X - Xb)/(Xo - Xb),
                     ifelse (X < Xc, k * (Xc - X)/(Xc - Xo), 0))) 
  Y
  }
