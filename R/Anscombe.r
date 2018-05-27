AnscombeTukey<-function(formula, data){

## model is a lm object
model <- lm(formula, data=data)
ENNE <- length(model$residuals)
NI <- model$df.residual
premium <- 0.025
prob <- (premium * NI / ENNE)
z <- abs(qnorm(prob))
K <- 1.4 + 0.85 * z
c <- K * (1 - (K^2 - 2) / (4 * NI)) * sqrt(NI / ENNE)
devst <- sqrt((t(model$residuals)%*%model$residuals)/model$df.residual)
maxteor = c * devst

## Verifica se esistono residui superiori al massimo ammesso
ResMax <- max(model$residuals)
Pos <- which.max(model$residuals)
if (abs(ResMax) < maxteor ) {
  cat(paste("No outliers were found", "\n")) }
else {
  cat(paste("The data in position n. ", Pos, "may be considered an outlier", "\n"))
  new <- data[Pos,]
  data <- rbind(data[1:Pos-1,], data[(Pos+1):length(data[,1]), ])
  model2 <- update(model)  
  Correct <- predict(model2, new)
  cat(paste("The correct data is ", Correct, "\n"))  
} 
## diff<-ResMax-maxteor
#Find correct data



}