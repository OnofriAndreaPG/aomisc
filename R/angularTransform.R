#' Angular transformation for percentages
#' 
#' @param percentage A number.
#' @returns A number.
#' @examples
#' angularTransform(25)
#' angularTransform(35.2)

angularTransform <- function(percentage){
  # angular transformation for proportions
  if(any(percentage < 0) | any(percentage > 100)) stop("This transformation can only be used with percentages from 0 to 100%")
  result <- asin(sqrt(percentage/100))*(180/pi)
  return(result)
}
invAngularTransform <- function(angle){
  # angular transformation for proportions
  # if(angle < 0 | percentage > 100) stop("This transformation can only be used with percentages from 0 to 100%")
  percentage <- sin(angle/180*pi)^2*100
  return(percentage)
}

angular <- function(){
  # angular transformation for proportions
  list(linkfun = function(mu) angularTransform(mu),
       linkinv = function(eta) invAngularTransform(eta),
       mu.eta = function(eta) 2 * (cos(eta/180 * pi) * (1/180 * pi) * sin(eta/180 * pi)) * 100 
        )
}

