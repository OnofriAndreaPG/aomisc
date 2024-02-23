# This function extends the anova method to the objects
# fitted with 'aov()'
anova.aovlist <- function(object, ...){
  summary(object)
}