# This function extends the anova method to the objects
# fitted with 'aov()'
anova.aovList <- function(object){
  summary(object)
}