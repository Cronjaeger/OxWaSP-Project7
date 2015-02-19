# ## This is a test script
# 
# 
# 
# squareVector <- function(X){
#   # Squares a vector...
#   #library(parallelComparisonpkg)
#   n = length(X)
#   returnVector = vector(mode = "numeric",length = n)
#   returnVector = .C("square",N = as.integer(n), x = as.double(X))$x
#   return(returnVector)
# }
# 
# ## we test
# 
# vec = 1:10
# print(squareVector(vec))