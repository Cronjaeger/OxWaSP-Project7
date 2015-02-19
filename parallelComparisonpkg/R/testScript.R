## This is a test script for checking if we can figure out how to integrate .c code into an R-package


#' A wrapper function for calling square.c and parsing output.
squareVector <- function(X){
  # Squares a vector...
  n = length(X)
  returnVector = vector(mode = "numeric",length = n)
  returnVector = .C("square",N = as.integer(n), x = as.double(X))$x
  return(returnVector)
}

## we test if things work
testSq <- function(){
  vec = squareVector(1:10)
  print(vec)
}