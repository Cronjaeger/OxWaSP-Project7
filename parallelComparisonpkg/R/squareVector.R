## This is a test script for checking if we can figure out how to integrate .c code into an R-package

#' A wrapper function for calling square.c and parsing output.
squareVector <- function(X){
  # Squares a vector...
  n = length(X)
  returnVector = vector(mode = "numeric",length = n)
  returnVector = .C("square",N = as.integer(n), x = as.double(X))$x
  return(returnVector)
}

#' Test if square.c has been compiled and is callable.
testSq <- function(n=10){
  vec = squareVector(1:n)
  print(tail(vec))
}