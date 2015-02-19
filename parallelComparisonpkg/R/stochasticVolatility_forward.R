## Simulate a stochastic volatility modell forwards in time

## Default Values
steps_D <- 200
x0_D <- rep(0.0,3)
B_D <- matrix(data = c(1,0.5,0.5,0.2,0.8,0,1,0.5,0.6,0.7,0,0,1,0.3,0.5),nrow = 5,ncol = 3)
U_D <- matrix(data = c(0.5,0.2,0.1,0.2,0.5,0.2,0.1,0.2,0.5), nrow = 3, ncol = 3)
psi_D <- rep(0.5,5)
phi_D <- rep(0.9,3)

#' Simulate a stochastic volatitlity modell forwards
#' 
#' @author Mathias C. Cronj\"ager
#' @description Simulates A stochastic volatilityu model as outlined in section 4.2 of (Lee et al. 2010).
#' 
#' x[t] ~ Normal(diag(phi)*x[t-1], U)
#' 
#' y[t] ~ Normal(0, B * diag(exp(x[t])) * B' + diag(psi)).
#' 
#' All arguments default to the arguments used in (Lee et al. 2010).
#' 
#' @param steps Number of timesteps to be taken
#' @param x0 initial conditions
#' @param B load matrix.
#' @param U covariance matrix of x-update
#' @param psi weights added to covariance-matrix of observations
#' @param phi x-rescaling before updating 
#' 
#' @return Returns a list of two matrices "X" and "Y". The column-index of each matrix corresponds to time.
#' 
#' @export simulateForwards
simulateForwards <- function(steps = steps_D,
                             x0 = x0_D,
                             B = B_D,
                             U = U_D,
                             psi = psi_D,
                             phi = phi_D){
 
  require(MASS) # mvrnorm is in this package.

  M <- dim(B)[1]
  K <- dim(B)[2]
  Psi <- diag(psi)
  Phi <- diag(phi)

  # TODO: check if dimensions mach (or give warnings if they dont?)
  
  #initialize objects
  X <- matrix(data = c(x0,rep(NA,(steps)*K)), nrow = K, ncol = steps+1)
  Y <- matrix(data = NA, nrow = M, ncol = steps + 1)
  H <- diag(NA,nrow = K,ncol = K)
  zero_M <- rep(0.0,M)
  
  #loop for computing values of x[2] , ... , x[steps+1]
  for(t in 2:(steps + 1)){
    X[,t] <- mvrnorm(mu = Phi %*% X[,(t-1)], Sigma = U)
  }
  
  #loop for computing values of y[1] , ... , y[steps+1]
  for(t in 1:(steps +1 )){
    H <- diag(exp(X[,t]))
    Y[,t] <- mvrnorm(mu = zero_M, Sigma = B %*% H %*% t(B) + Psi)
  }
  
  return(list("X" =X,"Y" = Y))
}