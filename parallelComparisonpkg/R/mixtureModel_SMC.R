# Use an SMC algorithm as outlined in [1] and [2] to infer the modes of a mixture of normals.
#
# Bibliogrphy:
# [1] M. A. Suchard, Q. Wang, C. Chan, J. Frelinger, A. Cron, and M. West,
# “Understanding GPU Programming for Statistical Computation: Studies in
# Massively Parallel Massive Mixtures,” J. Comput. Graph. Stat., vol. 19,
# no. February, pp. 419–438, 2010.
#
# [2] P. Del Moral, A. Doucet, and A. Jasra, “Sequential Monte Carlo samplers,”
# J. R. Stat. Soc. Ser. B (Statistical Methodol., vol. 68, no. 3,
# pp. 411–436, Jun. 2006.

#' A function for sampling from a mixture of normals.
#'
#' @param n Number of samples
#' @param Mu a vector of means
#' @param Sigma a vector of standard deviations
#' @param Weightts a vector of weights (should sum to 1)
#' @return A vector of length \code{n} containing n samples.
#'
#' @author Mathias C. Cronjager
sampleMM <- function(n = 1,
                     Mu = c(-3,0,3,6),
                     Sigma = rep(0.55,4),
                     Weights = rep(2^-2,4)){

  # Compute which distribution is sampled from in each of the n instances
  indices <- sample(1:length(Weights),n,replace = TRUE,prob = Weights)

  #Sample from the distributions and return the apropriate value.
  return( (rnorm(n) * Sigma[indices]) + Mu[indices] )
}

#referenceData <- sampleMM(100)

# MCMCSteps <- function(theta, steps){
#
# }

#' computes the likelihood-function p(y|Mu,Sigma,Weights) (up to normalization)
#' @param y a vector of observations
#' @param Mu a vector of means
#' @param Sigma a vector of standard deviations
#' @param Weightts a vector of weights (should sum to 1)
#' @param log Logical variable indicating if the likelihood or the
#'  log-likelihood should be computed. For reasons of numarical stability, the
#'  log-likelihood is preferable when y is very large.
#' @return a number corresponding to the (log) likelihood of observing y given
#'  Mu,Sigma,Weights
likelihood <- function(y, Mu, Sigma = 0.55, Weights = 2^-2, log = FALSE){

  if(max(abs(Mu)) > 10){
    if(log){
      return(-Inf)
    } else{
      return(0.0)
    }
  }
  #an auxiliary function for computing p(y[i] | Mu, Sigma, Weights)
  if(log){
    marginal_likelihood <- function(y2){
      log(sum(dnorm(x = y2, mean = Mu, sd = Sigma) * Weights))
    }
  } else {
    marginal_likelihood <- function(y2){
        sum(dnorm(x = y2, mean = Mu, sd = Sigma) * Weights)
    }
  }

  if(log){ #We return the log-likelihood.
    return(sum(sapply(y,FUN = marginal_likelihood)))
  } else { # We return the likelihood.
    return(prod(sapply(y,FUN = marginal_likelihood)))
  }
}

#' an auxiliary function for computing pi_n(x_(n-1)) / pi_(n-1)(x_(n-1))
#' \eqn{\propto} p( x_(n-1) | y ) ^ ((2n -1)/M^2). Used when updating weights.
#' @param y Observed data
#' @param x Our proposed value for Mu
pi_updateStep <- function(y,x,n,M = 200){
#   if(max(abs(x)) > 10){
#     return(0)
#   } else {
    return(exp( ( (2*n -1)/(M^2) ) * likelihood(y,Mu = x,log = TRUE)))
#   }
}

#' Returns TRUE if ESS(W) > c * N; N=length(W), c should be between 0 and 1.
#' The closer c is to 1, the more frequently we will resample.
ESS_is_largeEnough <- function(W,c=0.5){
  #ESS = 1/sum(W^2)
  return( 1 >= sum(W^2) * c * length(W))
}

#' Auxiliary function for performing 10 MH-steps with a symmetric Gaussian
#' proposal. The where the likelihood of the target-distribution is given by p( x_(n-1) | y ) ^ ((n/M)^2)
#' @param y Observed Data
#' @param x particle encoding the possible modes of the mixture to be approximated
#' @param n encodes the numerator of the exponent of the likelihood (the denominator is 200)
propagate_particle <- function(x, y, n, M = 200, steps = 10, sigma_MH = 1){

  #Metropolis-Hastings Steps
  for(i in 1:steps){

    #sample proposal
    xNew <- rnorm(n = length(x), mean = x, sd = sigma_MH)

    #compute acceptance-ratio
    alpha <- min(1,exp( (n/M)^2
                        * ( likelihood(y,Mu = xNew,log = TRUE)
                            - likelihood(y,Mu = x,log = TRUE))))

    #Choose wether to accept
    if(runif(1) < alpha){
      x <- xNew
    }

  }
  return(x)
}

## Auxiliary function for resampling
# resample <- function(X,W,N){
#   indices <- sample(N,size = N,prob = W,replace = TRUE)
#   return(X[,indices])
# }

#' Sample from the posterior p(Mu | y), using an SMC sampler as outlined in [1]
#' With simulated annealing.
#'
#' @param y A vector of observations
#' @param N the number of particles to be sampled
#' @param M the number of steps in the annealing process
#' @param c A constant between 0 and 1. We resample whenever the
#' effective sample size drops below c*N. Hence c close to 0 indicates that we will
#' never resample (effectively running N independent copies of the MH algorithm)
#' , whereas c close to 1 indicates that we will resample very often.
#' @param mh_steps how many MH steps are made when propagating particles forwards
#' @param mh_sigma the standard deviation of the proposal distribution to be used
#' in each MH step when propagating.
#'
#' @return an 4 x N matrix storing the samples columnwise
#'
#' @references
#' [1] P. Del Moral, A. Doucet, and A. Jasra, “Sequential Monte Carlo samplers,”
#'  J. R. Stat. Soc. Ser. B (Statistical Methodol., vol. 68, no. 3,
#'  pp. 411–436, Jun. 2006.
SMC_sampler <- function(y,N,M=200,c=0.5,mh_steps = 10, mh_sigma =1){

#   y <- as.array(y)
  
  #Initialize
  X <- matrix(data = runif(n = 4*N,min = -10,max = 10) )
#               nrow = 4,
#               ncol = N)
  dim(X) <- c(4,N)

  W <- rep(1/N,N)

  #Iterate from 1 to M
  n <- 1
  while(TRUE){

    # Resample if nessecary
    if(!ESS_is_largeEnough(W,c)){
      #We have to resample particles
      X <- X[,sample(N,size = N,prob = W,replace = TRUE)]
      W <- rep(1/N,N)
    }

    #Update particle Weights (can be done before propagation in this case)
    W <- W * apply(X,2,function(x) pi_updateStep(y,x,n,M))
    W <- W / sum(W)

    # Break out of loop when the M-th set of weights has been computed. We do
    # not propagate the last particles forward
    if(n==M) break
    
#    if(n==1) print(head(X))
    # Propagate Particles
    X <- apply(X,2,function(x) propagate_particle(x,y,n,M,steps = mh_steps,sigma_MH = mh_sigma))
#    if(n==1) print(head(X))
    
    #increment counter
    n <- n+1
  }
  return(list(X = X, W = W))
}

generatePlot <- function(N,N_y = 100){
  y <- sampleMM(N_y)
  smcData <- SMC_sampler(y,N)
  Ncolours <- 100
  blue <- 4/6
  yellow <- 1/6
  green <- 2/6
  red <- 1
  colVec <- rainbow(Ncolours,start = blue,end = green)[cut(smcData$W,breaks = Ncolours,labels=FALSE)]
  par(cex = 0.75 ,fg = "lightgrey", bg="black" ,col = "lightgrey", col.axis = "lightgrey", col.lab = "lightgrey",col.main = "lightgrey",col.sub = "lightgrey")
  plot(smcData$X[1,],smcData$X[2,],col = colVec, pch=".",ps=50)
  return(list(y,smcData,colVec))
}

SMC_sampler_FAST <- function(y,N){
  X <- vector(mode = "numeric",length = 4*N)
  W <- vector(mode = "numeric",length = N)

  outputRaw <- .C("smc_sampler_for_R",
                    yObs = as.double(y),
                    N_yObs = as.integer(length(y)),
                    N_particles = as.integer(N),
                    X_vec = as.double(X),
                    W = as.double(W),
                    t_SMC = as.double(0))

  X <- outputRaw$X_vec
  dim(X) <- c(N,4)
  W <- outputRaw$W
  runTime <- outputRaw$t_SMC
  
  return( list(X = X , W = W, runTime = runTime) )
}