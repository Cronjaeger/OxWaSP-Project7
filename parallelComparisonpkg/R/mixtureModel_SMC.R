# Use an SMC algorithm as outlined in [1] and [2] to infer the modes of a mixture of normals.
# 
# Bibliogrphy:
# [1] M. A. Suchard, Q. Wang, C. Chan, J. Frelinger, A. Cron, and M. West, “Understanding GPU Programming for Statistical Computation: Studies in Massively Parallel Massive Mixtures,” J. Comput. Graph. Stat., vol. 19, no. February, pp. 419–438, 2010.
# [2] P. Del Moral, A. Doucet, and A. Jasra, “Sequential Monte Carlo samplers,” J. R. Stat. Soc. Ser. B (Statistical Methodol., vol. 68, no. 3, pp. 411–436, Jun. 2006.

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