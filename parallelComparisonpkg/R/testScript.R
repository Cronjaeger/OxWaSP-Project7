## a series of tests for benchmarking


#' A function for calling SMC_sampler_FAST with a renge of inputs and keeping track of runtimes
#' 
#' @param N_seq a sequence of integers. Each element is passed to SMC_sampler_FAST().
#' @param y the samoke of observations to be passed to SMC_sampler_FAST() defaults to a new sample of 100 observations.
#' @param returnEverything a logical variable indicating if we want all output produced
#' (particles, weights,...) returnsed. If true the method returns a list of lists. Each list
#' contains the output of 1 run. If false, only a vector of runtimes is returned.
#' @param verbose wether to print output when generated. The ith line printed will hvae the form
#' @param checkpointing If set to true, the output will be saved to "benchmarkCheckpoint.RData"
#' after each run.
#' 
#' "N = N_seq[i], t = SMC_sampler_FAST(y,N[i])$runTime"
benchmark_SMC <- function(N_seq = 2^(1:15),y = sampleMM(100),returnEverything = FALSE,verbose = TRUE,checkpointing=FALSE,useR=FALSE){
  require(parallelComparisonpkg)
  outList <- list()
  if(useR){
    sampler <- function(y,N) SMC_sampler(y,N)
  } else {
    sampler <- function(y,N) SMC_sampler_FAST(y,N)
  }
  for(N in N_seq){
    out <- sampler(y,N)
    outList <- c(outList,list(out))
    outStr <- paste("N =",N,"t =",out$runTime)
    if(verbose) print(outStr)
    if(checkpointing) save(outList,file="benchmarkCheckpoint.RData")
  }
  if(returnEverything) return(outList) #return the whole shebang
  else return(sapply(outList,function(x) x$runTime)) #return only the runtimes
}


#' A function for generating a log-log-plot of number of particles Vs. runtime
#' @param N_seq the sequence of run-times to be investigated.
generateRuntimePlot <- function(N_seq = 2^(4:13)){
  outTimes <- benchmark_SMC(N_seq,verbose = FALSE)
  plot(x = N_seq, 
       y = outTimes,
       type = "p",
       #bty="n",
       pch = 3,
       log="xy",
       xlab = "",
       ylab = "",
       main = "runtime (sec) Vs #particles (log-log)")
}