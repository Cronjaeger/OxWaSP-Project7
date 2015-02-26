## a series of tests for benchmarking

benchmark_SMC <- function(N_seq = 2^(4:15),y = sampleMM(100),returnEverything = FALSE,verbose = TRUE,checkpointing=FALSE){
  require(parallelComparisonpkg)
  outList <- list()
  for(N in N_seq){
    out <- SMC_sampler_FAST(y,N)
    outList <- c(outList,list(out))
    outStr <- paste("N =",N,"t =",out$runTime)
    print(outStr)
    if(checkpointing) save(outList,file="benchmarkCheckpoint.RData")
  }
  if(returnEverything) return(outList) #return the whole shebang
  else return(sapply(outList,function(x) x$runTime)) #return only the runtimes
}

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