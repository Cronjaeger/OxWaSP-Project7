## ----sampleY-------------------------------------------------------------
sampleMM(10)

y <- sampleMM(2^15)
hist(y,breaks = 100,freq = F,
     col="grey",border="grey",bty="n")
rug(y[1:100])

## ----FastVsSlow,eval=FALSE-----------------------------------------------
#  y <- sampleMM(100)
#  N_particles <- 1000
#  out_slow <- SMC_sampler(y,N_particles)
#  out_FAST <- SMC_sampler_FAST(y,N_particles)
#  print(paste("t_slow =",out_slow$runTime,"sec"))
#  print(paste("t_fast =",out_FAST$runTime,"sec"))

## ----scatterplot,cache=TRUE,warning=FALSE--------------------------------
y <- sampleMM(100)
N_particles <- 1000
out_FAST <- SMC_sampler_FAST(y,N_particles)
generatePlot(SMC_Data = out_FAST)

## ----benchmarking,cache=TRUE---------------------------------------------
#prints time in secconds for each N in N_seq. stores runtimes as a vector
times <- benchmark_SMC(N_seq = 20*1:5)
times

#runs the above and prints a log-log-plot of runtime vs number of particles
generateRuntimePlot(N_seq = 20*1:5)

