# basic source script function
mySourceScript <- function(x)

{
  source("/home/jason/R/runoptGPP/testing/Null_extent_03_gpp_pcm_opt.R")
}

#mySourceScript()

Sys.time()
# ClusterMQ-Call
clustermq::Q(fun = mySourceScript,
             x=1,
             n_jobs = 1,
             template = list(n_cpus = 32,
                             memory = 120000,
                             partition = "all",
                             log_file = "/home/jason/R/sedconnect/clustermq_debug.log", # adapt here
                             job_name = "GrdOpt")
)

Sys.time()

fx = function(x) {
  library(doParallel)
  registerDoParallel(4)
  foo <- foreach(1:4) %dopar% {
    foo = runif(1000000000)
  }
  stopImplicitCluster()

}

clustermq::Q(
  fun = fx, x = 1, n_jobs = 1,
  template = list(n_cpus = 4, partition = "all", memory = 10000))
