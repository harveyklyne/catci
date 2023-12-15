library(parallelly)
library(parallel)
library(future)
library(future.apply)
library(progressr)
library(tidyverse)

library(rjson)
library(reshape2)

options(future.wait.interval=0L)

### Initialise progress bar
handlers(handler_progress(format="[:bar] :percent :eta :message"))


### Initialise computing cluster
worker_names <- NULL # put your worker names here
working_directory <- NULL # put your working directory here


# # Wake workers
for (worker in worker_names){
  system(paste0("/alt/bin/wake ", worker))
}
Sys.sleep(60)

nodes_per_worker <- 1

my_cluster <- parallelly::makeClusterPSOCK(
  rep(worker_names, nodes_per_worker),
  outfile="", # this option ensures that error messages and status messages are printed to the console of the host (i.e. the computer the running the script)
  homogeneous=FALSE) # homogeneous = FALSE is crucial if the operating system of the host (i.e. the computer running the script) differs from the operating system of the workers


devtools::load_all()

clusterEvalQ(my_cluster, library(devtools))
clusterExport(my_cluster, c("working_directory"), envir = .GlobalEnv)
clusterEvalQ(my_cluster, setwd(working_directory))
clusterEvalQ(my_cluster, devtools::load_all())

plan(cluster, workers=my_cluster)

#####

settings <- c("sin", "sig", "lin", "vee", "hat")

n_tr <- 800
n_te <- 5000
num_class <- 8
reps <- 100

for (setting in settings){

  xgb <- tune_xgb(target = "x",
                  n_tr = n_tr,
                  n_te = n_te,
                  num_class = num_class,
                  setting = setting,
                  reps = reps)

  output <- list("xgb" = xgb)

  jsonData <- toJSON(output)

  write(jsonData, paste0("tune_", setting, "_results.json"))

}



stopCluster(my_cluster)



