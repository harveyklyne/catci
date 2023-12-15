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


### Define simulations

devtools::load_all()

reps <- 1000

ns <- c(1000)
ds <- c(8)
xysettings <- rbind(#c("sin", "sin"),
                    #c("sin", "sig"),
                    #c("sig", "sig"),
                    c("lin", "lin"),
                    c("lin", "vee"),
                    c("lin", "hat"),
                    c("vee", "vee"),
                    c("vee", "hat"),
                    c("hat", "hat"))

for (xyi in seq_len(nrow(xysettings))){
  xsetting <- xysettings[xyi, 1]
  ysetting <- xysettings[xyi, 2]

  print(paste0("Starting x = ", xsetting, ", y = ", ysetting))

  param_grid <- expand.grid(n = ns,
                            d = ds,
                            xsetting = xsetting,
                            ysetting = ysetting)
  sim_df <- dplyr::slice(param_grid, rep(1:dplyr::n(), each = reps))
  sim_df$rep <- rep(1:reps, nrow(param_grid))

  with_progress( {
    prog_bar <- progressor(along=1:(nrow(sim_df)))
    sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, simplify=FALSE, FUN = function(x) {
      prog_bar()
      n <- as.numeric(x["n"])
      d <- as.numeric(x["d"])
      xsetting <- as.character(x["xsetting"])
      ysetting <- as.character(x["ysetting"])
      n_boot <- 100

      methods <- "tree" # c("tree", "ordinal", "max", "euclid", "mGCM")


      data <- simulate_data(n = n,
                            xnum_class = d,
                            ynum_class = d,
                            xsetting = xsetting,
                            ysetting = ysetting,
                            strength = 0,
                            intsetting = "step",
                            permute = FALSE)

      xparams <- rjson::fromJSON(file = paste0("data-raw/tuning/n", n, "_numclass", d, "/tune_", xsetting, "_results.json"))$xgb
      yparams <- rjson::fromJSON(file = paste0("data-raw/tuning/n", n, "_numclass", d, "/tune_", ysetting, "_results.json"))$xgb

      stats <- formulate_statistics(data = data,
                                    xnum_class = d,
                                    ynum_class = d,
                                    method = "xgb",
                                    xparams = xparams,
                                    yparams = yparams,
                                    nfolds = 5,
                                    normalise = FALSE)

      values <- evaluate_sim(data = stats,
                             dx = d,
                             dy = d,
                             n_boot = n_boot,
                             methods = methods)

      return(values)
    })})

  sim_res_df <- cbind(sim_df, data.table::rbindlist(sim_res, fill = TRUE))
  write.csv(sim_res_df, paste0("size_", xsetting, "_", ysetting, ".csv"), row.names=FALSE)

}

stopCluster(my_cluster)
