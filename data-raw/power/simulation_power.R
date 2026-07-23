library(parallelly)
library(parallel)
library(future)
library(future.apply)
library(progressr)
library(dplyr)

library(rjson)

options(future.wait.interval=0L)

### Initialise progress bar
handlers(handler_progress(format="[:bar] :percent :eta :message"))


### Initialise computing cluster
# Local multicore backend. Each worker is a fresh R session that loads the
# installed catci package. Install it first with:  R CMD INSTALL .
library(catci)
# Cap at 5 (of 6 cores) to leave headroom for the OS / main process on an 8 GB machine.
plan(multisession, workers = 5)


### Define simulations

dir.create("results", showWarnings = FALSE)

reps <- 200

ns <- c(1000)
ds <- c(8)
xysettings <- rbind(c("lin", "lin"), # ordered cheapest-first (sin/sig settings are slower to fit)
                    c("lin", "vee"),
                    c("lin", "hat"),
                    c("vee", "vee"),
                    c("vee", "hat"),
                    c("hat", "hat"),
                    c("sig", "sig"),
                    c("sin", "sig"),
                    c("sin", "sin"))
intsettings <- c("binary_tree", "step")
strengths <- seq(0.2, 1.8, by = 0.2) # seq(0, 0.01, length.out = 11)

### Per-block timing log
timing_log <- file.path("results", "timing_power.log")
cat(sprintf("# power run started %s\n", format(Sys.time())), file = timing_log, append = TRUE)

for (xyi in seq_len(nrow(xysettings))){
  xsetting <- xysettings[xyi, 1]
  ysetting <- xysettings[xyi, 2]

  for (intsetting in intsettings){

    # Resume support: skip a block whose output CSV already exists.
    out_csv <- file.path("results", paste0("power_", xsetting, "_", ysetting, "_", intsetting, ".csv"))
    if (file.exists(out_csv)) {
      print(paste0("Skipping (already done): ", out_csv))
      next
    }

    print(paste0("Starting x = ", xsetting, ", y = ", ysetting, ", int = ", intsetting))
    block_t0 <- Sys.time()

    param_grid <- expand.grid(n = ns,
                              d = ds,
                              xsetting = xsetting,
                              ysetting = ysetting,
                              intsetting = intsetting,
                              strength = strengths)
    sim_df <- dplyr::slice(param_grid, rep(1:dplyr::n(), each = reps))
    sim_df$rep <- rep(1:reps, nrow(param_grid))


    with_progress( {
      prog_bar <- progressor(along=1:(nrow(sim_df)))
      sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, simplify=FALSE, FUN = function(x) {
        prog_bar()
        tryCatch({
        n <- as.numeric(x["n"])
        d <- as.numeric(x["d"])
        xsetting <- as.character(x["xsetting"])
        ysetting <- as.character(x["ysetting"])
        intsetting <- as.character(x["intsetting"])
        strength <- as.numeric(x["strength"])
        n_boot <- 100

        methods <- c("max", "euclid", "mGCM")
        if (intsetting == "binary_tree"){
          methods <- c("tree", methods)
        }
        if (intsetting == "step"){
          methods <- c("ordinal", methods)
        }

        data <- simulate_data(n = n,
                              xnum_class = d,
                              ynum_class = d,
                              xsetting = xsetting,
                              ysetting = ysetting,
                              strength = strength,
                              intsetting = intsetting,
                              permute = FALSE)

        xparams <- rjson::fromJSON(file = paste0("data-raw/tuning/n", n, "_numclass", d, "/tune_", xsetting, "_results.json"))$xgb
        yparams <- rjson::fromJSON(file = paste0("data-raw/tuning/n", n, "_numclass", d, "/tune_", ysetting, "_results.json"))$xgb

        stats <- catci:::formulate_statistics(data = data,
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
        }, error = function(e) list(error = conditionMessage(e)))
      })})

    sim_res_df <- cbind(sim_df, data.table::rbindlist(sim_res, fill = TRUE))
    write.csv(sim_res_df, out_csv, row.names=FALSE)

    block_min <- as.numeric(difftime(Sys.time(), block_t0, units = "mins"))
    n_err <- if ("error" %in% names(sim_res_df)) sum(!is.na(sim_res_df$error)) else 0L
    log_line <- sprintf("%s  x=%-3s y=%-3s int=%-11s  rows=%-5d errors=%-4d  %6.1f min",
                        format(Sys.time()), xsetting, ysetting, intsetting, nrow(sim_df), n_err, block_min)
    cat(log_line, "\n", file = timing_log, append = TRUE, sep = "")
    print(log_line)

  }
}

plan(sequential)  # shut down workers
