# One-off: re-run the lin_lin_binary_tree power setting with the FIXED tree
# method (CODE_REVIEW.md finding #1) to see whether Figure 2 changes.
#
# - Uses devtools::load_all() so the FIXED R/ source is what runs, NOT the old
#   installed catci package.
# - Adds "ordinal" alongside "tree" so tree-vs-ordinal is measured on identical
#   data within one run (the power sim normally runs only "tree" for binary_tree).
# - Writes to a NEW filename; never touches results/power_lin_lin_binary_tree.csv
#   (the existing buggy-tree baseline).
# Run from the project root.

library(future)
library(future.apply)
library(progressr)
library(rjson)

options(future.wait.interval = 0L)
options(parallelly.fork.enable = TRUE)   # allow forking from Rscript
handlers(handler_progress(format = "[:bar] :percent :eta :message"))

devtools::load_all(quiet = TRUE)   # FIXED source
# multicore (forked workers) so each worker inherits the load_all'd namespace
# without needing the package installed. fit_xgboost uses nthread = 1, so
# forking is safe here.
plan(multicore, workers = 5)

out_csv <- file.path("results", "power_lin_lin_binary_tree_FIXED.csv")
stopifnot(!file.exists(out_csv))   # refuse to clobber

reps <- 200
n <- 1000
d <- 8
xsetting <- "lin"; ysetting <- "lin"; intsetting <- "binary_tree"
strengths <- seq(0.2, 1.8, by = 0.2)

param_grid <- expand.grid(n = n, d = d, xsetting = xsetting, ysetting = ysetting,
                          intsetting = intsetting, strength = strengths)
sim_df <- dplyr::slice(param_grid, rep(1:dplyr::n(), each = reps))
sim_df$rep <- rep(1:reps, nrow(param_grid))

t0 <- Sys.time()
with_progress({
  prog_bar <- progressor(along = 1:nrow(sim_df))
  sim_res <- future_apply(sim_df, MARGIN = 1, future.seed = TRUE, simplify = FALSE,
    FUN = function(x) {
      prog_bar()
      tryCatch({
        n <- as.numeric(x["n"]); d <- as.numeric(x["d"])
        xsetting <- as.character(x["xsetting"]); ysetting <- as.character(x["ysetting"])
        intsetting <- as.character(x["intsetting"]); strength <- as.numeric(x["strength"])
        n_boot <- 100

        # tree AND ordinal on the same data, plus the non-adaptive comparators.
        methods <- c("tree", "ordinal", "max", "euclid", "mGCM")

        data <- simulate_data(n = n, xnum_class = d, ynum_class = d,
                              xsetting = xsetting, ysetting = ysetting,
                              strength = strength, intsetting = intsetting,
                              permute = FALSE)

        xparams <- rjson::fromJSON(file = paste0("data-raw/tuning/n", n, "_numclass", d, "/tune_", xsetting, "_results.json"))$xgb
        yparams <- rjson::fromJSON(file = paste0("data-raw/tuning/n", n, "_numclass", d, "/tune_", ysetting, "_results.json"))$xgb

        stats <- formulate_statistics(data = data, xnum_class = d, ynum_class = d,
                                      method = "xgb", xparams = xparams, yparams = yparams,
                                      nfolds = 5, normalise = FALSE)

        evaluate_sim(data = stats, dx = d, dy = d, n_boot = n_boot, methods = methods)
      }, error = function(e) list(error = conditionMessage(e)))
    })
})

sim_res_df <- cbind(sim_df, data.table::rbindlist(sim_res, fill = TRUE))
write.csv(sim_res_df, out_csv, row.names = FALSE)
cat(sprintf("DONE %s  rows=%d  %.1f min -> %s\n",
            format(Sys.time()), nrow(sim_df),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), out_csv))
plan(sequential)
