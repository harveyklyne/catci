tune_xgb <- function(target, n_tr, n_te, num_class, setting, reps){

  print(paste0("Tuning ", target, " XGB ", setting))

  sim_df <- expand.grid(reps=1:reps)

  eta <- 0.01
  depths <- 1:4
  gammas <- seq(0,3, by=0.5)
  maxrounds <- 1000

  if (setting %in% c("lin", "vee", "hat")){
    intsetting <- "step"
  }
  if (setting %in% c("sin", "sig")){
    intsetting <- "binary_tree"
  }

  with_progress( {
    prog_bar <- progressor(along=1:(nrow(sim_df)))
    sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
      prog_bar()

      out <- suppressWarnings(xgb_tune_grid(target = target,
                                            n_tr = n_tr,
                                            n_te = n_te,
                                            xnum_class = num_class,
                                            ynum_class = num_class,
                                            xsetting = setting,
                                            ysetting = setting,
                                            strength = 0.5,
                                            intsetting = intsetting,
                                            eta = eta,
                                            depths = depths,
                                            gammas = gammas,
                                            maxrounds = maxrounds))
      return(out)
    })
  })

  mean_mse <- Reduce("+", sim_res) / length(sim_res)
  mean_mse_out <- melt(mean_mse)
  mean_mse_out[,2] <- gammas[mean_mse_out[,2]]
  names(mean_mse_out) = c("depth", "gamma", "nround", "mse")

  # write.csv(mean_mse_out, paste0("xgb_tune_", setting, "_results.csv"), row.names=FALSE)

  tmp <- arrayInd(which.min(mean_mse), dim(mean_mse))
  depth <- depths[tmp[1]]
  gamma <- gammas[tmp[2]]
  nrounds <- seq(maxrounds)[tmp[3]]

  result <- list("eta" = eta,
                 "max.depth" = depth,
                 "gamma" = gamma,
                 "nrounds" = nrounds)

  return(result)

}

xgb_tune_grid <- function(target,
                          n_tr,
                          n_te,
                          xnum_class,
                          ynum_class,
                          xsetting,
                          ysetting,
                          strength,
                          intsetting,
                          eta,
                          depths,
                          gammas,
                          maxrounds){

  data <- simulate_train_test(n_tr = n_tr,
                              n_te = n_te,
                              xnum_class = xnum_class,
                              ynum_class = ynum_class,
                              xsetting = xsetting,
                              ysetting = ysetting,
                              strength = strength,
                              intsetting = intsetting)

  if (target == "x"){
    test_error <- xgb_tune_grid_evaluate(train = xgboost::xgb.DMatrix(data = data$train$z,
                                                                      label = data$train$x - 1),
                                         test = xgboost::xgb.DMatrix(data = data$test$z,
                                                                     label = data$test$x - 1),
                                         num_class = xnum_class,
                                         eta = eta,
                                         depths = depths,
                                         gammas = gammas,
                                         maxrounds = maxrounds)
  }
  if (target == "y"){
    test_error <- xgb_tune_grid_evaluate(train = xgboost::xgb.DMatrix(data = data$train$z,
                                                                      label = data$train$y - 1),
                                         test = xgboost::xgb.DMatrix(data = data$test$z,
                                                                     label = data$test$y - 1),
                                         num_class = ynum_class,
                                         eta = eta,
                                         depths = depths,
                                         gammas = gammas,
                                         maxrounds = maxrounds)
  }

  return(test_error)
}

xgb_tune_grid_evaluate <- function(train, test, num_class, eta, depths, gammas, maxrounds){

  test_error <- array(NA, c(length(depths), length(gammas), maxrounds))

  for (i in seq(length(depths))){
    for (j in seq(length(gammas))){
      fit <- xgboost::xgb.train(data = train,
                                watchlist = list(train = train, test = test),
                                eta = eta,
                                max.depth = depths[i],
                                gamma = gammas[j],
                                nrounds = maxrounds,
                                objective = "multi:softprob",
                                eval_metric = "mlogloss",
                                num_class = num_class,
                                nthread = 2,
                                verbose = 0)
      test_error[i,j,] <- fit$evaluation_log$test_mlogloss
    }
  }

  return(test_error)
}
