fit_xgboost <- function(X, y, num_class = NULL, params) {
  #' Fit pre-tuned XGBoost regression for use in simulations.
  #'
  #' @param X matrix of covariates.
  #' @param y integer vector of responses in {1, ..., num_class}.
  #' @param num_class integer number of levels of desired output.
  #' @param params XGBoost hyperparameters.
  #'
  #' @return list containing a function "fit" which takes matrix input of the
  #'     same width as X, and returns a matrix of probability estimates, with each
  #'     column corresponding to a label of y.

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package \"xgboost\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  nrounds <- params$nrounds
  params$nrounds <- NULL

  xgb <- xgboost::xgboost(data = X,
                          label = y - 1, # xgboost requires labels to be in {0, ..., num_class - 1}.
                          params = params,
                          nrounds = nrounds,
                          verbose = 0,
                          nthread = 2,
                          objective = "multi:softprob",
                          eval_metric = "mlogloss",
                          num_class = num_class)

  fit <- function(X, n = NULL){
    pred <- stats::predict(xgb, newdata = as.matrix(X))
    return(matrix(pred, nrow = dim(X)[1], byrow = TRUE))
  }

  return(list("fit" = fit))

}

crossfit <- function(data,
                     xnum_class = NULL,
                     ynum_class = NULL,
                     method = "xgb",
                     xparams,
                     yparams,
                     nfolds = 5,
                     verbose = FALSE){

  # Check inputs

  if (nfolds < 1){
    warning("nfolds must be a positive integer. Setting nfolds = 1.")
    nfolds <- 1L
  }

  if (round(nfolds) != nfolds){
    warning("nfolds must be an integer. Rounding.")
    nfolds <- round(nfolds)
  }

  # Randomise folds

  n <- length(data$y)
  foldid <- sample(rep(seq(nfolds), length.out = n))

  f_pred <- matrix(NA, nrow = n, ncol = xnum_class)
  g_pred <- matrix(NA, nrow = n, ncol = ynum_class)
  out <- list()
  if (verbose){ out[["foldid"]] <- foldid }

  for (fold in seq(nfolds)){

    te <- (foldid == fold)

    if (nfolds > 1){
      tr <- !te
    } else{
      tr <- te
    }

    if (is.null(data$z)) { trz <- tez <- NULL } else {
      trz <- matrix(data$z[tr,], ncol = ncol(data$z))
      tez <- matrix(data$z[te,], ncol = ncol(data$z))
    }

    train <- list("z" = trz,
                  "x" = data$x[tr],
                  "y" = data$y[tr])
    test <- list("z" = tez,
                 "x" = data$x[te],
                 "y" = data$y[te])

    if (method == "xgb"){
      f_fit <- fit_xgboost(X = train$z,
                           y = train$x,
                           num_class = xnum_class,
                           params = xparams)$fit
      g_fit <- fit_xgboost(X = train$z,
                           y = train$y,
                           num_class = ynum_class,
                           params = yparams)$fit
    }

    if (verbose){
      out[[paste0("f_fit", fold)]] <- f_fit
      out[[paste0("g_fit", fold)]] <- g_fit
    }

    f_pred[te, ] <- f_fit(X = test$z, n = length(test$x))
    g_pred[te, ] <- g_fit(X = test$z, n = length(test$x))

  } # End of fold

  out[["f"]] <- f_pred
  out[["g"]] <- g_pred

  return(out)

}
