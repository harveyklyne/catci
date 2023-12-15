bootstrap_T <- function(Sigma, n_boot){
  #' Generate bootstrap samples T ~ N(0, Sigma)
  #'
  #' @param Sigma covariance matrix with dimension dx*dy by dx*dy.
  #' @param n_boot integer number of bootstrap samples.
  #'
  #' @return dx*dy by n_boot matrix, each column of which is distributed as N(0,Sigma).

  sqrt_Sigma <- matrix_sqrt(Sigma)
  rnk <- dim(sqrt_Sigma)[1]
  return(t(sqrt_Sigma) %*% matrix(stats::rnorm(n_boot * rnk), nrow = rnk, ncol=n_boot))
}

matrix_sqrt <- function(Sigma){
  sqrt_Sigma <- suppressWarnings(chol(Sigma, pivot=TRUE))
  pivot <- attr(sqrt_Sigma, "pivot")
  rnk <- attr(sqrt_Sigma, "rank")
  sqrt_Sigma <- sqrt_Sigma[seq_len(rnk), order(pivot)]
  if (rnk == 1){
    return(matrix(sqrt_Sigma, nrow = 1))
  } else {return(sqrt_Sigma)}
}


double_bootstrap_pvalue <- function(metrics, metrics_boot){
  if (!is.vector(metrics)){ metrics <- as.vector(metrics) }
  if (length(metrics) == 1){
    if (is.matrix(metrics_boot)){
      if (nrow(metrics_boot)!=1){
        stop("Incompatible dimensions.")
      } else { metrics_boot <- as.vector(metrics_boot) }
    } else{
      return(1 - bootstrap_cdf(metrics, metrics_boot))
    }
  } else{

    max_q <- max_boot_cdf(metrics, metrics_boot)
    q_boot <- t(apply(metrics_boot, 1, function(x){
      (rank(x, ties.method = "random") - stats::runif(length(x))) / length(x)
    }))
    max_q_boot <- apply(q_boot, 2, max)
    return(1 - bootstrap_cdf(max_q, max_q_boot))
  }
}

max_boot_cdf <- function(metrics, metrics_boot){
  if (length(metrics) != nrow(metrics_boot)){
    stop("Incompatible dimensions.")
  }
  boot_cdfs <- sapply(X = seq_len(length(metrics)),
                      FUN = function(j){
                        bootstrap_cdf(metrics[j], metrics_boot[j,])
                      })
  return(max(boot_cdfs))
}

bootstrap_cdf <- function(x, x_boot){
  if (!is.vector(x)){ x <- as.vector(x) }
  if (length(x) != 1){
    stop("Incompatible dimensions.")
  }
  B <- length(x_boot)
  r <- sum(x_boot < x) + stats::runif(1) * (1 + sum(x_boot == x))
  return(r/(B+1))
}


