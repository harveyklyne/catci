merge <- function(T_vector, Sigma, dx, dy,
                  metric, xsearch, ysearch, colsample_bylevel,
                  trees = list(NULL, NULL), categories= list(NULL, NULL)){
  #' Compute p-values after merging each possible pair of X and Y labels.
  #'
  #' X-merges are only considered if dx>2, similarly for Y.
  #'
  #' @param T_vector vector of length dx*dy.
  #' @param Sigma matrix of dimension dx*dy by dx*dy.
  #' @param dx integer dimension of X.
  #' @param dy integer dimension of Y.
  #' @param metric str controlling the metric used for searching.
  #' @param xsearch str "greedy", "ordinal", or "tree" controlling the X search procedure.
  #' @param ysearch str "greedy", "ordinal", or "tree" controlling the Y search procedure.
  #' @param colsample_bylevel float in (0,1]. The subsample ratio of potential merges to be considered.
  #'     Higher values result in more regularisation.
  #' @param trees list containing tree structures (only used if xsearch or ysearch == tree)
  #' @param categories list containing current merging structures (only used if xsearch or ysearch == tree)
  #'
  #' @return matrix with 4 named columns. "dim" equals 1 for X-label merges and 2 for Y-label merges.
  #' ind1 < ind2 correspond to the indices of the labels being merged. "value"
  #' is a measure of belief in the null hypothesis of conditional independence, with smaller
  #' values favouring rejection.

  num_rows <- get_num_levels(xsearch, dx, trees[[1]], categories[[1]]) + get_num_levels(ysearch, dy, trees[[2]], categories[[2]])

  out <- matrix(NA, nrow = num_rows, ncol=4)
  colnames(out) = list("dim","ind1","ind2","value")

  metric_list <- metric_lookup(metric)
  metric_fun <- metric_list$fun
  metric_args <- compute_metric_args(T_vector, Sigma, metric_list$args)

  row_num <- 0 # count row number
  permitted_rows <- sort(sample(num_rows, ceiling(num_rows * colsample_bylevel), replace = FALSE))

  # Consider possible merges
  for (dimension in c(1,2)){ # dimension 1 is X, dimension 2 is Y
    d <- c(dx, dy)[dimension]
    search <- c(xsearch, ysearch)[dimension]

    if (d > 2){
      for (j1 in seq(d - 1)){
        index1 <- get_index(dimension, j1, dx, dy)

        for (j2 in get_ind2(ind1 = j1, d = d, search = search,
                            tree = trees[[dimension]],
                            category = categories[[dimension]])) {
          row_num <- row_num + 1

          if (row_num %in% permitted_rows){
            index2 <- get_index(dimension, j2, dx, dy)

            new_metric_args <- update_metric_args(T_vector = T_vector,
                                                  Sigma = Sigma,
                                                  index1 = index1,
                                                  index2 = index2,
                                                  current_args = metric_args)

            value <- metric_fun(new_metric_args)
          } else { value <- NA }

          out[row_num,] <- c(dimension, j1, j2, value)
        }
      }
    }
  }

  return(out)
}


get_index <- function(dimension, j, dx, dy){
  if (dimension == 1){ return(X_index(j, dx, dy)) }
  if (dimension == 2){ return(Y_index(j, dx, dy)) }
}

X_index <- function(j, dx, dy){
  #' Return the indices whose x-entry corresponds to element j.
  #'
  #' For a vector of length dx*dy indexed (1,1), ..., (dx, 1), (1, 2), ...,
  #' (dx, 2), ..., (1, dy), ..., (dx, dy), this function returns the indices
  #' corresponding to (j,k) for k = 1, ..., dy.
  #'
  #' @param j integer index 1 <= j <= dx.
  #' @param dx integer dimension of X.
  #' @param dy integer dimension of Y.
  #'
  #' @return boolean vector of length dx*dy with kth TRUE entry corresponding
  #' to index (j,k).

  if ((j<1) | (j>dx)){stop("Need 1 <= j <= dx.")}

  return((seq(dx*dy) %% dx) == (j %% dx))

}

Y_index <- function(k, dx, dy){
  #' Return the indices whose y-entry corresponds to element k.
  #'
  #' For a vector of length dx*dy indexed (1,1), ..., (dx, 1), (1, 2), ...,
  #' (dx, 2), ..., (1, dy), ..., (dx, dy), this function returns the indices
  #' corresponding to (j,k) for j = 1, ..., dx.
  #'
  #' @param k integer index 1 <= k <= dy.
  #' @param dx integer dimension of X.
  #' @param dy integer dimension of Y.
  #'
  #' @return boolean vector of length dx*dy with jth TRUE entry corresponding
  #' to index (j,k).

  if ((k<1) | (k>dy)){stop("Need 1 <= k <= dy.")}

  return((seq(0,dx*dy-1) %/% dx) == (k - 1))

}

update_T <- function(T_vector, index1, index2){
  #' Update T vector by merging pairs of indices.
  #'
  #' @param T_vector vector of length p.
  #' @param index1 boolean vector of length p, with d TRUE entries.
  #' @param index2 boolean vector of length p, with d TRUE entries.
  #'
  #' @return vector of length p-d corresponding to merging the d pairs of
  #' labels indexed by index1 and index2.

  if (sum(index1)!=sum(index2)){
    stop("Need index1 and index2 to have the name number of TRUE entries.")
  }

  new_T <- T_vector
  new_T[index1] <- T_vector[index1] + T_vector[index2]
  new_T <- new_T[!index2]

  return(new_T)
}

update_Sigma <- function(Sigma, index1, index2){
  #' Update Sigma matrix by merging pairs of indices.
  #'
  #' @param Sigma matrix of dimension p by p.
  #' @param index1 boolean vector of length p, with d TRUE entries.
  #' @param index2 boolean vector of length p, with d TRUE entries.
  #'
  #' @return matrix of dimension (p-d) by (p-d) corresponding to merging the d
  #' pairs of labels indexed by index1 and index2.

  if (is.null(Sigma)){ return(NULL) } else{

    if (sum(index1)!=sum(index2)){
      stop("Need index1 and index2 to have the name number of TRUE entries.")
    }

    new_Sigma <- Sigma
    new_Sigma[index1,] <- new_Sigma[index1, ] + new_Sigma[index2,]
    new_Sigma[,index1] <- new_Sigma[,index1] + new_Sigma[,index2]
    new_Sigma <- new_Sigma[!index2, !index2]

    return(new_Sigma)
  }
}

update_normsq <- function(normsq, T_vector, index1, index2){
  #' Update \|T\|^2_2 after merging pairs of indices.
  #'
  #' @param normsq float equal to \|T_vector\|^2_2.
  #' @param T_vector vector of length p.
  #' @param index1 boolean vector of length p, with d TRUE entries.
  #' @param index2 boolean vector of length p, with d TRUE entries.
  #'
  #' @return float equal to the squared norm of T after merging the pairs of
  #' labels.

  if (sum(index1)!=sum(index2)){
    stop("Need index1 and index2 to have the name number of TRUE entries.")
  }

  return(normsq + 2 * sum(T_vector[index1] * T_vector[index2]))
}

update_tr <- function(tr, Sigma, index1, index2){
  #' Update tr(Sigma) after merging pairs of indices.
  #'
  #' @param tr float equal to tr(Sigma).
  #' @param Sigma matrix of dimension p by p.
  #' @param index1 boolean vector of length p, with d TRUE entries.
  #' @param index2 boolean vector of length p, with d TRUE entries.
  #'
  #' @return float equal to the trace of Sigma after merging the pairs of
  #' labels.

  if (sum(index1)!=sum(index2)){
    stop("Need index1 and index2 to have the name number of TRUE entries.")
  }

  return(tr + 2*sum(diag(as.matrix(Sigma[index1,index2]))))
}

update_tr2 <- function(tr2, Sigma, index1, index2){
  #' Update tr(Sigma^2) after merging pairs of indices.
  #'
  #' @param tr2 float equal to tr(Sigma^2).
  #' @param Sigma matrix of dimension p by p.
  #' @param index1 boolean vector of length p, with d TRUE entries.
  #' @param index2 boolean vector of length p, with d TRUE entries.
  #'
  #' @return float equal to the trace of Sigma^2 after merging the pairs of
  #' labels.

  if (sum(index1)!=sum(index2)){
    stop("Need index1 and index2 to have the name number of TRUE entries.")
  }

  return(tr2 + 4 * sum(Sigma[,index1] * Sigma[,index2]) +
           2 * sum(Sigma[index1,index1] * Sigma[index2,index2] +
                     Sigma[index1,index2] * Sigma[index2,index1]))
}



get_ind2 <- function(ind1, d, search, tree, category){
  if (search == "greedy"){
    return(seq(from=ind1+1, to=d))
  }
  if (search == "ordinal"){
    return(ind1+1)
  }
  if (search == "tree"){
    return(get_tree_ind2(ind1 = ind1, tree = tree, category = category))
  }
}

get_num_levels <- function(method, d, tree, category){
  if (method == "greedy"){
    return((d > 2) * choose(d, 2))
  }
  if (method == "ordinal"){
    return((d > 2) * (d - 1))
  }
  if (method == "tree"){
    return(get_tree_levels(tree, category))
  }
}

metric_lookup <- function(metric){
  fun <- NULL
  if (metric == "approx_chi"){
    fun <- function(args){
      approx_chi_metric(normsq = args$normsq,
                        tr = args$tr,
                        tr2 = args$tr2)}
    args <- names(formals(approx_chi_metric))
  }
  # if (metric == "l2"){
  #   fun <- function(args){
  #     sum_sq_metric(normsq = args$normsq)}
  #   args <- names(formals(sum_sq_metric))
  # }

  if (is.null(fun)){
    stop("Metric not found.")
  } else{ return(list("fun" = fun,
                      "args" = args)) }
}


compute_metric_args <- function(T_vector, Sigma, args){
  inputs <- list()
  if ("normsq" %in% args){inputs["normsq"] <- sum(T_vector^2)}
  if ("tr" %in% args){inputs["tr"] <- sum(diag(Sigma))}
  if ("tr2" %in% args){inputs["tr2"] <- sum(Sigma^2)}
  return(inputs)
}

update_metric_args <- function(T_vector, Sigma, index1, index2, current_args){
  args <- names(current_args)
  new_inputs <- list()
  if ("normsq" %in% args){
    new_inputs["normsq"] <- update_normsq(normsq = current_args$normsq,
                                          T_vector = T_vector,
                                          index1 = index1,
                                          index2 = index2)
  }
  if ("tr" %in% args){
    new_inputs["tr"] <- update_tr(tr = current_args$tr,
                                  Sigma = Sigma,
                                  index1 = index1,
                                  index2 = index2)
  }
  if ("tr2" %in% args){
    new_inputs["tr2"] <- update_tr2(tr2 = current_args$tr2,
                                    Sigma = Sigma,
                                    index1 = index1,
                                    index2 = index2)
  }
  return(new_inputs)
}


approx_chi_metric <- function(normsq, tr, tr2){
  #' Evaluate the approximate chi-squared cdf using the method of Box (1954).
  #'
  #' Under the null T is asymptotically N(0, Sigma), and \|T\|_2^2 can be
  #' approximated by a rescaled chi-square with scaling and degrees of freedom
  #' determined by Sigma.
  #'
  #' @param normsq float equal to \|T_vector\|^2_2.
  #' @param tr float equal to tr(Sigma).
  #' @param tr2 float equal to tr(Sigma^2).
  #'
  #' @return float between 0 and 1

  g <- tr2 / tr
  h <- tr^2 / tr2

  return(stats::pchisq(normsq/g, df=h))

}

# sum_sq_metric <- function(normsq, ...){normsq}
