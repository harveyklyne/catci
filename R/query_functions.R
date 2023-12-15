greedy_query <- function(T_vector, Sigma, dx, dy,
                         metric, xsearch, ysearch, colsample_bylevel,
                         trees = list(NULL, NULL)){
  #' Compute metrics after successive label merging.
  #'
  #' Merge pairs of labels for each of X and Y so as to maximise the
  #' metric at each stage.
  #'
  #' @param T_vector vector of length dx*dy.
  #' @param Sigma matrix of dimension dx*dy by dx*dy.
  #' @param dx integer dimension of X.
  #' @param dy integer dimension of Y.
  #' @param metric str "approx_chi" or "norm" controlling the metric used for searching.
  #' @param xsearch str "greedy" or "ordinal" controlling the X search procedure.
  #' @param ysearch str "greedy" or "ordinal" controlling the Y search procedure.
  #' @param colsample_bylevel float in (0,1]. The subsample ratio of labels considered for merging at each level.
  #'     Subsampling occurs once for every new depth level reached.
  #' @param trees list containing tree structures (only used if xsearch or ysearch == tree)
  #'
  #' @return list containing vector of metrics and the corresponding merged labels.

  metric_list <- metric_lookup(metric)
  metric_fun <- metric_list$fun
  metric_args <- compute_metric_args(T_vector, Sigma, metric_list$args)

  all_metrics <- metric_fun(metric_args)

  categories <- list("x"=as.list(seq(dx)), "y"=as.list(seq(dy)))
  all_categories <- list(categories)

  dx_dy <- c(dx, dy)

  while(any(dx_dy > 2)){
    merge_out <- merge(T_vector = T_vector,
                       Sigma = Sigma,
                       dx = dx_dy[1],
                       dy = dx_dy[2],
                       metric = metric,
                       xsearch = xsearch,
                       ysearch = ysearch,
                       colsample_bylevel = colsample_bylevel,
                       trees = trees,
                       categories = categories)

    max_metric <- merge_out[which.max(merge_out[,"value"]),]
    # print(max_metric)
    all_metrics <- c(all_metrics, unname(max_metric["value"]))

    dimension <- unname(max_metric["dim"])
    lab1 <- unname(max_metric["ind1"])
    lab2 <- unname(max_metric["ind2"])
    index1 <- get_index(dimension = dimension,
                        j = lab1,
                        dx = dx_dy[1],
                        dy = dx_dy[2])
    index2 <- get_index(dimension = dimension,
                        j = lab2,
                        dx = dx_dy[1],
                        dy = dx_dy[2])

    categ <- c("x", "y")[dimension]
    categories[[categ]][[lab1]] <- c(categories[[categ]][[lab1]], categories[[categ]][[lab2]])
    categories[[categ]][[lab2]] <- NULL
    all_categories <- append(all_categories, list(categories))

    dx_dy[dimension] <- dx_dy[dimension] - 1
    T_vector <- update_T(T_vector, index1, index2)
    Sigma <- update_Sigma(Sigma, index1, index2)
  }

  return(list("values" = all_metrics,
              "categories" = all_categories))
}



query_lookup <- function(method){
  query <- NULL
  if (method == "ordinal"){
    query <- function(T_vector, Sigma, dx = NULL, dy = NULL){
      greedy_query(T_vector = T_vector,
                   Sigma = Sigma,
                   dx = dx,
                   dy = dy,
                   metric = "approx_chi",
                   xsearch = "ordinal",
                   ysearch = "ordinal",
                   colsample_bylevel = 1,
                   trees = list(NULL, NULL))$values}
  }
  if (method == "tree"){
    query <- function(T_vector, Sigma, dx = NULL, dy = NULL){
      greedy_query(T_vector = T_vector,
                   Sigma = Sigma,
                   dx = dx,
                   dy = dy,
                   metric = "approx_chi",
                   xsearch = "ordinal",
                   ysearch = "ordinal",
                   colsample_bylevel = 1,
                   trees = list(make_binary_tree(dx), make_binary_tree(dy)))$values}
  }
  if (method == "mGCM"){
    query <- function(T_vector, Sigma, dx = NULL, dy = NULL){
      max(abs(T_vector / sqrt(diag(Sigma))))
    }
  }
  if (method == "max"){
    query <- function(T_vector, Sigma, dx = NULL, dy = NULL){
      max(abs(T_vector))
    }
  }
  if (method == "euclid"){
    query <- function(T_vector, Sigma, dx = NULL, dy = NULL){
      sqrt(sum(T_vector^2))
    }
  }

  if (is.null(query)){stop("Method not found.")} else{ return(query) }
}
