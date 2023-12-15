ankan <- function(x, y, f, g){
  #' Compute the p-value for conditional independence using the method of Ankan and Textor 2022.
  #'
  #' @param x integer vector of length n taking values in {1, ..., dx}.
  #' @param y integer vector of length n taking values in {1, ..., dy}.
  #' @param f matrix of dimension n by dx. ith row is the estimate of (Pr(X=1|Zi), ..., Pr(X=dx|Zi)).
  #' @param g matrix of dimension n by dy. ith row is the estimate of (Pr(Y=1|Zi), ..., Pr(Y=dy|Zi)).
  #'
  #' @return float p-value.

  n <- length(x)
  dx <- dim(f)[2]
  dy <- dim(g)[2]
  if (any(round(c(x,y)) != c(x,y))){stop("x and y must be integer vectors.")}
  if (length(y) != n){stop("x and y must be of the same length.")}
  if ((dim(f)[1] != n) | (dim(g)[1] != n)){stop("f and g must have the same number of rows as x and y.")}
  if ((min(c(x,y)) < 1) | (max(x) > dx) | (max(y) > dy)){stop("Must have 1 <= x <= dx and 1 <= y <= dy.")}

  rx <- Li_residual(x, f)
  ry <- Li_residual(y, g)

  stat <- 1/sqrt(n) * sum(rx * ry) / stats::sd(rx * ry) # asymptotically standard normal under H0
  return(1 - stats::pchisq(stat^2, df=1))
}

Li_residual <- function(x, f){
  #' Compute the residual of Li and Shepherd 2012.
  #'
  #' @param x integer vector of length n taking values in {1, ..., dx}.
  #' @param f matrix of dimension n by dx. ith row is the estimate of (Pr(X=1|Zi), ..., Pr(X=dx|Zi)).
  #'
  #' @return vector of residuals.

  n <- length(x)
  dx <- dim(f)[2]
  xcdf <- t(apply(f, 1, cumsum))
  xcdf <- cbind(0, xcdf[,1:(dx-1)], 1)
  # ij entry is Pr(Xi<j|Zi)
  # 1 - ij entry is Pr(Xi >= j|Zi) = Pr(Xi > j-1|Zi)
  return(xcdf[cbind(1:n, x)] - (1 - xcdf[cbind(1:n, x+1)]))
}

multinomial <- function(x, y, z){
  #' Compute the p-value for conditional independence using nested multinomial models.
  #'
  #' @param x integer vector of length n taking values in {1, ..., dx}.
  #' @param y integer vector of length n taking values in {1, ..., dy}.
  #' @param z matrix with n rows.
  #'
  #' @return float p-value.

  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("Package \"nnet\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  if (is.null(z)){
    mm1 <- nnet::multinom(y ~ 1, trace = FALSE)
    # mm2 <- nnet::multinom(y ~ x, trace = FALSE)
    mm2 <- nnet::multinom(y ~ as.factor(x), trace = FALSE)
  } else{
    mm1 <- nnet::multinom(y ~ z, trace = FALSE)
    # mm2 <- nnet::multinom(y ~ x + z, trace = FALSE)
    mm2 <- nnet::multinom(y ~ as.factor(x) + z, trace = FALSE)
  }
  return(stats::anova(mm1, mm2)$`Pr(Chi)`[2])
}

chi_sq <- function(T_vector, Sigma, dx, dy){
  #' Evaluate the chi-squared p-value for T_vector given Sigma.
  #'
  #' Under the null T is asymptotically N(0, Sigma), and so T^T %*% Sigma^+ %*% T
  #' follows a chi-square distribution with (dx-1)(dy-1) degrees of freedom, where
  #' Sigma^+ is the pseudoinverse of Sigma.
  #'
  #' @param T_vector vector of length dx*dy.
  #' @param Sigma matrix of dimension dx*dy by dx*dy.
  #' @param dx integer dimension of X.
  #' @param dy integer dimension of Y.
  #'
  #' @return float between 0 and 1

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  return(1 - stats::pchisq(as.numeric(t(T_vector) %*% MASS::ginv(Sigma) %*% T_vector), df=(dx-1)*(dy-1)))
}
