simulate_data <- function(n, xnum_class, ynum_class,
                          xsetting, ysetting,
                          strength, intsetting, permute){
  #' Generate an i.i.d. (X,Y,Z) dataset.
  #'
  #' @param n integer number of samples.
  #' @param xnum_class integer number of classes for discrete variable X.
  #' @param ynum_class integer number of classes for discrete variable Y.
  #' @param xsetting string determining X | Z setting.
  #' @param ysetting string determining Y | Z setting.
  #' @param strength float between 0 and 1 controlling the conditional dependence of Y on X given Z.
  #' @param intsetting string determining conditional dependence setting.
  #' @param permute logical determining if X, Y labels should be randomly permuted (violating ordinal structure)
  #'
  #' @return list containing (X, Y, Z) samples and potentially also corresponding population (f(Z), g(Z)).
  #' @export

  z <- z_correlated_normal(n = n, p = 5, corr = 0.5)

  if ((xsetting %in% c("lin", "vee", "hat")) & (ysetting %in% c("lin", "vee", "hat"))){

    qu <- 1/(1+exp(-3*z[,1]))
    qv <- 1/(1+exp(-3*z[,2]))

    u <- stats::rbinom(n, 1, qu)
    v <- stats::rbinom(n, 1, qv)

    x <- y <- rep(NA, n)
    for (u_case in c(0,1)){
      for (v_case in c(0,1)){
        ind <- (u == u_case) & (v == v_case)
        num <- sum(ind)
        if (num > 0){
          x_pdf <- get_pdf(setting = xsetting, d = xnum_class, z = NULL, modify = (u_case == 1))
          y_pdf <- get_pdf(setting = ysetting, d = ynum_class, z = NULL, modify = (v_case == 1))
          indep_pdf <- outer(x_pdf, y_pdf)
          dep_pdf <- indep_pdf + (strength/(8*xnum_class*ynum_class)) * get_int(intsetting, xnum_class, ynum_class)
          xy <- apply(stats::rmultinom(num, 1, dep_pdf), 2, function(x){which(x == 1)})
          x[ind] <- (xy - 1) %% xnum_class + 1
          y[ind] <- (xy - 1) %/% xnum_class + 1
        }
      }
    }

    f <- outer(1-qu, get_pdf(setting = xsetting, d = xnum_class, z = NULL, modify = FALSE)) + outer(qu, get_pdf(setting = xsetting, d = xnum_class, z = NULL, modify = TRUE))
    g <- outer(1-qv, get_pdf(setting = ysetting, d = ynum_class, z = NULL, modify = FALSE)) + outer(qv, get_pdf(setting = ysetting, d = ynum_class, z = NULL, modify = TRUE))
  }

  if (xsetting %in% c("sin", "sig")){
    f <- get_pdf(setting = xsetting, d = xnum_class, z = z[,1], modify = FALSE)
  }
  if (ysetting %in% c("sin", "sig")){
    g <- get_pdf(setting = ysetting, d = ynum_class, z = z[,2], modify = FALSE)
  }
  indep_pdf <- t(sapply(seq_len(n), function(i){as.vector(outer(f[i,], g[i,]))}))
  dep_pdf <- indep_pdf + (strength * min(indep_pdf) / 2) * matrix(as.vector(get_int(intsetting, xnum_class, ynum_class)), nrow = n, ncol = xnum_class*ynum_class, byrow = TRUE)

  xy <- sapply(seq_len(n), function(i){which(stats::rmultinom(1, 1, dep_pdf[i,]) == 1)})
  x <- (xy - 1) %% xnum_class + 1
  y <- (xy - 1) %/% xnum_class + 1


  if (any(is.na(c(x,y,as.vector(f), as.vector(g))))){stop("Data not generated correctly.")}

  if (permute){
    xperm <- permute_labels(x, f)
    x <- xperm$x
    f <- xperm$f
    yperm <- permute_labels(y, g)
    y <- yperm$x
    g <- yperm$f
  }

  return(list("x" = x,
              "y" = y,
              "z" = z,
              "f" = f,
              "g" = g
  ))

}


z_correlated_normal <- function(n, p, corr){
  #' Generate Z jointly gaussian Z variables.
  #'
  #' Z ~ N_{p}(0,Sigma),
  #' where Sigma_{jj} = 1, Sigma_{jk} = corr for all j not equal to k.
  #'
  #' @param n integer number of samples.
  #' @param p integer number of dimensions.
  #' @param corr float correlation in (-1,1).
  #'
  #' @return n by p matrix.

  Sigma <- matrix(corr, nrow=p, ncol=p)
  diag(Sigma) <- 1
  sqrt_Sigma <- chol(Sigma)
  return(matrix(stats::rnorm(n*(p)), nrow=n) %*% sqrt_Sigma)
}

lin_pdf <- function(d){
  pdf_unscale <- seq(1, 3, length.out = d)
  return(pdf_unscale / sum(pdf_unscale))
}

vee_pdf <- function(d){
  pdf_half <- seq(1, 3, length.out = ceiling(d/2))
  pdf_unscale <- c(pdf_half, rev(pdf_half[seq_len(floor(d/2))]))
  return(pdf_unscale / sum(pdf_unscale))
}

hat_pdf <- function(d){
  pdf_unscale <- rep(1, d)
  pdf_unscale[seq(floor(d/3) + 1, ceiling(2*d/3))] <- 3
  return(pdf_unscale / sum(pdf_unscale))
}

sin_pdf <- function(z, d){
  n <- length(z)
  f0 <- sinosoidal(z, a = 1)
  std <- 2
  if ((d == 8) & (std == 2)) {
    br <- c(-4.53, -2.30, -1.20, -0.60,  0.00,  0.60,  1.20, 2.30, 4.53)
  } else{stop("sin_pdf not defined")}
  cdf <- (matrix(br[2:d], nrow = n, ncol = d - 1, byrow = TRUE) - matrix(f0, nrow = n, ncol = d - 1, byrow = FALSE)) / std
  cdf <- stats::pnorm(cdf)
  cdf <- cbind(0, cdf, 1)
  return(cdf[,2:(d + 1)] - cdf[,1:d])
}

sinosoidal <- function(x,a){exp(-x^2/2)*sin(a*x)}

sig_pdf <- function(z, d){
  n <- length(z)
  f0 <- sigmoidal(z, s = 3)
  std <- 2
  if ((d == 8) & (std == 2)) {
    br <- c(-5.50, -2.00, -1.00,  0.00,  0.50,  1.00,  2.00, 3.00,  6.50)
  } else{stop("sig_pdf not defined")}
  cdf <- (matrix(br[2:d], nrow = n, ncol = d - 1, byrow = TRUE) - matrix(f0, nrow = n, ncol = d - 1, byrow = FALSE)) / std
  cdf <- stats::pnorm(cdf)
  cdf <- cbind(0, cdf, 1)
  return(cdf[,2:(d + 1)] - cdf[,1:d])
}

sigmoidal <- function(x,s){1/(1+exp(-s*x))}

get_pdf <- function(setting, d, z, modify){
  out <- NULL
  if(setting == "lin"){out <- lin_pdf(d)}
  if(setting == "vee"){out <- vee_pdf(d)}
  if(setting == "hat"){out <- hat_pdf(d)}
  if(setting == "sin"){out <- sin_pdf(z, d)}
  if(setting == "sig"){out <- sig_pdf(z, d)}
  if(is.null(out)){stop("Setting not recognised.")}


  if (modify){out <- 2/d - out}

  return(out)
}

get_int <- function(setting, xnum_class, ynum_class){
  out <- NULL
  if(setting == "step"){
    outcol <- c(rep(-1, floor(xnum_class/2)),
                rep(0, xnum_class %% 2),
                rep(1, floor(xnum_class/2)))
    outall <- c(rep(outcol, times = floor(ynum_class/4)),
                rep(rev(outcol), times = floor(ynum_class/4)),
                rep(0, xnum_class * (ynum_class %% 4)),
                rep(rev(outcol), times = floor(ynum_class/4)),
                rep(outcol, times = floor(ynum_class/4)))
    out <- matrix(outall, nrow = xnum_class, ncol = ynum_class)
  }
  if(setting == "alt"){
    outcol <- rep(c(-1,1), times = floor(xnum_class/2))
    if (xnum_class %% 2 == 1) {
      outcol <- c(outcol[1:floor(xnum_class/2)], 0, outcol[-(1:floor(xnum_class/2))])
    }
    outall <- rep(c(outcol, -outcol), times = floor(ynum_class/2))
    if (ynum_class %% 2 == 1) {
      outall <- c(outall[1:(xnum_class*floor(ynum_class/2))], rep(0, xnum_class), outall[-(1:(xnum_class*floor(ynum_class/2)))])
    }
    out <- matrix(outall, nrow = xnum_class, ncol = ynum_class)
  }
  if (setting == "binary_tree"){
    if ((xnum_class != 8) | (ynum_class != 8)){
      stop("Binary tree interaction only coded for xnum_class = ynum_class = 8.")
    }

    outall <- c(-1.1,  -1.1,  -0.9,  -0.9,  0.9,   0.9,   1.1,   1.1,   -1.1,  -1.1,  -0.9,
                -0.9,  0.9,   0.9,   1.1,   1.1,   -0.9,  -0.9,  -1.1,  -1.1,  1.1,   1.1,
                0.9,   0.9,   -0.9,  -0.9,  -1.1,  -1.1,  1.1,   1.1,   0.9,   0.9,   0.9,
                0.9,   1.1,   1.1,   -1.1,  -1.1,  -0.9,  -0.9,  0.9,   0.9,   1.1,   1.1,
                -1.1,  -1.1,  -0.9,  -0.9,  1.1,   1.1,   0.9,   0.9,   -0.9,  -0.9,  -1.1,
                -1.1,  1.1,   1.1,   0.9,   0.9,   -0.9,  -0.9,  -1.1,  -1.1)
    out <- matrix(outall, nrow = xnum_class, ncol = ynum_class)
  }
  if(is.null(out)){stop("Setting not recognised.")}
  return(out)
}

permute_labels <- function(x, f){
  xnum_class <- dim(f)[2]
  perm <- sample(seq_len(xnum_class))
  ord <- order(perm)
  x <- perm[x]
  f <- f[, ord]
  return(list("x" = x,
              "f" = f))
}

formulate_statistics <- function(data,
                                 xnum_class = NULL,
                                 ynum_class = NULL,
                                 method = "xgb",
                                 xparams,
                                 yparams,
                                 nfolds = 5,
                                 normalise){

  pred <- crossfit(data = data,
                   xnum_class = xnum_class,
                   ynum_class = ynum_class,
                   method = method,
                   xparams = xparams,
                   yparams = yparams,
                   nfolds = nfolds,
                   verbose = FALSE)

  out <- form_T_Sigma(x = data$x,
                      y = data$y,
                      f = pred$f,
                      g = pred$g,
                      normalise = normalise)


  # Multinomial Y ~ (X,Z)
  out[["multinomial"]] <- multinomial(x = data$x,
                                      y = data$y,
                                      z = data$z)

  # Ankan and Textor
  out[["ankan"]] <- ankan(x = data$x,
                          y = data$y,
                          f = pred$f,
                          g = pred$g)

  # Chi square test using pseudo inverse
  out$chi <- chi_sq(T_vector = out$T_vector,
                    Sigma = out$Sigma,
                    dx = xnum_class,
                    dy = ynum_class)

  return(out)
}

form_T_Sigma <- function(x, y, f, g, normalise){
  #' Form T_vector and Sigma from observed labels x and y and probability estimate
  #' matrices f and g.
  #'
  #' @param x integer vector of length n taking values in {1, ..., dx}.
  #' @param y integer vector of length n taking values in {1, ..., dy}.
  #' @param f matrix of dimension n by dx. ith row is the estimate of (Pr(X=1|Zi), ..., Pr(X=dx|Zi)).
  #' @param g matrix of dimension n by dy. ith row is the estimate of (Pr(Y=1|Zi), ..., Pr(Y=dy|Zi)).
  #' @param normalise logical determining if normalised or unnormalised generalised covariance is to be used.
  #'
  #' @return list containing T_vector and Sigma matrix.

  n <- length(x)
  dx <- dim(f)[2]
  dy <- dim(g)[2]
  if (any(round(c(x,y)) != c(x,y))){stop("x and y must be integer vectors.")}
  if (length(y) != n){stop("x and y must be of the same length.")}
  if ((dim(f)[1] != n) | (dim(g)[1] != n)){stop("f and g must have the same number of rows as x and y.")}
  if ((min(c(x,y)) < 1) | (max(x) > dx) | (max(y) > dy)){stop("Must have 1 <= x <= dx and 1 <= y <= dy.")}


  # convert integer vectors to design matrices
  X <- stats::model.matrix(~ factor(x, levels = seq(dx)) -1 ) - f
  Y <- stats::model.matrix(~ factor(y, levels = seq(dy)) -1 ) - g

  prod_mat <- matrix(NA, nrow = n, ncol = dx * dy)
  for (row in seq(n)){
    prod_mat[row, ] <- Y[row,] %x% X[row,] # (xij-E(Xj|Zi))(yik-E(Yk|Zi))
  }
  T_vector <- sqrt(n)*colMeans(prod_mat)
  Sigma <- stats::var(prod_mat) # Consider computing using structure under the null.

  if (normalise){
    T_scale <- sqrt(diag(Sigma))
    T_vector <- T_vector / T_scale
    Sigma <- Sigma / (T_scale %o% T_scale)
  }

  return(list("T_vector" = T_vector, "Sigma" = Sigma))

}


evaluate_sim <- function(data, dx, dy, n_boot, methods){
  #' Wrapper function to evaluate simulation data.
  #'
  #' @param data list containing T_vector and Sigma.
  #' @param dx integer number of x labels.
  #' @param dy integer number of y labels.
  #' @param n_boot integer number of bootstrap samples.
  #' @param methods vector of strings specifying the query functions.
  #'
  #' @return vector of p-values corresponding to different methods.
  #' @export

  # Save p-values computed during data simulation
  out <- data
  out$T_vector <- NULL
  out$Sigma <- NULL

  T_all <- cbind(data$T_vector, bootstrap_T(data$Sigma, n_boot))

  for (method in methods){

    query <- query_lookup(method)
    metrics_all <- matrix(apply(T_all, 2, function(T_vec){
      query(T_vector = T_vec,
            Sigma = data$Sigma,
            dx = dx,
            dy = dy)
    }), ncol = 1 + n_boot)

    out[[method]] <- double_bootstrap_pvalue(metrics_all[,1],
                                             metrics_all[,-1])
  }

  return(out)

}


simulate_train_test <- function(n_tr, n_te,
                                xnum_class,
                                ynum_class,
                                xsetting, ysetting,
                                strength,
                                intsetting){
  data <- simulate_data(n = n_tr + n_te,
                        xnum_class = xnum_class,
                        ynum_class = ynum_class,
                        xsetting = xsetting,
                        ysetting = ysetting,
                        strength = strength,
                        intsetting = intsetting,
                        permute = FALSE)

  smpl <- sample(length(data$x), size = n_tr, replace = FALSE)
  train <- list("x" = data$x[smpl],
                "y" = data$y[smpl],
                "z" = data$z[smpl,])
  test <- list("x" = data$x[-smpl],
               "y" = data$y[-smpl],
               "z" = data$z[-smpl,])

  return(list("train" = train,
              "test" = test))

}

