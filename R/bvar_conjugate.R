# todo: 
# make parallel computations in estimate 

#' Build Y matrix from supplied data
#' 
#' Build Y matrix from supplied data
#' 
#' Build Y matrix from supplied data
#'
#' @param p number of lags
#' @param Y_in multivariate time series [T_in x m]
#' @return Y [T x m] matrix of left hand side endogeneous variables, T = T_in - p
#' @export
#' @examples 
#' data(Yraw)
#' Y <- bvar_build_Y(Yraw, p=4)
bvar_build_Y <- function(Y_in, p=1) {
  
  # transform Y_in from data.frame to matrix
  # so that cbind(NULL, Y_block) will work 
  Y_in <- as.matrix(Y_in)
  
  # p first observations of Y_in are used only as regressors
  Y <- tail(Y_in, -p)
  
  return(Y)
}

#' Build X matrix from supplied data
#' 
#' Build X matrix from supplied data
#' 
#' Build X matrix from supplied data
#' 
#' @param p number of lags
#' @param Y_in endogeneous variables [T_in x m]
#' @param Z_in optional exogeneous variables [T_in x ...]
#' @param constant whether the constant should be included, default is TRUE
#' @return X [T x k] matrix of right hand side regressors: endogeneous and exogeneous variables,
#' T = T_in - p, k = m*p + d, d is the number of exogeneous variables including constant if present
#' @export
#' @examples 
#' data(Yraw)
#' X <- bvar_build_X(Yraw, constant=TRUE, p=4)
bvar_build_X <- function(Y_in, Z_in=NULL, constant=TRUE, p=1) {
  
  # if Z_in is provided it should have the same number of rows that Y_in
  if (!is.null(Z_in)) 
    if (!nrow(Y_in)==nrow(Z_in))
      stop("Number of rows in Y_in and Z_in should be equal. 
           The first p rows of Z_in are not used and may be filled with NA")
  
  # if requested add constant to exogeneous regressors
  if (constant) {
    Z_in <- cbind(rep(1, nrow(Y_in)), Z_in)
    colnames(Z_in)[1] <- "const"
  }
  
  # transform Y_in from data.frame to matrix
  # so that cbind(NULL, Y_block) will work 
  Y_in <- as.matrix(Y_in)
  
  # p first observations of Y_in are used only as regressors
  Y <- tail(Y_in, -p)
  # p first observations of Z_in are not used only at all
  Z <- tail(Z_in, -p)
  
  # number of observations supplied
  T_in <- nrow(Y_in)
  
  # create X matrix
  X <- NULL
  for (j in 1:p) {
    Y_block <- Y_in[(p+1-j):(T_in-j),]
    colnames(Y_block) <- paste0(colnames(Y), ", l=", j)
    X <- cbind(X, Y_block) 
  }
  X <- cbind(X,Z)
  # X = [endo | const | other exo]
  
  return(X)
}

#' Recover original Y_in from Y, X and number of lags p
#' 
#' Recover original Y_in from Y, X and number of lags p
#' 
#' Recover original Y_in from Y, X and number of lags p
#' 
#' @param p number of lags
#' @param X [T x k] matrix of right hand side regressors: endogeneous and exogeneous variables
#' @param Y [T x m] matrix of left hand side endogeneous variables
#' @return original supplied Y_in
#' @export
#' @examples 
#' data(Yraw)
#' Y <- bvar_build_Y(Yraw, p=4)
#' X <- bvar_build_X(Yraw, constant=TRUE, p=4)
#' Y_in <- bvar_get_Y_in(Y,X,p=4)
bvar_get_Y_in <- function(Y,X,p) {
  # Y contains all obs except the first p
  m <- ncol(Y)
  Y_in <- rbind(X[1:p, (m*(p-1)+1):(m*p)] , Y)
  colnames(Y_in) <- colnames(Y)
  rownames(Y_in) <- NULL
  return(Y_in)
}


#' Calculate hyperparameters from artificial observations  
#' 
#' Calculate hyperparameters from artificial observations 
#' 
#' Calculate hyperparameters from artificial observations. Can be used both for prior and posterior
#' hyperparameters. In the case of prior hyperparameters X_star usually consists of X_cniw, X_io and X_sc. 
#' In the case of posterior hyperparameters X_star usually consists of X, X_cniw, X_io and X_sc. 
#' 
#' @param X_star [T_star x k] matrix of right hand side regressors: endogeneous and exogeneous variables
#' @param Y_star [T_star x m] matrix of left hand side endogeneous variables
#' @return list of hyperparameters: Omega, Omega_root, Phi, S
#' @export
#' @examples 
#' data(Yraw)
#' X <- bvar_build_X(Yraw, constant=TRUE, p=4)
#' Y <- bvar_build_Y(Yraw, p=4)
#' post <- bvar_conj_dummy2hyper(Y,X) # no dummy observation, just for demo
bvar_conj_dummy2hyper <- function(Y_star, X_star) {
  
  X_star_svd <- svd(X_star) # risky operation
  
  diag_inv <- diag(1/X_star_svd$d)
  
  if (all(is.finite(diag_inv))) {
    Omega_root <- X_star_svd$v %*% diag_inv %*% t(X_star_svd$v)
    
    # following operations are riskless (no inverse, no root, + and * only)
    Omega <- Omega_root %*% Omega_root
    
    Phi_star <- Omega %*% crossprod(X_star, Y_star)
    E_star <- Y_star - X_star %*% Phi_star
    S <- crossprod(E_star)
    
    endo_varnames <- bvar_get_endo_varnames(Y_star)
    X_varnames <- colnames(X_star)
  } else { # we'll do our best to deal with infinite values...
    
  }
  
  


  colnames(Omega) <- X_varnames
  rownames(Omega) <- X_varnames
  colnames(Omega_root) <- X_varnames
  rownames(Omega_root) <- X_varnames
  colnames(S) <- endo_varnames
  rownames(S) <- endo_varnames
  colnames(Phi_star) <- paste0("eq_",endo_varnames)
  rownames(Phi_star) <- X_varnames
  
  hyper <- list(Omega_root=Omega_root, Omega=Omega, Phi=Phi_star, S=S)
  return(hyper)
}


#' Calculate artificial observations from hyperparameters
#' 
#' Calculate artificial observations from hyperparameters
#' 
#' Calculate artificial observations from hyperparameters. 
#' Usually used to calculate artificial observations
#' from prior hyperparameters. 
#' 
#' @param Omega [k x k] hyperparameter, scale matrix for covariance of coefficients
#' @param S [m x m] hyperparameter, scale matrix for IW-distribution
#' @param Phi [k x m] hyperparameter, expected values of coefficients
#' @return list of artificial observations: X_plus, Y_plus
#' @export
#' @examples 
#' data(Yraw)
#' X <- bvar_build_X(Yraw, constant=TRUE, p=4)
#' Y <- bvar_build_Y(Yraw, p=4)
#' hyper <- bvar_conj_dummy2hyper(Y,X) # no dummy observation, just for demo
#' dummy <- bvar_conj_hyper2dummy(hyper$Omega, hyper$S, hyper$Phi)
bvar_conj_hyper2dummy <- function(Omega, S, Phi) {
  # http://math.stackexchange.com/questions/712993/cholesky-decomposition-of-the-inverse-of-a-matrix
  X_block_1 <- solve(chol(Omega)) # very risky operation: both chol and solve may fail 
  Y_block_1 <- solve(X_block_1, t(X_block_1) %*% X_block_1 %*% Phi)
  
  X_block_2 <- matrix(0, nrow=ncol(S), ncol=ncol(Omega))
  Y_block_2 <- chol(S)
  
  X_plus <- rbind(X_block_1, X_block_2)
  Y_plus <- rbind(Y_block_1, Y_block_2)
  
  rownames(X_plus) <- NULL
  rownames(Y_plus) <- NULL
  colnames(Y_plus) <- colnames(S)
  
  dummy <- list(X_plus=X_plus,Y_plus=Y_plus)
  return(dummy)
}


#' Create prior hyperparameters from lambdas 
#' 
#' Create prior hyperparameters from lambdas
#' 
#' Create prior hyperparameters from lambdas. 
#' Lambdas specification is based on Carriero p. 52-53
#'
#' @param p number of lags
#' @param Y_in multivariate time series
#' @param lambda vector = (l_1, l_lag, l_sc, l_io, l_const, l_exo), the l_kron is set to 1 automatically for 
#' conjugate N-IW prior. Short summary valid for NO sc/io case:
#' sd(const in eq i) = l_const * sigma_i
#' sd(exo in eq i)= l_exo * sigma_i
#' sd(coef for var j lag l in eq i) = l_1*sigma_i/sigma_j/l^l_lag
#' lambdas may be Inf
#' l_io or l_sc equal to NA means no corresponding dummy observations
#' @param Z_in exogeneous variables
#' @param constant logical, default is TRUE, whether the constant should be included
#' @param s2_lag number of lags in AR() model used to estimate s2 (equal to p by default)
#' Carriero uses 1 in his matlab code
#' @param delta vector [m x 1] or scalar or "AR1". Are used for prior Phi_1 and in sc/io dummy observations
#' Scalar value is replicated m times. If set to "AR1" then deltas will be estimated as AR(1) coefficients (but not greater than one).
#' Diagonal of Phi_1 is equal to delta. y_bar is multiplied by delta componentwise.
#' By default delta is equal to 1.
#' @param y_bar_type (either "all" or "initial"). Determines how y_bar for sc and io dummy is calculated.
#' "all": y_bar is mean of y for all observations, "initial": p initial observations
#' Carriero: all, Sim-Zha: initial
#' @param carriero_hack logical, if TRUE sigma^2 will be estimated using biased estimator
#' and supposed error with no square roots in dummy observations will be reproduced
#' FALSE by default
#' @return dummy list containing:
#' Omega, Omega_root
#' S, Phi
#' @export
#' @examples 
#' data(Yraw)
#' prior <- bvar_conj_lambda2hyper(Yraw, p = 4, lambda = c(0.2,1,1,1,100,100))
bvar_conj_lambda2hyper <- function(Y_in, Z_in=NULL, constant=TRUE, p=4, 
                                   lambda=c(0.2,1,1,1,100,100), 
                                   delta=1,
                                   s2_lag=NULL, 
                                   y_bar_type = c("initial", "all"),
                                   carriero_hack = FALSE) {
 
  
  prior <- list()
  return(prior)
}

#' Create [m x 1] vector of deltas from delta description
#' 
#' Create [m x 1] vector of deltas from delta description
#' 
#' Create [m x 1] vector of deltas from delta description
#' 
#' @param Y_in multivariate time series
#' @param delta description of delta. May be: one number, "AR1", or ready to go [m x 1] vector
#' @export
#' @return [m x 1] vector of deltas
#' @examples 
#' data(Yraw)
#' dummy <- bvar_conj_delta(Yraw, delta="AR1") 
bvar_conj_delta <- function(Y_in, delta) {
  m <- ncol(Y_in)
  if (delta[1] == "AR1") { # set deltas as AR(1) coefficients but no more than 1
    delta <- rep(1,m) # reserve space
    for (j in 1:m) {
      y_uni <- Y_in[,j] # univariate time series
      # more robust version: fails only in the case of  severe multicollinearity
      AR_1 <- ar.ols(y_uni, aic=FALSE, order.max = 1) # AR(1) model
      delta[j] <- AR_1$ar
      if (delta[j]>1) delta[j] <- 1
    }
  }
  if (length(delta)==1) delta <- rep(delta, m) # copy scalar values
  if (! length(delta)==m ) stop("Length of delta should be equal to 1 or m")
  return(delta)
}



#' Create [m x 1] vector of sigma^2 from supplied time series
#' 
#' Create [m x 1] vector of sigma^2 from supplied time series
#' 
#' Create [m x 1] vector of sigma^2 from supplied time series
#' 
#' @param Y_in multivariate time series
#' @param s2_lag number of lags to estimate sigma^2
#' @param carriero_hack logical, if TRUE sigma^2 will be estimated using biased estimator
#' and supposed error with no square roots in dummy observations will be reproduced
#' FALSE by default
#' @export
#' @return [m x 1] vector of sigmas^2
#' @examples 
#' data(Yraw)
#' sigma2 <- bvar_conj_sigma2(Yraw) 
bvar_conj_sigma2 <- function(Y_in, s2_lag = 1, carriero_hack = FALSE) {
  m <- ncol(Y_in)

  sigmas_sq <- rep(NA, m)
  for (j in 1:m) {
    y_uni <- Y_in[,j] # univariate time series
    # more robust version: fails only in the case of  severe multicollinearity
    AR_p <- ar.ols(y_uni, aic=FALSE, order.max = s2_lag) # AR(p) model
    resid <- tail(AR_p$resid,-s2_lag) # omit first p NA in residuals
    sigmas_sq[j] <- sum(resid^2)/(length(resid)-s2_lag-1)
    if (carriero_hack) { # Carriero uses biased sigma^2 estimator in his code
      sigmas_sq[j] <- sum(resid^2)/length(resid)
    }
  }
  
  return(sigmas_sq)
}




#' Create dummy observations from lambdas 
#' 
#' Create dummy observations from lambdas
#' 
#' Create dummy observations from lambdas. 
#' Lambdas specification is based on Carriero p. 52-53
#'
#' @param p number of lags
#' @param Y_in multivariate time series
#' @param lambda vector = (l_1, l_lag, l_sc, l_io, l_const, l_exo), the l_kron is set to 1 automatically for 
#' conjugate N-IW prior. Short summary valid for NO sc/io case:
#' sd(const in eq i) = l_const * sigma_i
#' sd(exo in eq i)= l_exo * sigma_i
#' sd(coef for var j lag l in eq i) = l_1*sigma_i/sigma_j/l^l_lag
#' lambdas may be Inf
#' l_io or l_sc equal to NA means no corresponding dummy observations
#' @param Z_in exogeneous variables
#' @param constant logical, default is TRUE, whether the constant should be included
#' @param s2_lag number of lags in AR() model used to estimate s2 (equal to p by default)
#' Carriero uses 1 in his matlab code
#' @param delta vector [m x 1] or scalar or "AR1". Are used for prior Phi_1 and in sc/io dummy observations
#' Scalar value is replicated m times. If set to "AR1" then deltas will be estimated as AR(1) coefficients (but not greater than one).
#' Diagonal of Phi_1 is equal to delta. y_bar is multiplied by delta componentwise.
#' By default delta is equal to 1.
#' @param y_bar_type (either "all" or "initial"). Determines how y_bar for sc and io dummy is calculated.
#' "all": y_bar is mean of y for all observations, "initial": p initial observations
#' Carriero: all, Sim-Zha: initial
#' @param carriero_hack logical, if TRUE sigma^2 will be estimated using biased estimator
#' and supposed error with no square roots in dummy observations will be reproduced
#' FALSE by default
#' @return dummy list containing:
#' X_cniw,  Y_cniw
#' X_sc, Y_sc
#' X_io, Y_io
#' X_plus, Y_plus binding all corresponging Xs and Ys
#' @export
#' @examples 
#' data(Yraw)
#' dummy <- bvar_conj_lambda2dummy(Yraw, p = 4, lambda = c(0.2,1,1,1,100,100))
bvar_conj_lambda2dummy <- function(Y_in, Z_in=NULL, constant=TRUE, p=4, 
                                   lambda=c(0.2,1,1,1,100,100), 
                                   delta=1,
                                   s2_lag=NULL, 
                                   y_bar_type = c("initial", "all"),
                                   carriero_hack = FALSE) {
  
  l_1 <- lambda[1]
  # l_kron <- 1 # for the reader of this code. l_kron is always 1 for conjugate NIW
  l_lag <- lambda[2]
  l_sc <- lambda[3]
  l_io <- lambda[4]
  l_const <- lambda[5]
  l_exo <- lambda[6]
  
  y_bar_type <- match.arg(y_bar_type)
  
  Y_in <- as.matrix(Y_in) # to clear tbl_df if present :)
  m <- ncol(Y_in)
  
  # get variable names
  endo_varnames <- bvar_get_endo_varnames(Y_in)
  exo_varnames <- bvar_get_exo_varnames(Z_in, constant=constant)
  
  # calculate d, the number of exogeneous regressors
  if (is.null(Z_in)) {
    d <- 1*constant
  } else {
    d <- ncol(Z_in) + 1*constant
  }
  
  # if requested add constant to exogeneous regressors
  # constant to the left of other exo variables
  if (constant) Z <- cbind(rep(1, nrow(Y_in)), Z_in)
  
  k <- m*p+d
  
  ##### create vector of delta from Y_in and delta descriptions
  delta <- bvar_conj_delta(Y_in, delta)
  
  
  ######  estimate sigma^2 from univariate AR(p) processes
  # Carriero matlab code: always AR(1)! 
  if (is.null(s2_lag)) s2_lag <- p
  sigmas_sq <- bvar_conj_sigma2(Y_in, s2_lag, carriero_hack = carriero_hack)
    

  
  ###### calculate y_bar 
  if (y_bar_type=="initial") sc_io_numrows <- p
  if (y_bar_type=="all") sc_io_numrows <- nrow(Y_in)
  
  
  y_bar <- apply(as.matrix(Y_in[1:sc_io_numrows,],nrow=sc_io_numrows), 2, mean) # vector [m x 1] of mean values of each endo-series
  if (is.null(Z)) { # special case of no constant and no exo vars
    z_bar <- NULL 
  }
  if (!is.null(Z)) { # "as.matrix" is needed to avoid errors for sc_io_numrows=1 or d=1
    z_bar <- 
       apply(as.matrix(Z[1:sc_io_numrows,],nrow=sc_io_numrows), 2, mean) # vector [d x 1] of mean values of each exo-series
  }
  
  
  # sc, sum of coefficients prior
  Y_sc <- NULL
  X_sc <- NULL
  if (!is.na(l_sc)) {
    Y_sc <- matrix(0, m, m) # zero matrix [m x m]
    colnames(Y_sc) <- endo_varnames
    
    diag(Y_sc) <- delta * y_bar / l_sc
    
    exo_dummy <- matrix(0, m, d)
    
    X_sc <- cbind( kronecker(matrix(1, 1, p), Y_sc) , exo_dummy)  # zero matrix [m x k]
    colnames(X_sc) <- bvar_create_X_colnames(endo_varnames, exo_varnames, p)
  }
  
  # io, dummy initial observation
  Y_io <- NULL
  X_io <- NULL
  if (!is.na(l_io)) {
    Y_io <- matrix(delta * y_bar/l_io, nrow=1)
    colnames(Y_io) <- endo_varnames
    
    X_io <- matrix(c(rep(delta * y_bar/l_io, p), z_bar/l_io), nrow=1)
    colnames(X_io) <- bvar_create_X_colnames(endo_varnames, exo_varnames, p)
  }
  
  # dummy cNIW = conjugate Normal Inverse Wishart
  y_cniw_block_1 <- diag(sqrt(sigmas_sq)*delta)/l_1
  y_cniw_block_2 <- matrix(0, nrow=m*(p-1), ncol=m)
  y_cniw_block_3 <- diag(sqrt(sigmas_sq))
  y_cniw_block_4 <- matrix(0, nrow=1, ncol=m)
  
  if (carriero_hack) { # wrong formulas, but useful for comparison with carriero code
    y_cniw_block_1 <- diag(sigmas_sq*delta)/l_1
    y_cniw_block_3 <- diag(sigmas_sq)
  }
  
  
  Y_cniw <- rbind(y_cniw_block_1, y_cniw_block_2, y_cniw_block_3, y_cniw_block_4) 
  colnames(Y_cniw) <- endo_varnames
  
  
  x_cniw_block_1 <- cbind( kronecker(diag((1:p)^l_lag), diag(sqrt(sigmas_sq)) )/l_1, matrix(0, nrow=m*p, ncol=d))
  x_cniw_block_2 <- matrix(0, nrow=m, ncol=k)
  x_cniw_block_3 <- c( rep(0, m*p), rep(1/l_const, constant), rep(1/l_exo, d-constant) )
  
  if (carriero_hack) { # wrong formulas, but useful for comparison with carriero code
    x_cniw_block_1 <- cbind( kronecker(diag((1:p)^l_lag), diag(sigmas_sq) )/l_1, matrix(0, nrow=m*p, ncol=d))
  }
  
  
  X_cniw <- rbind(x_cniw_block_1, x_cniw_block_2, x_cniw_block_3)
  rownames(X_cniw) <- NULL
  colnames(X_cniw) <- bvar_create_X_colnames(endo_varnames, exo_varnames, p)
  
  
  
  X_plus <- rbind(X_cniw, X_sc, X_io)
  Y_plus <- rbind(Y_cniw, Y_sc, Y_io)
  
  dummy <- list(X_sc=X_sc, Y_sc=Y_sc,
                X_io=X_io, Y_io=Y_io,
                X_cniw=X_cniw, Y_cniw=Y_cniw,
                X_plus=X_plus, Y_plus=Y_plus,
                sigma_sq=sigmas_sq)
  return(dummy)
}

#' Get endogeneous variable names from supplied data
#' 
#' Get endogeneous variable names from supplied data
#' 
#' Get endogeneous variable names from supplied data
#' 
#' @param Y_in [T_in x m] multivariate time series
#' @export
#' @return vector of variable names
#' @examples 
#' data(Yraw)
#' bvar_get_endo_varnames(Yraw)
bvar_get_endo_varnames <- function(Y_in) {
  endo_varnames <- colnames(Y_in)
  if (is.null(endo_varnames)) endo_varnames <- paste0("endo_",1:ncol(Y_in))
  return(endo_varnames)
}

#' Get exogeneous variable names from supplied data
#' 
#' Get exogeneous variable names from supplied data
#' 
#' Get exogeneous variable names from supplied data
#' 
#' @param constant logical, default is TRUE, whether the constant should be included
#' @param Z_in [T_in x m] multivariate time series
#' @export
#' @return vector of variable names
#' @examples 
#' data(Yraw)
#' bvar_get_exo_varnames(Yraw)
bvar_get_exo_varnames <- function(Z_in, constant=TRUE) {
  exo_varnames <- colnames(Z_in)
  if (is.null(exo_varnames) & !is.null(Z_in)) exo_varnames <- paste0("exo_",1:ncol(Z_in))
  
  exo_varnames <- c(rep("const",constant),exo_varnames)
  return(exo_varnames)
}

#' Get X column names from exo_varnames and endo_varnames
#' 
#' Get X column names from exo_varnames and endo_varnames
#' 
#' Get X column names from exo_varnames and endo_varnames
#' 
#' @param p number of lags
#' @param exo_varnames exogeneous variable names
#' @param endo_varnames endogeneous variable names
#' @export
#' @return vector of variable names
#' @examples 
#' data(Yraw)
#' exo_names <- bvar_get_exo_varnames(Yraw)
#' endo_names <- bvar_get_endo_varnames(Yraw)
#' bvar_create_X_colnames(endo_names, exo_names, p=2)
bvar_create_X_colnames <- function(endo_varnames, exo_varnames, p) {
  X_colnames <- paste0(rep(endo_varnames,p),", l=",rep(1:p,each=length(endo_varnames)))
  X_colnames <- c(X_colnames, exo_varnames)
  return(X_colnames)
}

#' Create model setup from lambdas 
#' 
#' Create model setup from lambdas 
#' 
#' Create model setup from lambdas 
#' Lambdas specification is based on Carriero p. 52-53
#'
#' @param p number of lags
#' @param Y_in multivariate time series
#' @param lambda vector = (l_1, l_lag, l_sc, l_io, l_const, l_exo), the l_kron is set to 1 automatically for 
#' conjugate N-IW prior. Short summary valid for NO sc/io case:
#' sd(const in eq i) = l_const * sigma_i
#' sd(exo in eq i)= l_exo * sigma_i
#' sd(coef for var j lag l in eq i) = l_1*sigma_i/sigma_j/l^l_lag
#' lambdas may be Inf
#' l_io or l_sc equal to NA means no corresponding dummy observations
#' @param Z_in exogeneous variables
#' @param constant logical, default is TRUE, whether the constant should be included
#' @param s2_lag number of lags in AR() model used to estimate s2 (equal to p by default)
#' Carriero uses 1 in his matlab code
#' @param delta vector [m x 1] or scalar or "AR1". Are used for prior Phi_1 and in sc/io dummy observations
#' Scalar value is replicated m times. If set to "AR1" then deltas will be estimated as AR(1) coefficients (but not greater than one).
#' Diagonal of Phi_1 is equal to delta. y_bar is multiplied by delta componentwise.
#' @param y_bar_type (either "all" or "initial"). Determines how y_bar for sc and io dummy is calculated.
#' "all": y_bar is mean of y for all observations, "initial": p initial observations
#' Carriero: all, Sim-Zha: initial
#' @param carriero_hack logical, if TRUE sigma^2 will be estimated using biased estimator
#' and supposed error with no square roots in dummy observations will be reproduced
#' FALSE by default
#' @param v_prior prior value of hyperparameter nu, m+2 by default, 
#' may use character formula involving m, k, p, T_dummy, constant, something like "T_dummy+2"
#' @return setup list containing
#' X, Y, X_plus, Y_plus,
#' p, v_prior
#' @export
#' @examples 
#' data(Yraw)
#' setup <- bvar_conj_setup(Yraw, p = 4, lambda = c(0.2,1,1,1,100,100))
bvar_conj_setup <- function(Y_in, Z_in=NULL, constant=TRUE, p=4, 
                            lambda=c(0.2,1,1,1,100,100), 
                            delta=1,
                            s2_lag=NULL, 
                            y_bar_type = c("initial", "all"),
                            v_prior=NULL,
                            carriero_hack = FALSE) {
  
  # general philosophy: calculate everything but use bare minimum
  
  m <- ncol(Y_in)
  
  Y <- bvar_build_Y(Y_in=Y_in, p=p)
  X <- bvar_build_X(Y_in=Y_in, Z_in=Z_in, constant=constant, p=p)
  
  dummy <- bvar_conj_lambda2dummy(Y_in=Y_in, Z_in=Z_in, constant=constant, p=p,
                                  lambda=lambda, delta=delta, s2_lag=s2_lag,
                                  y_bar_type=y_bar_type,
                                  carriero_hack = carriero_hack)
    
  if (is.null(v_prior)) v_prior <- m + 2
  
  if (is.character(v_prior)) { # parse character formula for v_prior
    k <- ncol(dummy$X_plus)
    T_dummy <- nrow(dummy)
    v_prior <- eval(parse(text=v_prior))
  }
  
  # these are needed only for bvar_conj_mdd() and not for forecasting/estimation:
  prior <- bvar_conj_lambda2hyper(Y_in=Y_in, Z_in=Z_in, constant=constant, p=p,
                                   lambda=lambda, delta=delta, s2_lag=s2_lag,
                                   y_bar_type=y_bar_type,
                                   carriero_hack = carriero_hack)  
  
  setup <- list(Y = Y, X = X, Y_plus = dummy$Y_plus, X_plus = dummy$X_plus, v_prior = v_prior, p=p,
                constant = constant, 
                Omega_prior = prior$Omega, 
                S_prior = prior$S,
                Omega_prior_root = prior$Omega_root,
                Phi_prior = prior$Phi,
                S_prior = prior$S,
                Y_in = Y_in,
                Z_in = Z_in)
  return(setup)
}

#' Estimate bvar conjugate model from setup
#' 
#' Estimate bvar conjugate model from setup
#'  
#' Estimate bvar conjugate model from setup
#' 
#' @param setup list containing:
#' X [T x k]
#' Y [T x m]
#' X_plus [T_plus x k] 
#' Y_plus [T_plus x m]
#' v_prior 
#' p number of lags
#' @param keep the number of simulations
#' If keep is zero only posteriors are calculated
#' @param verbose TRUE will show some debugging messages, FALSE by default
#' @return model list containing
#' X, Y, X_plus, Y_plus,
#' p, v_prior,
#' sample (if keep>0)
#' Omega_post, v_post, S_post, Phi_post
#' @export
#' @examples 
#' data(Yraw)
#' setup <- bvar_conj_setup(Yraw, p = 4, lambda = c(0.2,1,1,1,100,100))
#' model <- bvar_conj_estimate(setup, keep=10)
bvar_conj_estimate <- function(setup, keep=2000, verbose=FALSE) {
  
  T <- nrow(setup$Y) # number of observations
  
  # calculate hyperpar of posterior
  X_star <- rbind(setup$X_plus, setup$X)
  Y_star <- rbind(setup$Y_plus, setup$Y)
  
  if (verbose) message("Calculating posterior hyperparameters...")
  post <- bvar_conj_dummy2hyper(Y_star, X_star)
  v_post <- setup$v_prior + T
  
  # simulate
  if (verbose) message("Simulating...")
  if (keep>0) setup$sample <- 
    bvar_conj_simulate(v_post = v_post, 
                       Omega_post_root = post$Omega_root, 
                       S_post = post$S, Phi_post = post$Phi, 
                       verbose = verbose, keep = keep) 
  
  setup$v_post <- v_post
  setup$S_post <- post$S
  setup$Omega_post <- post$Omega
  setup$Phi_post <- post$Phi
                                       
  return(setup)
}


#' Set conjugate N-IW priors from lambdas as in Carriero
#' 
#' Set conjugate N-IW priors from lambdas as in Carriero
#' 
#' Set conjugate N-IW priors from lambdas as in Carriero
#' Based on Carriero p. 52-53
#'
#' @param p number of lags
#' @param Y_in [T_in x m] multivariate time series
#' @param lambdas vector = (l_1, l_lag, l_sc, l_io, l_const, l_exo), the l_kron is set to 1 automatically for 
#' conjugate N-IW prior. Short summary:
#' sd(const in eq i) = l_const * sigma_i
#' sd(exo in eq i)= l_exo * sigma_i
#' sd(coef for var j lag l in eq i) = l_1*sigma_i/sigma_j/l^l_lag
#' lambdas may be Inf
#' l_io or l_sc equal to NA means no corresponding dummy observations
#' @param Z_in [T_in x ...] exogeneous variables
#' @param constant logical, default is TRUE, whether the constant should be included
#' @param s2_lag number of lags in AR() model used to estimate s2 (equal to p by default)
#' Carriero uses 1 in his matlab code
#' @param delta vector [m x 1] or scalar or "AR1". Are used for prior Phi_1 and in sc/io dummy observations
#' Scalar value is replicated m times. If set to "AR1" then deltas will be estimated as AR(1) coefficients (but not greater than one).
#' Diagonal of Phi_1 is equal to delta. y_0_bar is multiplied by delta componentwise.
#' @param y_0_bar_type (either "all" or "initial"). Determines how y_bar0 for sc and io dummy is calculated.
#' "all": y_bar0 is mean of y for all observations, "initial": p initial observations
#' Carriero: all, Sim-Zha: initial
#' @return priors list containing Phi_prior [k x m], Omega_prior [k x k], S_prior [m x m], v_prior [1x1],
#' where k = mp+d
#' @export
#' @examples 
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4, lambdas = c(0.2,1,1,1,100,100))
#' model <- bvar_conjugate0(priors = priors, keep=100)
#' # if something goes wrong then we need info!
#' model <- bvar_conjugate0(priors = priors, keep=10, verbose=TRUE)
Carriero_priors <- function(Y_in, Z_in=NULL, constant=TRUE, p=4, 
                            lambdas=c(0.2,1,1,1,100,100), 
                            delta=1,
                            s2_lag=NULL, 
                            y_0_bar_type = c("initial", "all") ) {
  l_1 <- lambdas[1]
  l_kron <- 1
  l_lag <- lambdas[2]
  l_sc <- lambdas[3]
  l_io <- lambdas[4]
  l_const <- lambdas[5]
  l_exo <- lambdas[6]
  
  y_0_bar_type <- match.arg(y_0_bar_type)
  
  Y_in <- as.matrix(Y_in) # to clear tbl_df if present :)
  
  m <- ncol(Y_in)

  
  # get variable names
  endo_varnames <- colnames(Y_in)
  if (is.null(endo_varnames)) endo_varnames <- paste0("endo_",1:m)
  

  exo_varnames <- NULL  
  
  # calculate d, the number of exogeneous regressors
  if (is.null(Z_in)) {
    d <- 1*constant
    exo_varnames <- "const"
  } else {
    d <- ncol(Z_in) + 1*constant

    exo_varnames <- colnames(Z_in)
    if (is.null(exo_varnames)) exo_varnames <- paste0("exo_",1:ncol(Z_in))
    if (constant) exo_varnames <- c("const",exo_varnames)
  }
  
  # if requested add constant to exogeneous regressors
  # constant to the left of other exo variables
  if (constant) Z <- cbind(rep(1, nrow(Y_in)), Z_in)
  
  k <- m*p+d

  
  if (delta=="AR1") { # set deltas as AR(1) coefficients but no more than 1
    delta <- rep(1,m) # reserve space
    for (j in 1:m) {
      y_uni <- Y_in[,j] # univariate time series
      # more robust version: fails only in the case of  severe multicollinearity
      AR_1 <- ar.ols(y_uni, aic=FALSE, order.max = 1) # AR(1) model
      delta[j] <- AR_1$ar
      if (delta[j]>1) delta[j] <- 1
    }
    
  }
  if (length(delta)==1) delta <- rep(delta, m) # copy scalar values
  if (! length(delta)==m ) stop("Length of delta should be equal to 1 or m")
  
  #VAR_in <- match.arg(VAR_in)
  
  if (!l_kron==1) warning("Conjugate N-IW is impossible for lambda_2 <> 1")
  
  # Litterman takes 6 lags in AR(p)
  
  # estimate sigma^2 from univariate AR(p) processes
  # Carriero matlab code: always AR(1)! 
  if (is.null(s2_lag)) s2_lag <- p
  
  sigmas_sq <- rep(NA, m)
  for (j in 1:m) {
    
    y_uni <- Y_in[,j] # univariate time series

    
    # old version: it fails when ML estimation fails :)
    #AR_p <- forecast::Arima(y_uni, order = c(p,0,0), method="ML") # AR(p) model
    #sigmas_sq[j] <- AR_p$sigma2
    
    
   
    
    # more robust version: fails only in the case of  severe multicollinearity
    AR_p <- ar.ols(y_uni, aic=FALSE, order.max = s2_lag) # AR(p) model
    resid <- tail(AR_p$resid,-s2_lag) # omit first p NA in residuals
    sigmas_sq[j] <- sum(resid^2)/(length(resid)-s2_lag-1)
  }
  
  # set Phi_prior
  Phi_1 <- matrix(0, m,m)
  diag(Phi_1) <- delta
  Phi_prior <- t( cbind(Phi_1, matrix(0, nrow=m, ncol=k-m)) )
  
  colnames(Phi_prior) <- paste0("eq_",endo_varnames)
  rownames(Phi_prior) <- c(paste0(rep(endo_varnames, p),", l=",rep(1:p,each=m)) , exo_varnames)
  
  S_prior <- diag(sigmas_sq)
  v_prior <- m+2
  
  colnames(S_prior) <- endo_varnames
  rownames(S_prior) <- endo_varnames
  
  # set Omega_prior
  # the diagonal of Omega_prior begins with endogeneous part:
  endo_diagonal <- l_1^2/rep(sigmas_sq, p)/rep((1:p)^(2*l_lag), each=m)
  
  # and ends with exogeneous part:
  exo_diagonal <- c( rep(l_const^2,constant), rep(l_exo^2,d-constant) )
  
  Omega_diagonal <- c(endo_diagonal, exo_diagonal)
  # and set zero prior covariances
  Omega_prior <- diag(Omega_diagonal)
  
  
  
  
  # create dummy observations
  if (y_0_bar_type=="initial") sc_io_numrows <- p
  if (y_0_bar_type=="all") sc_io_numrows <- nrow(Y_in)
  
  
  y_0_bar <- apply(as.matrix(Y_in[1:sc_io_numrows,],nrow=sc_io_numrows), 2, mean) # vector [m x 1] of mean values of each endo-series
  if (is.null(Z)) z_bar <- NULL # special case of no constant and no exo vars
  if (!is.null(Z)) z_bar <- 
    apply(as.matrix(Z[1:sc_io_numrows,],nrow=sc_io_numrows), 2, mean) # vector [d x 1] of mean values of each exo-series
  # "as.matrix" above is needed to avoid errors for sc_io_numrows=1 or d=1
  
  
  
  # sum of coefficients prior
  Y_dummy_sc <- NULL
  X_dummy_sc <- NULL
  if (!is.na(l_sc)) {
    Y_dummy_sc <- matrix(0, m, m) # zero matrix [m x m]
    diag(Y_dummy_sc) <- delta * y_0_bar / l_sc
  
    exo_dummy <- matrix(0, m, d)
    
    X_dummy_sc <- cbind( kronecker(matrix(1, 1, p), Y_dummy_sc) , exo_dummy)  # zero matrix [m x k]
  }
  
  # dummy initial observation
  Y_dummy_io <- NULL
  X_dummy_io <- NULL
  if (!is.na(l_io)) {
    Y_dummy_io <- matrix(delta * y_0_bar/l_io, nrow=1)
    X_dummy_io <- matrix(c(rep(delta * y_0_bar/l_io, p), z_bar/l_io), nrow=1)
  }
  
  # dummy cNIW = conjugate Normal Inverse Wishart
  y_cniw_block_1 <- diag(sqrt(sigmas_sq)*delta)/l_1
  y_cniw_block_2 <- matrix(0, nrow=m*(p-1), ncol=m)
  y_cniw_block_3 <- diag(sqrt(sigmas_sq))
  y_cniw_block_4 <- matrix(0, nrow=1, ncol=m)
  Y_dummy_cniw <- rbind(y_cniw_block_1, y_cniw_block_2, y_cniw_block_3, y_cniw_block_4) 
  
  x_cniw_block_1 <- cbind( kronecker(diag((1:p)^l_lag), diag(sqrt(sigmas_sq)) )/l_1, matrix(0, nrow=m*p, ncol=d))
  x_cniw_block_2 <- matrix(0, nrow=m, ncol=k)
  x_cniw_block_3 <- c( rep(0, m*p), rep(1/l_const, constant), rep(1/l_exo, d-constant) )
  X_dummy_cniw <- rbind(x_cniw_block_1, x_cniw_block_2, x_cniw_block_3)
    
    
    
  
  # robust calculation of (Omega_prior)^{-1/2}
  # the diagonal of Omega_prior begins with endogeneous part:
  endo_diagonal_m05 <- 1/l_1*rep(sqrt(sigmas_sq), p)/rep((1:p)^(l_lag), each=m)
  
  # and ends with exogeneous part:
  exo_diagonal_m05 <- rep(1/l_exo,d)
  if (constant) exo_diagonal_m05[1] <- 1/l_const
  
  Omega_diagonal_m05 <- c(endo_diagonal_m05, exo_diagonal_m05)
  # and set zero prior covariances
  Omega_prior_m05 <- diag(Omega_diagonal_m05)
  # end of robust calculation of Omega_diagonal_m05
  
  
  # order of dummies does not matter
  X_dummy <- rbind(X_dummy_io, X_dummy_sc)
  Y_dummy <- rbind(Y_dummy_io, Y_dummy_sc)
  
    
  priors <- list(v_prior=v_prior, S_prior=S_prior, 
                 Phi_prior=Phi_prior, Omega_prior=Omega_prior, 
                 Y_dummy=Y_dummy, X_dummy=X_dummy,
                 X_dummy_io=X_dummy_io, 
                 Y_dummy_io=Y_dummy_io,
                 X_dummy_sc=X_dummy_sc,
                 Y_dummy_sc=Y_dummy_sc,
                 X_dummy_cniw=X_dummy_cniw,
                 Y_dummy_cniw=Y_dummy_cniw,
                 sc_io_numrows=sc_io_numrows,
                 Y_in=Y_in, Z_in=Z_in, p=p, # to avoid duplicating
                 sigmas_sq = sigmas_sq,
                 endo_varnames=endo_varnames,
                 exo_varnames=exo_varnames,
                 Omega_prior_m05=Omega_prior_m05,
                 lambdas=lambdas) # get more info from function
  
  return(priors)
}





#' Set conjugate N-IW priors as in matlab code of Koops-Korobilis
#' 
#' Set conjugate N-IW priors as in matlab code of Koops-Korobilis
#' 
#' Set conjugate N-IW priors as in matlab code of Koops-Korobilis. Mainly for testing.
#'
#' @param p number of lags
#' @param Y_in multivariate time series
#' @param constant logical, default is TRUE, whether the constant should be included
#' @param Z_in exogeneous variables
#' @return priors list containing Phi_prior [k x m], Omega_prior [k x k], S_prior [m x m], v_prior [1x1],
#' where k = mp+d
#' @export
#' @examples 
#' data(Yraw)
#' priors <- KK_code_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors, keep=100)
KK_code_priors <- function(Y_in, Z_in=NULL, constant=TRUE, p=4) {

  Y_in <- as.matrix(Y_in) # to clear tbl_df if present :)
  
  m <- ncol(Y_in)

  
  # get variable names
  endo_varnames <- colnames(Y_in)
  if (is.null(endo_varnames)) endo_varnames <- paste0("endo_",1:m)
  
  
  exo_varnames <- NULL  
  
  # calculate d, the number of exogeneous regressors
  if (is.null(Z_in)) {
    d <- 1*constant
    exo_varnames <- "const"
  } else {
    d <- ncol(Z_in) + 1*constant
    
    exo_varnames <- colnames(Z_in)
    if (is.null(exo_varnames)) exo_varnames <- paste0("exo_",1:ncol(Z_in))
    if (constant) exo_varnames <- c("const",exo_varnames)
  }
  
  # if requested add constant to exogeneous regressors
  # constant to the left of other exo variables
  if (constant) Z <- cbind(rep(1, nrow(Y_in)), Z_in)
  
  k <- m*p+d
  
  v_prior <- m + 1
  S_prior <- diag(m)
  Omega_prior <- 10*diag(k) 
  Phi_prior <- matrix(0, nrow=k, ncol=m)
    

  X_dummy <- NULL
  Y_dummy <- NULL
  
  
  priors <- list(v_prior=v_prior, S_prior=S_prior, 
                 Phi_prior=Phi_prior, 
                 Omega_prior=Omega_prior, 
                 Y_dummy=Y_dummy, X_dummy=X_dummy,
                 Y_in=Y_in, Z=Z_in, p=p,
                 endo_varnames=endo_varnames,
                 exo_varnames=exo_varnames) # get more info from function
  
  return(priors)
}


#' Simulate from conjugate N-IW posterior distribution
#' 
#' Simulate from conjugate N-IW posterior distribution
#' 
#' Simulate from conjugate N-IW posterior distribution
#' 
#' @param keep number of simulations
#' @param verbose TRUE will show some debugging messages, FALSE by default
#' @param Omega_post_root [k x k] "root" of Omega posterior matrix
#' @param v_post posterior nu, scalar
#' @param S_post [m x m] posterior scale matrix for noise covariance
#' @param Phi_post [k x m] posterior expected values for coefs
#' @export
#' @return coda object with simulations
#' @examples 
#' data(Yraw)
#' dummy <- bvar_conj_lambda2dummy(Yraw,p=2)
#' hyper <- bvar_conj_dummy2hyper(dummy$Y_cniw, dummy$X_cniw)
#' # these are priors but just for testing let's pretend that they are posteriors
#' post_sample <- bvar_conj_simulate(v_post=10,
#' hyper$Omega_root, hyper$S, hyper$Phi)
bvar_conj_simulate <- function(v_post, Omega_post_root, S_post, Phi_post, 
                               verbose=FALSE, keep=10) {
  k <- nrow(Phi_post)
  m <- ncol(Phi_post)
  
  answer <- matrix(0, nrow=keep, ncol = m*k + m*m)
  
  for (i in 1:keep) {
    if ((verbose) & (i %% 10^3 == 0)) message("Iteration ",i," out of ",keep)
    
    Sigma <- MCMCpack::riwish(v_post,S_post) 
    
    
    # fast way to generate matrix-variate normal
    V <- matrix(rnorm(k*m), ncol = m) # [(mp+d) x m] matrix of standard normal
    Phi <- Phi_post + Omega_post_root %*% V %*% chol(Sigma)
    
    Phi_vec <- as.vector(Phi)
    Sigma_vec <- as.vector(Sigma) # length = m x m
    
    answer[i,] <- c(Phi_vec, Sigma_vec)
  } # cycle from 1 to keep

  # save as mcmc object to make some good functions available
  answer <- coda::as.mcmc(answer)
  return(answer)
}


#' Compute inverse of symmetric positive definite matrix using Cholesky decomposition
#'
#' Compute inverse of symmetric positive definite matrix using Cholesky decomposition
#'  
#' Compute inverse of symmetric positive definite matrix using Cholesky decomposition
#' 
#' @param A symmetric matrix 
#' @return inverse of A
#' @export
#' @examples
#' A <- matrix(c(2,1,1,2),nrow=2)
#' sym_inv(A)
sym_inv <- function(A) {
  A_chol <- try(chol(A), silent=TRUE)
  if (class(A_chol)=="try-error") {
    inv_A <- MASS::ginv(A) # for more stable results in multicollinearity cases
    attr(inv_A,"Moore-Penrose") <- TRUE
  } else {
  inv_A <- chol2inv(A_chol)
  }
  return(inv_A)
}


#' check whether matrix is diagonal
#'
#' check whether matrix is diagonal
#'  
#' check whether matrix is diagonal
#' 
#' @param A symmetric matrix 
#' @param epsilon tolerance, default is 0. 
#' @return logical, TRUE/FALSE, TRUE for diagonal matrices
#' @export
#' @examples
#' A <- matrix(c(2,1,1,2),nrow=2)
#' sym_inv(A)
is.diagonal <- function(A, epsilon=0) {
  non_diag_elements <- !diag(nrow(A))
  answer <- all( abs(A[non_diag_elements]) <= epsilon )
  return(answer)
}

#' predict with conjugate Normal-Inverse-Wishart bayesian VAR model
#'
#' predict with conjugate Normal-Inverse-Wishart bayesian VAR model
#'  
#' predict with conjugate Normal-Inverse-Wishart bayesian VAR model
#' 
#' @param model estimated conjugate N-IW model
#' @param h number of periods for forecasting
#' @param level confidence levels for prediction intervals
#' @param Y_in (NULL by default) past values of endogeneous variables (shold have at least p observations).
#' If NULL, then Y_in supplied for estimation will be used. For out-of-sample forecast only last p values of Y_in are used
#' @param Z_f future values of exogeneous variables
#' @param type ("prediction" by default) type of interval: "prediction" incorporates uncertainty about
#' future shocks; "credible" deals only with parameter uncertainty.
#' @param output (default "long") --- long or wide table for mean/quantiles of forecasts
#' @param out_of_sample logical, default is TRUE, whether forecasts are out-of-sample or in-sample.
#' If forecasts are not out-of-sample, then parameter h is ignored
#' @param include (default is c("mean", "median", "sd", "raw")) what type of summary to provide
#' If "raw" is present raw forecasts will be reported. If only "raw" is present
#' then function will return coda mcmc object with raw forecasts. 
#' Otherwise raw forecasts will be attached as attribute.
#' @param fast_forecast logical, FALSE by default. If TRUE then only mean forecast is calculated,
#' posterior expected values of hyperparameters are used. No confidence intervals, no sd, no median. 
#' This mode is activated by default if there are no simulations in supplied model.
#' @param verbose (default FALSE) if true some messages will be printed
#' @export
#' @return forecast results
#' @examples 
#' data(Yraw)
#' setup <- bvar_conj_setup(Yraw, p = 4)
#' model <- bvar_conj_estimate(setup, keep=100)
#' bvar_conj_forecast(model, h=2, output="wide")
#' bvar_conj_forecast(model, out_of_sample = FALSE, include="mean")
bvar_conj_forecast <- function(model, 
                               Y_in=NULL, Z_f=NULL,
                               output=c("long","wide"),
                               h=1, out_of_sample=TRUE,
                               type=c("prediction","credible"), level=c(80,95),
                               include=c("mean","median","sd","interval", "raw"),
                               fast_forecast=FALSE,
                               verbose=FALSE) {
  
  
  # select type of prediction specified by user
  type <- match.arg(type)
  output <- match.arg(output)
  
  # simplify notation, extract params from attribute
  T <- nrow(model$Y) # number of observations in regression
  p <- model$p
  k <- ncol(model$X)
  m <- ncol(model$Y)
  d <- k - m*p
  T_dummy <- nrow(model$X_plus)
  
  keep <- 0
  if ("sample" %in% names(model)) keep <- nrow(model$sample) 
  
  constant <- model$constant 
  
  if ((keep==0) & (!fast_forecast)) {
    if (verbose) message("No simulations found in 'model', fast_forecast option is set to TRUE.")
    fast_forecast <- TRUE
  }
  
  if (fast_forecast) {
    level <- NULL
    include <- "mean"
    keep <- 1
    
    # type of confidence interval does not make sense for fast mean forecast
    # but setting type to 'credible' avoids simulations of future epsilon, so
    type <- "credible"
  }
  
  # if Y_in is not supplied take Y_in from estimation. Or is it better to save it in model list?
  if (is.null(Y_in)) Y_in <- bvar_get_Y_in(model$Y, model$X, p=p)
  
  
  # in case of in-sample forecast h is set to T, or is is better set to NA?
  if (!out_of_sample) h <- T
  
  
  # sanity check
  if (!is.null(Z_f))
    if (!nrow(Z_f)==h) stop("I need exactly h = ",h," observations for exogeneous variables.")
  
  if (nrow(Y_in)<p) stop("Model has ",p," lags. To predict I need at least ",p," observations, but only ",nrow(Y_in)," are provided.")
  
  
  # if requested add constant to exogeneous regressors
  if (constant) Z_f <- cbind(rep(1, h), Z_f)
  
  
  if (out_of_sample) {
    # take last p observations of Y_in (out-of-sample forecast)
    Y_in <- tail(Y_in, p)
  } else {
    # take first p observations of Y_in (in-sample forecast)
    Y_in <- head(Y_in, p)
  }
  
  e_t <- rep(0, m) # ok for bayesian credible intervals  
  
  # space to store all forecasted values
  forecast_raw <- matrix(0, nrow = keep, ncol = m*h)
  
  
  
  x_t <- rep(0, k)
  
  for (i in 1:keep) {
    # forecast h steps for given sampling of Phi
    
    if (fast_forecast) { # use posterior expected value of Phi
      Phi <- model$Phi_post
    } else { # use Phi^[s] and Sigma^[s]
      Phi <- matrix(model$sample[i,1:(k*m)], nrow=k)
      
      
      Sigma <- matrix(model$sample[i,(k*m+1):(m*k + m*m)],nrow=m) # Sigma [m x m]
      # find square root of draw from Sigma (code is part of mvtnorm function)
      ev <- eigen( Sigma, symmetric = TRUE)
      if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
        warning("Omega_post is numerically not positive definite")
      }
      # precalculate R to do less operations in case h>1
      R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    }
    
    Phi_transp <- t(Phi) # precalculate to do less operations in case h>1
    
    for (j in 1:h) {
      
      
      if (out_of_sample) { # out-of-sample forecast
        # fill exogeneous values in x_t
        x_t[(m * p + 1):(m * p + d)] <- Z_f[j,]
        
        # fill endogeneous values in x_t (first out-of-sample forecast)
        if (j == 1)
          x_t[1:(m * p)] <- as.vector(t(Y_in)[,p:1])
        
        # fill endogeneous values recursively (second+ out-of-sample forecast)
        if (j > 1) {
          if (p>1) x_t[(m + 1):(m * p)] <- x_t[1:(m * (p - 1))]
          x_t[1:m] <- y_t
        }
      } else { # in-sample forecast
        x_t <- model$X[j,]
      }
      
      if (type=="prediction") e_t <- R %*% rnorm(m) 
      # e_t is 0 for bayesian credible intervals
      
      y_t <- Phi_transp %*% x_t + e_t
      forecast_raw[i, (m*(j-1)+1):(m*j)] <- y_t
    } # end_for j in 1:h
  } # end_for i in 1:keep
  
  # save as mcmc object for standartisation
  forecast_raw <- coda::as.mcmc(forecast_raw)
  
  # we have m endogeneous variables and h forecasts for each
  id_block <- data.frame(variable=rep( bvar_get_endo_varnames(Y_in), h), h=rep(1:h, each=m))
  
  forecast_summary <- NULL
  
  # calculate mean
  if ("mean" %in% include) {
    what <- rep("mean", h*m)
    value <- apply(forecast_raw, 2, mean)
    block <- cbind(id_block, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }
  
  # calculate median
  if ("median" %in% include) {
    what <- rep("median", h*m)
    value <- apply(forecast_raw, 2, median)
    block <- cbind(id_block, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }
  
  # sd
  if ("sd" %in% include) {
    what <- rep("sd", h*m)
    value <- apply(forecast_raw, 2, sd)
    block <- cbind(id_block, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }
  
  # calculate quantiles
  if ("interval" %in% include) {
    for (lev in level) {
      # lower
      what <- rep(paste0("lower_",lev), h*m)
      value <- apply(forecast_raw, 2, function(x) quantile(x, probs=(1-lev/100)/2))
      block <- cbind(id_block, what, value) # block of information
      forecast_summary <- rbind(forecast_summary, block)
      
      # upper
      what <- rep(paste0("upper_",lev), h*m)
      value <- apply(forecast_raw, 2, function(x) quantile(x, probs=(1+lev/100)/2))
      block <- cbind(id_block, what, value) # block of information
      forecast_summary <- rbind(forecast_summary, block)
    }
  }
  
  rownames(forecast_summary)  <- NULL
  
  
  if ((output=="wide") & (!is.null(forecast_summary))) { # transform to wide format if requested
    forecast_summary <- reshape2::dcast(forecast_summary, variable+h~what)
  }
  
  
  if ("raw" %in% include) {
    if (length(include)==1) { # only raw forecasts are requested
      forecast_summary <- forecast_raw
    } else { # attach raw forecasts as attribute
      attr(forecast_summary, "forecast_raw") <- forecast_raw
    }
  }
  return(forecast_summary)
}





#' Estimate conjugate Normal-Inverse-Wishart bayesian VAR model
#'
#' Estimate conjugate  Normal-Inverse-Wishart bayesian VAR model
#'  
#' Estimate conjugate  Normal-Inverse-Wishart bayesian VAR model
#' 
#' @param Y_in the matrix or data.frame with endogeneous VAR variables
#' @param Z_in (NULL by default) the matrix or data.frame with exogeneous VAR variables
#' @param constant (TRUE by default) whether we should include constant
#' @param p (2 by default) the number of lags
#' @param keep (10000 by default) the number of Gibbs sampling replications to keep
#' Is ignored when the fast_forecast is TRUE.
#' @param verbose (FALSE by default)
#' @param priors the list containing at least Phi_prior [k x m], Omega_prior [k x k], 
#' S_prior [m x m], v_prior [1x1],
#' it may also contain
#' Y_dummy [T_dummy x m], X_dummy [T_dummy x k]
#' where k = mp+d
#' @param fast_forecast logical, FALSE by default. If TRUE then no simulations are done,
#' only posterior hyperparameters are calculated.
#' @param way_omega_post_root the way for (Omega_post)^{1/2} calculation: "svd" or "cholesky"
#' @return the list containing all results of bayesian VAR estimation
#' @export
#' @examples
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors, keep=100)
bvar_conjugate0 <-
  function(Y_in=NULL, Z_in=NULL, constant=TRUE, p=NULL, keep=10000, verbose=FALSE,
           priors=list(), fast_forecast=FALSE,
           way_omega_post_root = c("cholesky","svd") ) {

    exo_varnames <- priors$exo_varnames
    endo_varnames <- priors$endo_varnames
    
    way_omega_post_root <- match.arg(way_omega_post_root)
    
    if ( (is.null(Y_in)) & (!is.null(priors$Y_in)) ){
      Y_in <- priors$Y_in
      if (verbose) message("Y_in is inferred from priors data.")
    }
    
    if ( (is.null(Z_in)) & (!is.null(priors$Z_in)) ){
      Z_in <- priors$Z_in
      if (verbose) message("Z_in is inferred from priors data.")
    }

    if ( (is.null(p)) & (!is.null(priors$p)) ) {
      p <- priors$p
      if (verbose) message("Number of lags is inferred from priors data: p = ",p)
    }
    
    if ( (is.null(p)) & (is.null(priors$p)) ) {
      p <- 4
      if (verbose) message("Number of lags, p, is not specified inside and outside priors, set to p = ",p)
    }
    
    
    
    
    # if Z_in is provided it should have the same number of rows that Y_in
    if (!is.null(Z_in)) 
      if (!nrow(Y_in)==nrow(Z_in))
        stop("Number of rows in Y_in and Z_in should be equal. 
           The first p rows of Z_in are not used and may be filled with NA")
    
    
    # number of observations supplied
    T_in <- nrow(Y_in)
    
    # if requested add constant to exogeneous regressors
    if (constant) Z_in <- cbind(rep(1, nrow(Y_in)), Z_in)
    
    # transform Y_in from data.frame to matrix
    # so that cbind(NULL, Y_block) will work 
    Y_in <- as.matrix(Y_in)
    
    
    # p first observations of Y_in are used only as regressors
    Y <- tail(Y_in, -p)
    # p first observations of Z_in are not used only at all
    Z <- tail(Z_in, -p)
    
    
    # create X matrix
    X <- NULL
    for (j in 1:p) {
      Y_block <- Y_in[(p+1-j):(T_in-j),]
      X <- cbind(X, Y_block) 
    }
    X <- cbind(X,Z)
    
    # save X and Y to include them in output later:
    X_wo_dummy <- X
    Y_wo_dummy <- Y
    
    
    # dimension of deterministic regressors
    d <- ncol(Z)
    
    # here we add dummy observations
    T_dummy <- 0 
    if (!is.null(priors$Y_dummy)) {
      if (!nrow(priors$Y_dummy)==nrow(priors$X_dummy)) stop("X_dummy and Y_dummy should have the same number of rows")
      T_dummy <- nrow(priors$Y_dummy)
    }
    
    # order of dummy observation does not matter :)
    Y <- rbind(priors$Y_dummy, Y)
    X <- rbind(priors$X_dummy, X)
        
    # get dimensions, here T counts: T_in (supplied obs) - n_lag + T_dummy
    T <- nrow(Y)
    m <- ncol(Y)
    k <- m*p + d
    
    if (verbose) { 
      message("Number of lags, p = ", p)
      message("Number of endogeneos variables, m = ",m)
      message("Number of exogeneos variables (including constant), d = ",d)
      message("Number of parameters, k = mp + d = ",k)
      message("Initial number of observations, T_in = ",T_in)
      message("Number of dummy observations, T_dummy = ", T_dummy )
      message("Number of observations available for regression, T = T_in + T_dummy - p = ",T)
      message("'fast_forecast' is ", fast_forecast)
    }

    # set some bad priors for lazy guys if not supplied 
    if (is.null(priors$v_prior)) {
      priors$v_prior <- m + 1
      if (verbose) message("v_prior was not specified, set to (m+1)")
    }
    if (is.null(priors$S_prior)) {
      priors$S_prior <- diag(m)
      if (verbose) message("S_prior was not specified, set to I [m x m]")
    }
    if (is.null(priors$Omega_prior)) {
      priors$Omega_prior <- 10*diag(k) 
      if (verbose) message("Omega_prior was not specified, set to 10I [k x k]")
    }  
    if (is.null(priors$Phi_prior)) {
      priors$Phi_prior <- matrix(0, nrow=k, ncol=m)
      if (verbose) message("Phi_prior was not specified, set to 0 [k x m]")
    }
    
    
        
    # extract priors from list for simplier notation
    v_prior <- priors$v_prior
    S_prior <- priors$S_prior
    Omega_prior <- priors$Omega_prior
    Phi_prior <- priors$Phi_prior
    

    # Phi|Sigma ~ MN(Phi_prior,Sigma o Omega_prior)
    # Sigma ~IW(S_prior, v_prior)
    
    # convinient short-cuts
    XtX <- t(X) %*% X
    
    if (verbose) message("Calculating (XtX)^{-1}...")
    XtX_inv <- sym_inv(XtX)
    if (!is.null(attr(XtX_inv,"Moore-Penrose"))) {
      if (verbose) message("The XtX matrix is so ugly... kappa(XtX) = ",kappa(XtX),". I will use the Moore-Penrose inverse :)")
    }
    # calculate posterior hyperparameters
    v_post <- v_prior + T
    
    if (verbose) message("Calculating Omega_prior^{-1}...")
    if (is.diagonal(Omega_prior)) { # if Omega_prior is diagonal we may accept Inf on diagonal
      Omega_prior_inv <- matrix(0, nrow=k, ncol=k)
      diag(Omega_prior_inv) <- 1/diag(Omega_prior)
    } else { # sym_inv cannot deal with Inf on the diagonal
      Omega_prior_inv <- sym_inv(Omega_prior)
      if (!is.null(attr(Omega_prior_inv,"Moore-Penrose"))) {
        if (verbose) message("The Omega_prior matrix is so ugly... kappa(Omega_prior) = ",kappa(Omega_prior),". I will use the Moore-Penrose inverse :)")
      }
    }
    
    if (verbose) message("Calculating Omega_post...")
    Omega_post <- sym_inv(Omega_prior_inv+XtX)
    if (!is.null(attr(Omega_post,"Moore-Penrose"))) {
      if (verbose) message("The (Omega_prior_inv+XtX) matrix is so ugly... kappa(Omega_prior_inv+XtX) = ",kappa(Omega_prior_inv+XtX),". I will use the Moore-Penrose inverse :)")
    }

    
    
    # here was a mistake :)
    if (verbose) message("Calculating Phi_post...")
    Phi_post <- Omega_post %*% (Omega_prior_inv %*% Phi_prior + t(X) %*% Y)
    
    colnames(Phi_post) <- paste0("eq_",endo_varnames)
    rownames(Phi_post) <- c(paste0(rep(endo_varnames, p),", l=",rep(1:p,each=m)) , exo_varnames)
    
    
    Phi_hat <- XtX_inv %*% t(X) %*% Y 
    E_hat <- Y - X %*% Phi_hat
    
    # Karlsson, p 15
    if (verbose) message("Calculating S_post...")
    S_post <- S_prior + t(E_hat) %*% E_hat +
      t(Phi_prior - Phi_hat) %*% 
                sym_inv(Omega_prior + XtX_inv) %*% 
                         (Phi_prior - Phi_hat)
    
    
    # Carriero, p 51 (mistake should be S_post=S_0+...)
    # S_post <- S_prior + t(E_hat) %*% E_hat +
    #             t(Phi_hat) %*% XtX %*% Phi_hat + 
    #               t(Phi_prior) %*% solve(Omega_prior) %*% Phi_prior -
    #                  t(Phi_post) %*% solve(Omega_post) %*% Phi_post
    
    
    
    
    if (fast_forecast) { # no simulations
      answer <- "Option 'fast_forecast' is TRUE. Hyperparameters available in attr(.,'posterior')."
    } 
    
    if (!fast_forecast) { # simulate phi, sigma
      
      # reserve space for Gibbs sampling replications
      answer <- matrix(0, nrow=keep, ncol = m*k + m*m)
      
      # precalculate ~(Omega_post)^{1/2} for faster cycle
      # way 1:
      if (way_omega_post_root=="cholesky") {
        if (verbose) message("Calculating ~'Omega_post^{1/2}' using chol(Omega_post)...")
        Omega_post_root <- t(chol(Omega_post))
      }
      
      # way 2:
      # (Omega_post)^{1/2} is equal to (Omega_prior^{-1}+X'X)^{-1/2}
      # But Omega_prior^{-1}+X'X maybe replaced using augmented X* as:
      # X*'X* = Omega_prior^{-1}+X'X
      # so we nedd to find
      # (Omega_post)^{1/2} = (X*'X*)^{-1/2}
      # http://math.stackexchange.com/questions/106774
      if (way_omega_post_root=="svd") {
        if ("Omega_prior_m05" %in% names(priors)) {
          Omega_prior_m05 <- priors$Omega_prior_m05
        } else {
          Opm_eigen <- eigen(Omega_prior_inv)
          Omega_prior_m05 <- Opm_eigen$vectors %*% diag(Opm_eigen$values) %*% t(Opm_eigen$vectors)
        }
        X_star <- rbind(Omega_prior_m05,X)
        X_star_svd <- svd(X_star)
        Omega_post_root <- X_star_svd$v %*% diag(1/X_star_svd$d) %*% t(X_star_svd$v)
      }
      
      
      for (i in 1:keep) {
        if ((verbose) & (i %% 10^3 == 0)) message("Iteration ",i," out of ",keep)
        
        Sigma <- MCMCpack::riwish(v_post,S_post) 
        
        
        # fast way to generate matrix-variate normal
        V <- matrix(rnorm((m*p+d)*m), ncol = m) # [(mp+d) x m] matrix of standard normal
        Phi <- Phi_post + Omega_post_root %*% V %*% chol(Sigma)
        
        Phi_vec <- as.vector(Phi)
        Sigma_vec <- as.vector(Sigma) # length = m x m
        
        answer[i,] <- c(Phi_vec, Sigma_vec)
      } # cycle from 1 to keep
      
      # save as mcmc object to make some good functions available
      answer <- coda::as.mcmc(answer)
    } # simulations 
    
    # set prior attributes:
    attr(answer, "params")  <- data.frame(k=k,m=m,p=p,d=d, 
                                          T_in=T_in,T=T,T_dummy=T_dummy,
                                          constant=constant,
                                          keep=keep, fast_forecast=fast_forecast)
    

    attr(answer, "data") <- list(Y_in=Y_in, Z_in=Z_in, 
                                 X_dummy=priors$X_dummy, Y_dummy=priors$Y_dummy,
                                 X_wo_dummy=X_wo_dummy, Y_wo_dummy=Y_wo_dummy,
                                 endo_varnames=endo_varnames,
                                 exo_varnames=exo_varnames)
    
    priors$type <- "conjugate"
    attr(answer, "prior")       <- priors
    attr(answer, "posterior") <- list(Omega_post=Omega_post,
                                      v_post=v_post,
                                      Phi_post=Phi_post,
                                      S_post=S_post)

    return(answer)
}

#' predict with conjugate Normal-Inverse-Wishart bayesian VAR model
#'
#' predict with conjugate Normal-Inverse-Wishart bayesian VAR model
#'  
#' predict with conjugate Normal-Inverse-Wishart bayesian VAR model
#' 
#' @param model estimated conjugate N-IW model
#' @param h number of periods for forecasting
#' @param level confidence levels for prediction intervals
#' @param Y_in (NULL by default) past values of endogeneous variables (shold have at least p observations).
#' If NULL, then Y_in supplied for estimation will be used. For out-of-sample forecast only last p values of Y_in are used
#' @param Z_f future values of exogeneous variables
#' @param type ("prediction" by default) type of interval: "prediction" incorporates uncertainty about
#' future shocks; "credible" deals only with parameter uncertainty.
#' @param output (default "long") --- long or wide table for mean/quantiles of forecasts
#' @param out_of_sample logical, default is TRUE, whether forecasts are out-of-sample or in-sample.
#' If forecasts are not out-of-sample, then parameter h is ignored
#' @param include (default is c("mean", "median", "sd", "raw")) what type of summary to provide
#' If "raw" is present raw forecasts will be reported. If only "raw" is present
#' then function will return coda mcmc object with raw forecasts. 
#' Otherwise raw forecasts will be attached as attribute.
#' @param fast_forecast logical, FALSE by default. If TRUE then only mean forecast is calculated,
#' posterior expected values of hyperparameters are used. No confidence intervals, no sd, no median. 
#' This mode is activated by default if there are no simulations in supplied model.
#' @param verbose (default FALSE) if true some messages will be printed
#' @export
#' @return forecast results
#' @examples 
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors, keep=100)
#' forecast_conjugate(model, h=2, output="wide")
#' forecast_conjugate(model, out_of_sample = FALSE, include="mean", level=NULL, type = "credible")
forecast_conjugate <- function(model, 
                               Y_in=NULL, Z_f=NULL,
                               output=c("long","wide"),
                               h=1, out_of_sample=TRUE,
                               type=c("prediction","credible"), level=c(80,95),
                               include=c("mean","median","sd","interval", "raw"),
                               fast_forecast=FALSE,
                               verbose=FALSE) {

  
  # select type of prediction specified by user
  type <- match.arg(type)
  output <- match.arg(output)
  
  # simplify notation, extract params from attribute
  T <- attr(model,"params")$T # number of observations minus p plus T_dummy
  p <- attr(model,"params")$p
  k <- attr(model,"params")$k
  m <- attr(model,"params")$m
  d <- attr(model,"params")$d
  T_dummy <- attr(model,"params")$T_dummy
  
  keep <- attr(model, "params")$keep

  constant <- attr(model,"params")$constant
  
  if ((attr(model,"params")$fast_forecast) & (!fast_forecast) ){
    if (verbose) message("No simulations found in 'model', fast_forecast option is set to TRUE.")
    fast_forecast <- TRUE
  }
  
  if (fast_forecast) {
    level <- NULL
    include <- "mean"
    keep <- 1
    
    # type of confidence interval does not make sense for fast mean forecast
    # but setting type to 'credible' avoids simulations of future epsilon, so
    type <- "credible"
  }
  
  
  # if Y_in is not supplied take Y_in from estimation
  if (is.null(Y_in)) Y_in <- attr(model, "data")$Y_in
  
  
  # in case of in-sample forecast h is set to T-T_dummy, or is is better set to NA?
  if (!out_of_sample) h <- T-T_dummy
  
  
  # sanity check
  if (!is.null(Z_f))
    if (!nrow(Z_f)==h) stop("I need exactly h = ",h," observations for exogeneous variables.")
  
  if (nrow(Y_in)<p) stop("Model has ",p," lags. To predict I need at least ",p," observations, but only ",nrow(Y_in)," are provided.")
  
  
  # if requested add constant to exogeneous regressors
  if (constant) Z_f <- cbind(rep(1, h), Z_f)
  

  if (out_of_sample) {
    # take last p observations of Y_in (out-of-sample forecast)
    Y_in <- tail(Y_in, p)
  } else {
    # take first p observations of Y_in (in-sample forecast)
    Y_in <- head(Y_in, p)
  }
  
  e_t <- rep(0, m) # ok for bayesian credible intervals  
  
  # space to store all forecasted values
  forecast_raw <- matrix(0, nrow = keep, ncol = m*h)
  
  
  
  x_t <- rep(0, k)
  
  for (i in 1:keep) {
    # forecast h steps for given sampling of Phi
    
    if (fast_forecast) { # use posterior expected value of Phi
      Phi <- attr(model, "posterior")$Phi_post
    } else { # use Phi^[s] and Sigma^[s]
      Phi <- matrix(model[i,1:(k*m)], nrow=k)
      
      
      Sigma <- matrix(model[i,(k*m+1):(m*k + m*m)],nrow=m) # Sigma [m x m]
      # find square root of draw from Sigma (code is part of mvtnorm function)
      ev <- eigen( Sigma, symmetric = TRUE)
      if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
        warning("Omega_post is numerically not positive definite")
      }
      # precalculate R to do less operations in case h>1
      R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    }
    
    Phi_transp <- t(Phi) # precalculate to do less operations in case h>1
    
    for (j in 1:h) {
      
      
      if (out_of_sample) { # out-of-sample forecast
        # fill exogeneous values in x_t
        x_t[(m * p + 1):(m * p + d)] <- Z_f[j,]
        
        # fill endogeneous values in x_t (first out-of-sample forecast)
        if (j == 1)
          x_t[1:(m * p)] <- as.vector(t(Y_in)[,p:1])
        
        # fill endogeneous values recursively (second+ out-of-sample forecast)
        if (j > 1) {
          if (p>1) x_t[(m + 1):(m * p)] <- x_t[1:(m * (p - 1))]
          x_t[1:m] <- y_t
        }
      } else { # in-sample forecast
        x_t <- attr(model,"data")$X_wo_dummy[j,]
      }
      
      if (type=="prediction") e_t <- R %*% rnorm(m) 
      # e_t is 0 for bayesian credible intervals
      
      y_t <- Phi_transp %*% x_t + e_t
      forecast_raw[i, (m*(j-1)+1):(m*j)] <- y_t
    } # end_for j in 1:h
  } # end_for i in 1:keep
  
  # save as mcmc object for standartisation
  forecast_raw <- coda::as.mcmc(forecast_raw)
  
  # we have m endogeneous variables and h forecasts for each
  id_block <- data.frame(variable=rep( attr(model ,"data")$endo_varnames, h), h=rep(1:h, each=m))
  
  forecast_summary <- NULL
  
  # calculate mean
  if ("mean" %in% include) {
    what <- rep("mean", h*m)
    value <- apply(forecast_raw, 2, mean)
    block <- cbind(id_block, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }

  # calculate median
  if ("median" %in% include) {
    what <- rep("median", h*m)
    value <- apply(forecast_raw, 2, median)
    block <- cbind(id_block, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }
  
  # sd
  if ("sd" %in% include) {
    what <- rep("sd", h*m)
    value <- apply(forecast_raw, 2, sd)
    block <- cbind(id_block, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }

  # calculate quantiles
  if ("interval" %in% include) {
    for (lev in level) {
      # lower
      what <- rep(paste0("lower_",lev), h*m)
      value <- apply(forecast_raw, 2, function(x) quantile(x, probs=(1-lev/100)/2))
      block <- cbind(id_block, what, value) # block of information
      forecast_summary <- rbind(forecast_summary, block)
      
      # upper
      what <- rep(paste0("upper_",lev), h*m)
      value <- apply(forecast_raw, 2, function(x) quantile(x, probs=(1+lev/100)/2))
      block <- cbind(id_block, what, value) # block of information
      forecast_summary <- rbind(forecast_summary, block)
    }
  }
  
  rownames(forecast_summary)  <- NULL
  
  
  if ((output=="wide") & (!is.null(forecast_summary))) { # transform to wide format if requested
    forecast_summary <- reshape2::dcast(forecast_summary, variable+h~what)
  }
  
  
  if ("raw" %in% include) {
    if (length(include)==1) { # only raw forecasts are requested
      forecast_summary <- forecast_raw
    } else { # attach raw forecasts as attribute
      attr(forecast_summary, "forecast_raw") <- forecast_raw
    }
  }
  return(forecast_summary)
}


#' summary of a conjugate Normal-Inverse-Wishart bayesian VAR model
#'
#' summary of a conjugate Normal-Inverse-Wishart bayesian VAR model
#'  
#' summary of a conjugate Normal-Inverse-Wishart bayesian VAR model
#' 
#' @param model estimated conjugate N-IW model
#' @export
#' @return nothing
#' @examples 
#' data(Yraw)
#' setup <- bvar_conj_setup(Yraw, p = 4)
#' model <- bvar_conj_estimate(setup = setup, keep=100)
#' bvar_conj_summary(model)
bvar_conj_summary <- function(model) {
  T <- nrow(model$X) # number of observations minus p
  p <- model$p
  k <- ncol(model$X)
  m <- ncol(model$Y)
  constant <- model$constant
  d <- k - m*p
  
  keep <- 0 
  if ("sample" %in% names(model)) keep <- nrow(model$sample)
  T_in <- T + p
  T_dummy <- nrow(model$X_plus)

  message("Number of lags, p = ", p)
  message("Number of endogeneos variables, m = ",m)
  
  const_message <- "(no constant is present)"
  if (constant) const_message <- "(including constant)"
  message("Number of exogeneos variables ", const_message ,", d = ",d)
  
  message("Number of parameters, k = mp + d = ",k)
  message("Initial number of observations, T_in = ",T_in)
  message("Number of dummy observations, T_dummy = ", T_dummy )
  message("Number of observations available for classic VAR, T = T_in - p = ",T)
  
  message("Posterior mean of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
  print(model$Phi_post)
  
  message("Posterior nu = ",model$v_post)
  
  if (keep>0) {
    message("Number of mcmc simulations, keep = ", keep) 
    
    post_mean <- apply(model$sample, 2, mean)
    post_sd <- apply(model$sample, 2, sd)
    
    message("Posterior sample mean of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
    Phi_post_sample_mean <- matrix(head(post_mean, k*m), nrow=k)
    colnames(Phi_post_sample_mean) <- colnames(model$Phi_post)
    rownames(Phi_post_sample_mean) <- rownames(model$Phi_post)
    print(Phi_post_sample_mean)
    
    message("Posterior sample mean of Sigma (noise covariance matrix) [m = ",m," x m = ",m,"]:")
    Sigma_post_sample_mean <- matrix(tail(post_mean, m*m), nrow=m)
    colnames(Sigma_post_sample_mean) <- bvar_get_endo_varnames(model$Y)
    rownames(Sigma_post_sample_mean) <- bvar_get_endo_varnames(model$Y)
    print(Sigma_post_sample_mean)
    
    
    message("Posterior sample sd of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
    Phi_post_sample_sd <- matrix(head(post_sd, k*m), nrow=k)
    colnames(Phi_post_sample_sd) <- colnames(model$Phi_post)
    rownames(Phi_post_sample_sd) <- rownames(model$Phi_post)
    print(Phi_post_sample_sd)
    
    message("Posterior sample sd of Sigma (noise covariance matrix) [m = ",m," x m = ",m,"]:")
    Sigma_post_sample_sd <- matrix(tail(post_sd, m*m), nrow=m)
    colnames(Sigma_post_sample_sd) <- bvar_get_endo_varnames(model$Y)
    rownames(Sigma_post_sample_sd) <- bvar_get_endo_varnames(model$Y)
    print(Sigma_post_sample_sd)
    
    
    
  } else {
    message("No mcmc simulations were done, only posteriors were calculated.")
  }
}

#' summary of a conjugate Normal-Inverse-Wishart bayesian VAR model
#'
#' summary of a conjugate Normal-Inverse-Wishart bayesian VAR model
#'  
#' summary of a conjugate Normal-Inverse-Wishart bayesian VAR model
#' 
#' @param model estimated conjugate N-IW model
#' @export
#' @return nothing
#' @examples 
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors, keep=100)
#' summary_conjugate(model)
summary_conjugate <- function(model) {
  T <- attr(model,"params")$T # number of observations minus p
  p <- attr(model,"params")$p
  k <- attr(model,"params")$k
  m <- attr(model,"params")$m
  d <- attr(model,"params")$d
  keep <- attr(model, "params")$keep
  T_in <- attr(model, "params")$T_in
  T_dummy <- attr(model, "params")$T_dummy
  fast_forecast <- attr(model, "params")$fast_forecast
  
  message("Number of lags, p = ", p)
  message("Number of endogeneos variables, m = ",m)
  message("Number of exogeneos variables (including constant), d = ",d)
  message("Number of parameters, k = mp + d = ",k)
  message("Initial number of observations, T_in = ",T_in)
  message("Number of dummy observations, T_dummy = ", T_dummy )
  message("Number of observations available for regression, T = T_in + T_dummy - p = ",T)
  
  message("Posterior mean of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
  print(attr(model,"posterior")$Phi_post)
  
  message("Posterior nu = ",attr(model,"posterior")$v_post)
  
  if (fast_forecast) {
    message("'fast_forecast' option is TRUE. Only posterior hyperparameters are calculated.")
  } else {
    message("Number of mcmc simulations, keep = ", keep) 
    post_mean <- apply(model, 2, mean)
    post_sd <- apply(model, 2, sd)
    
    message("Posterior sample mean of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
    Phi_post_sample_mean <- matrix(head(post_mean, k*m), nrow=k)
    colnames(Phi_post_sample_mean) <- colnames(attr(model,"posterior")$Phi_post)
    rownames(Phi_post_sample_mean) <- rownames(attr(model,"posterior")$Phi_post)
    print(Phi_post_sample_mean)
    
    message("Posterior sample mean of Sigma (noise covariance matrix) [m = ",m," x m = ",m,"]:")
    Sigma_post_sample_mean <- matrix(tail(post_mean, m*m), nrow=m)
    colnames(Sigma_post_sample_mean) <- attr(model,"data")$endo_varnames
    rownames(Sigma_post_sample_mean) <- attr(model,"data")$endo_varnames
    print(Sigma_post_sample_mean)

    
    message("Posterior sample sd of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
    Phi_post_sample_sd <- matrix(head(post_sd, k*m), nrow=k)
    colnames(Phi_post_sample_sd) <- colnames(attr(model,"posterior")$Phi_post)
    rownames(Phi_post_sample_sd) <- rownames(attr(model,"posterior")$Phi_post)
    print(Phi_post_sample_sd)

    message("Posterior sample sd of Sigma (noise covariance matrix) [m = ",m," x m = ",m,"]:")
    Sigma_post_sample_sd <- matrix(tail(post_sd, m*m), nrow=m)
    colnames(Sigma_post_sample_sd) <- attr(model,"data")$endo_varnames
    rownames(Sigma_post_sample_sd) <- attr(model,"data")$endo_varnames
    print(Sigma_post_sample_sd)

  }
    
}


#' Multivariate log-gamma-function
#'
#' Multivariate log-gamma-function
#'  
#' Multivariate log-gamma-function. Wikipedia-style parametrisation (p,a)
#' 
#' @param p dimension 
#' @param a argument 
#' @export
#' @return log G_p(a), log of multivariate gamma-function
#' @examples 
#' lmvgamma(2,3)
lmvgamma <- function(p,a) { # wikipedia parameters
  res <- p*(p-1)/4*log(pi) + sum( lgamma(a+(1-(1:p))/2) )
  return(res)
}


#' Calculate log marginal data density
#'
#' Calculate log marginal data density
#'  
#' Calculate log marginal data density. Formula from Carriero p. 55
#' 
#' @param model estimated conjugate N-IW model
#' @export
#' @return log of marginal data density
#' @examples 
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors, fast_forecast = TRUE)
#' marginal_data_density(model)
marginal_data_density <- function(model) {

  params <- attr(model,"params")
  data <- attr(model,"data")
  prior <- attr(model,"prior")
  post <- attr(model,"post")
  
  p <- params$p
  k <- params$k
  m <- params$m
  d <- params$d
  keep <- params$keep
  T_in <- params$T_in
  T_dummy <- params$T_dummy
  fast_forecast <- params$fast_forecast
  Omega_prior <- prior$Omega_prior
  Phi_prior <- prior$Phi_prior
  S_prior <- prior$S_prior
  v_prior <- prior$v_prior
  
  
  # what are exactly T, X, Y and v_post (it uses T). Shall we include dummy variables?!
  X <- data$X_wo_dummy # here we may experiment :)
  Y <- data$Y_wo_dummy
  
  T <- nrow(X)
  v_post <- v_prior + T 
  
 
  I_XoX <- diag(T) + X %*% Omega_prior %*% t(X)  # diag() = Identity matrix
  I_XoX_inv <- sym_inv(I_XoX)
  e_prior <- Y-X %*% Phi_prior    
  
  # just to avoid a very long line:
  line_1 <- -T*m/2 * log(pi) - m/2 * log(det(I_XoX)) + v_prior/2 * log(det(S_prior)) 
  line_2 <- -v_post/2 * log(det( t(e_prior) %*% I_XoX_inv %*% e_prior ))
  res <- lmvgamma(m, v_post/2) - lmvgamma(m, v_prior/2) + line_1 + line_2
  return(res)
  
}



#' Calculate log marginal data density
#'
#' Calculate log marginal data density
#'  
#' Calculate log marginal data density. Formula from Carriero p. 55
#' 
#' @param model estimated conjugate N-IW model
#' @export
#' @return log of marginal data density
#' @examples 
#' data(Yraw)
#' setup <- bvar_conj_setup(Yraw, p = 4)
#' model <- bvar_conj_estimate(setup = setup, keep=100)
#' bvar_conj_mdd(model)
bvar_conj_mdd <- function(model) {
  
  T <- nrow(model$X) # number of observations minus p
  p <- model$p
  k <- ncol(model$X)
  m <- ncol(model$Y)
  constant <- model$constant
  d <- k - m*p
  
  
  # old: # prior <- bvar_conj_dummy2hyper(model$Y_plus, model$X_plus)
  # model was estimated using additional dummy observations
  # but we use Omega_prior, Phi_prior and S_prior directly obtained from lambdas
  # this is more robust, than obtaining them from calculated dummy observations
  Omega_prior <- model$Omega_prior # prior$Omega
  Phi_prior <- model$Phi_prior # prior$Phi
  S_prior <- model$S_prior # prior$S
  
  v_prior <- model$v_prior
  
  X <- model$X 
  Y <- model$Y

  v_post <- model$v_post
  
  
  I_XoX <- diag(T) + X %*% Omega_prior %*% t(X)  # diag() = Identity matrix
  I_XoX_inv <- sym_inv(I_XoX)
  e_prior <- Y-X %*% Phi_prior    
  
  # just to avoid a very long line:
  line_1 <- -T*m/2 * log(pi) - m/2 * log(det(I_XoX)) + v_prior/2 * log(det(S_prior)) 
  line_2 <- -v_post/2 * log(det(S_prior + t(e_prior) %*% I_XoX_inv %*% e_prior )) 
  res <- lmvgamma(m, v_post/2) - lmvgamma(m, v_prior/2) + line_1 + line_2
  return(res)
  
}


