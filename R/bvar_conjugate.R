#' Set conjugate N-IW priors from lambdas as in Carriero
#' 
#' Set conjugate N-IW priors from lambdas as in Carriero
#' 
#' Set conjugate N-IW priors from lambdas as in Carriero
#' Based on Carriero p. 52-53
#'
#' @param p number of lags
#' @param Y_in multivariate time series
#' @param lambdas vector = (l_1, l_power, l_sc, l_io, l_const, l_exo), the l_2 is set to 1 automatically for 
#' conjugate N-IW prior. Short summary:
#' sd(const in eq i) = l_const * sigma_i
#' sd(exo in eq i)= l_exo * sigma_i
#' sd(coef for var j lag l in eq i) = l_1*sigma_i/sigma_j/l^l_power
#' @param Z_in exogeneous variables
#' @param s2_lag number of lags in AR() model used to estimate s2 (equal to p by default)
#' Carriero uses 1 in his matlab code
#' @param VAR_in (either "levels" or "growth rates"). 
#' If set to "levels" (default) Phi_1 matrix is identity, if "growth rates" then to zero matrix.
#' @param dummy_sc whether to include "sum of coefficients" dummies, logical (TRUE by default)
#' @param dummy_io whether to include "initial observation" dummies, logical (TRUE by default)
#' @return priors list containing Phi_prior [k x m], Omega_prior [k x k], S_prior [m x m], v_prior [1x1],
#' where k = mp+d
#' @export
#' @examples 
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4, lambdas = c(0.2,1,1,1,100,100))
#' model <- bvar_conjugate0(priors = priors)
Carriero_priors <- function(Y_in, Z_in=NULL, constant=TRUE, p=4, 
                            lambdas=c(0.2,1,1,1,100,100), 
                            VAR_in=c("levels","growth rates"), 
                            s2_lag=NULL, 
                            dummy_io=TRUE, dummy_sc=TRUE) {
  l_1 <- lambdas[1]
  l_2 <- 1
  l_power <- lambdas[2]
  l_sc <- lambdas[3]
  l_io <- lambdas[4]
  l_const <- lambdas[5]
  l_exo <- lambdas[6]
  
  Y_in <- as.matrix(Y_in) # to clear tbl_df if present :)
  
  # calculate d, the number of exogeneous regressors
  if (is.null(Z_in)) {
    d <- 1*constant
  } else {
    d <- ncol(Z_in) + 1*constant
  }
  
  # if requested add constant to exogeneous regressors
  # constant to the left of other exo variables
  if (constant) Z <- cbind(rep(1, nrow(Y_in)), Z_in)
  
  
  m <- ncol(Y_in)
  k <- m*p+d
  
  VAR_in <- match.arg(VAR_in)
  
  if (!l_2==1) warning("Conjugate N-IW is impossible for lambda_2 <> 1")
  
  # Litterman takes 6 lags in AR(p)
  
  # estimate sigma^2 from univariate AR(p) processes
  if (is.null(s2_lag)) s2_lag <- p
  
  sigmas_sq <- rep(NA, m)
  for (j in 1:m) {
    
    y_uni <- Y_in[,j] # univariate time series

    
    # old version: it fails when ML estimation fails :)
    #AR_p <- forecast::Arima(y_uni, order = c(p,0,0), method="ML") # AR(p) model
    #sigmas_sq[j] <- AR_p$sigma2
    
    
    # Carriero matlab code: always AR(1)! make an option?
    
    # more robust version: fails only in the case of  severe multicollinearity
    AR_p <- ar.ols(y_uni, aic=FALSE, order.max = s2_lag) # AR(p) model
    resid <- tail(AR_p$resid,-s2_lag) # omit first p NA in residuals
    sigmas_sq[j] <- sum(resid^2)/(length(resid)-s2_lag-1)
  }
  
  # set Phi_prior
  if (VAR_in=="levels") Phi_1 <- diag(m)
  if (VAR_in=="growth rates") Phi_1 <- matrix(0, m,m)
  Phi_prior <- t( cbind(Phi_1, matrix(0, nrow=m, ncol=k-m)) )
  
  S_prior <- diag(sigmas_sq)
  v_prior <- m+2
  
  # set Omega_prior
  # the diagonal of Omega_prior begins with endogeneous part:
  endo_diagonal <- l_1^2*rep(1/sigmas_sq, p)/rep(1/(1:p)^(2*l_power), each=m)
  
  # and ends with exogeneous part:
  exo_diagonal <- rep(l_exo^2,d)
  if (constant) exo_diagonal[1] <- l_const^2
  
  Omega_diagonal <- c(endo_diagonal, exo_diagonal)
  # and set zero prior covariances
  Omega_prior <- diag(Omega_diagonal)
  
  
  
  
  # create dummy observations
  y_0_bar <- apply(as.matrix(Y_in[1:p,],nrow=p), 2, mean) # vector [m x 1] of mean values of each endo-series
  if (is.null(Z)) z_bar <- NULL # special case of no constant and no exo vars
  if (!is.null(Z)) z_bar <- 
    apply(as.matrix(Z[1:p,],nrow=p), 2, mean) # vector [d x 1] of mean values of each exo-series
  # "as.matrix" above is needed to avoid errors for p=1 or d=1
  
  
  
  # sum of coefficients prior
  Y_dummy_sc <- NULL
  X_dummy_sc <- NULL
  if (dummy_sc) {
    Y_dummy_sc <- matrix(0, m, m) # zero matrix [m x m]
    diag(Y_dummy_sc) <- y_0_bar / l_sc
  
    X_dummy_sc <- matrix(0, m, k) # zero matrix [m x k]
    # X_dummy_sc is not a square matrix, 
    # but diag() will correctly fill "diagonal" elements, X_dummy[i,i]
    diag(X_dummy_sc) <- y_0_bar / l_sc
  }
  
  # dummy initial observation
  Y_dummy_io <- NULL
  X_dummy_io <- NULL
  if (dummy_io) {
    Y_dummy_io <- matrix(y_0_bar/l_io, nrow=1)
    X_dummy_io <- matrix(c(rep(y_0_bar/l_io, p), z_bar/l_io), nrow=1)
  }
  
  # order of dummies does not matter
  X_dummy <- rbind(X_dummy_io, X_dummy_sc)
  Y_dummy <- rbind(Y_dummy_io, Y_dummy_sc)
  
    
  priors <- list(v_prior=v_prior, S_prior=S_prior, 
                 Phi_prior=Phi_prior, Omega_prior=Omega_prior, 
                 Y_dummy=Y_dummy, X_dummy=X_dummy,
                 Y_in=Y_in, Z_in=Z_in, p=p, # to avoid duplicating
                 sigmas_sq = sigmas_sq ) # get more info from function
  
  return(priors)
}





#' Set conjugate N-IW priors as in matlab code of Koops-Korobilis
#' 
#' Set conjugate N-IW priors as in matlab code of Koops-Korobilis
#' 
#' Set conjugate N-IW priors as in matlab code of Koops-Korobilis
#'
#' @param p number of lags
#' @param Y_in multivariate time series
#' @param Z_in exogeneous variables
#' @return priors list containing Phi_prior [k x m], Omega_prior [k x k], S_prior [m x m], v_prior [1x1],
#' where k = mp+d
#' @export
#' @examples 
#' data(Yraw)
#' priors <- KK_code_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors)
KK_code_priors <- function(Y_in, Z_in=NULL, constant=TRUE, p=4) {

  # calculate d, the number of exogeneous regressors
  if (is.null(Z_in)) {
    d <- 1*constant
  } else {
    d <- ncol(Z_in) + 1*constant
  }
  
  # if requested add constant to exogeneous regressors
  if (constant) Z <- cbind(rep(1, nrow(Y_in)), Z_in)
  
  
  m <- ncol(Y_in)
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
                 Y_in=Y_in, Z=Z_in, p=p)
  
  return(priors)
}




#' Set conjugate N-IW priors from lambdas and mus as in Sim Zha
#' 
#' Set conjugate N-IW priors from lambdas and mus as in Sim Zha
#' 
#' Set conjugate N-IW priors from lambdas and mus as in Sim Zha
#' Should be compatible with szbvar function. Maybe error!!!! 
#' MAYBE lambda should be in denominator!!!!
#'
#' @param p number of lags
#' @param Y_in multivariate time series
#' @param lambdas vector = (l0, l1, l3, l4, l5), l2 is set to 1 for conjugate N-IW
#' @param mu56 vector = (mu5, mu6)
#' @param Z_in exogeneous variables
#' @param VAR_in (either "levels" or "growth rates")
#' @return priors list containing Phi_prior [k x m], Omega_prior [k x k], S_prior [m x m], v_prior [1x1],
#' where k = mp+d
#' @export
#' @examples 
#' data(Yraw)
#' priors <- szbvar_priors(Yraw, p = 4, lambdas = c(1,0.2,1,1,1), mu56=c(1,1))
#' model <- bvar_conjugate0(priors = priors)
szbvar_priors <- function(Y_in, Z_in=NULL, constant=TRUE, p=4, 
                          lambdas=c(1,0.2,1,1,1,1), mu56=c(1,1),
                            VAR_in=c("levels","growth rates")) {
  l0 <- lambdas[1]
  l1 <- lambdas[2]
  l2 <- 1
  l3 <- lambdas[3]
  l4 <- lambdas[4]
  l5 <- lambdas[5]
  mu5 <- mu56[1]
  mu6 <- mu56[2]
  
  Y_in <- as.matrix(Y_in)
  # calculate d, the number of exogeneous regressors
  if (is.null(Z_in)) {
    d <- 1*constant
  } else {
    d <- ncol(Z_in) + 1*constant
  }
  
  # if requested add constant to exogeneous regressors
  if (constant) Z <- cbind(rep(1, nrow(Y_in)), Z_in)
  
  
  m <- ncol(Y_in)
  k <- m*p+d
  
  VAR_in <- match.arg(VAR_in)
  
  
  # create dummy observations
  y_0_bar <- apply(as.matrix(Y_in[1:p,],nrow=p), 2, mean) # vector [m x 1] of mean values of each endo-series
  
  z_bar <- NULL
  if (!is.null(Z)) z_bar <- # avoid error with null Z and no constant
    apply(as.matrix(Z[1:p,],nrow=p), 2, mean) # vector [d x 1] of mean values of each exo-series
  # "as.matrix" above is needed to avoid errors for p=1 or d=1
  
  
  # sum of coefficients prior
  Y_dummy_sc <- matrix(0, m, m) # zero matrix [m x m]
  diag(Y_dummy_sc) <- y_0_bar * mu5
  
  X_dummy_sc <- matrix(0, m, k) # zero matrix [m x k]
  # X_dummy_sc is not a square matrix, 
  # but diag() will correctly fill "diagonal" elements, X_dummy[i,i]
  diag(X_dummy_sc) <- y_0_bar * mu5
  
  # dummy initial observation
  Y_dummy_io <- matrix(y_0_bar * mu6, nrow=1)
  X_dummy_io <- matrix(c(rep(y_0_bar * mu6, p), z_bar * mu6), nrow=1)
  
  
  # order of dummies??? does not matter??? # dummy_sc, dummy_io in MSBVAR::szbvar
  X_dummy <- rbind(X_dummy_sc, X_dummy_io)
  Y_dummy <- rbind(Y_dummy_sc, Y_dummy_io)
  
  
  message("MAYBE a bug and lambdas should be inverted!!!")
  
  
  if (!l2==1) warning("Conjugate N-IW is impossible for lambda_2 <> 1")
  
  # Litterman takes 6 lags in AR(p)
  
  # estimate sigma^2 from univariate AR(p) processes
  sigmas_sq <- rep(NA, m)
  
  Y <- rbind(Y_dummy,tail(Y_in,-p))
  
  for (j in 1:m) {
    # y_uni <- Y_in[,j] # univariate time series (original Y_in)
    # VERY strange, but:
    y_uni <- Y[,j] # univariate time series (Y instead of Y_in)
    
    # here MSBVAR::szbvar uses Y (= Y_in without first p obs and augmented with dummy obs)
    
    # AR_p <- forecast::Arima(y_uni, order = c(p,0,0)) # AR(p) model
    # sigmas_sq[j] <- AR_p$sigma2
    
    # in MSBVAR:
    sigmas_sq[j] <- ar.ols(y_uni, aic = FALSE, order.max = p, 
           intercept = TRUE, demean = FALSE)$var.pred
    
  }
  
  # set Phi_prior
  if (VAR_in=="levels") Phi_1 <- diag(m)
  if (VAR_in=="growth rates") Phi_1 <- matrix(0, m,m)
  Phi_prior <- t( cbind(Phi_1, matrix(0, nrow=m, ncol=k-m)) )
  
  # S_prior <- diag(m) # identity matrix [m x m]
  S_prior <- diag(sigmas_sq)*l0^2 # as in MSBVAR::szbvar
  
  v_prior <- m+1
  
  lag_power <- 1/(1:p)^l3
  
  # set Omega_prior
  Omega_diagonal <- c(kronecker(1/lag_power^2, 1/sigmas_sq), l0^2*l4^2, rep(l0^2*l5^2, d-1))
  # and set zero prior covariances
  Omega_prior <- diag(Omega_diagonal)
  
  
  
  

  
  
  priors <- list(v_prior=v_prior, S_prior=S_prior, 
                 Phi_prior=Phi_prior, Omega_prior=Omega_prior, 
                 Y_dummy=Y_dummy, X_dummy=X_dummy,
                 Y_in=Y_in, Z_in=Z_in, p=p, # to avoid duplicating
                 sigmas_sq=sigmas_sq, Y=Y) # get more info from function
  
  return(priors)
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
  A_chol <- chol(A)
  inv_A <- chol2inv(A_chol)
  return(inv_A)
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
#' @param verbose (FALSE by default)
#' @param priors the list containing Phi_prior [k x m], Omega_prior [k x k], 
#' S_prior [m x m], v_prior [1x1],
#' Y_dummy [T_dummy x m], X_dummy [T_dummy x k]
#' where k = mp+d
#' @return the list containing all results of bayesian VAR estimation
#' @export
#' @examples
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors)
bvar_conjugate0 <-
  function(Y_in=NULL, Z_in=NULL, constant=TRUE, p=NULL, keep=10000, verbose=FALSE,
           priors=list(Phi_prior=NULL, Omega_prior=NULL, S_prior=NULL, v_prior=NULL, 
                       Y_dummy=NULL, X_dummy=NULL, Y_in=NULL, Z_in=NULL, p=NULL) ) {

    if ( (is.null(Y_in)) & (!is.null(priors$Y_in)) ){
      Y_in <- priors$Y_in
      message("Y_in is inferred from priors data.")
    }
    
    if ( (is.null(Z_in)) & (!is.null(priors$Z_in)) ){
      Z_in <- priors$Z_in
      message("Z_in is inferred from priors data.")
    }

    if ( (is.null(p)) & (!is.null(priors$p)) ) {
      p <- priors$p
      message("Number of lags is inferred from priors data: p = ",p)
    }
    
    if ( (is.null(p)) & (is.null(priors$p)) ) {
      p <- 4
      message("Number of lags, p, is not specified inside and outside priors, set to p = ",p)
    }
    
    
    
    
    # if Z_in is provided it should have the same number of rows that Y_in
    if (!is.null(Z_in)) 
      if (!nrow(Y_in)==nrow(Z_in))
        stop("Number of rows in Y_in and Z_in should be equal. 
           The p rows of Z_in are not used and may be filled with NA")
    
    
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
        
    # get dimensions
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
    }

    # set some bad priors for lazy guys if not supplied 
    if (is.null(priors$v_prior)) {
      priors$v_prior <- m + 1
      message("v_prior was not specified, set to (m+1)")
    }
    if (is.null(priors$S_prior)) {
      priors$S_prior <- diag(m)
      message("S_prior was not specified, set to I [m x m]")
    }
    if (is.null(priors$Omega_prior)) {
      priors$Omega_prior <- 10*diag(k) 
      message("Omega_prior was not specified, set to 10I [k x k]")
    }  
    if (is.null(priors$Phi_prior)) {
      priors$Phi_prior <- matrix(0, nrow=k, ncol=m)
      message("Phi_prior was not specified, set to 0 [k x m]")
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
    XtX_inv <- try(solve(XtX), silent=TRUE)
    if (class(XtX_inv)=="try-error") {
      message("The XtX matrix is so ugly... :( \n I will use the Moore-Penrose inverse :) \n kappa(XtX) = ",kappa(XtX))
      XtX_inv <- MASS::ginv(XtX) # solve(XtX) # for more stable results in multicollinearity cases
    }
    # calculate posterior hyperparameters
    v_post <- v_prior + T
    Omega_post <- sym_inv(sym_inv(Omega_prior)+XtX)
    
    # here was a mistake :)
    Phi_post <- Omega_post %*% (sym_inv(Omega_prior) %*% Phi_prior + t(X) %*% Y)
    
    Phi_hat <- XtX_inv %*% t(X) %*% Y 
    E_hat <- Y - X %*% Phi_hat
    
    # Karlsson, p 15
    S_post <- S_prior + t(E_hat) %*% E_hat +
      t(Phi_prior - Phi_hat) %*% 
                sym_inv(Omega_prior + XtX_inv) %*% 
                         (Phi_prior - Phi_hat)
    
    
    # Carriero, p 51 (mistake should be S_post=S_0+...)
    # S_post <- S_prior + t(E_hat) %*% E_hat +
    #             t(Phi_hat) %*% XtX %*% Phi_hat + 
    #               t(Phi_prior) %*% solve(Omega_prior) %*% Phi_prior -
    #                  t(Phi_post) %*% solve(Omega_post) %*% Phi_post
    
    # reserve space for Gibbs sampling replications
    answer <- matrix(0, nrow=keep, ncol = m*k + m*m)
    
    # precalculate chol(Omega_post) for faster cycle
    chol_Omega_post <- chol(Omega_post)
    
    for (i in 1:keep) {
      if ((verbose) & (i %% 10^3 == 0)) message("Iteration ",i," out of ",keep)
      
      Sigma <- MCMCpack::riwish(v_post,S_post) 
      
      # slow way
      # Phi_vec <- mvtnorm::rmvnorm(n = 1, # vec(Phi) ~ N(vec(Phi_post), Sigma o Omega_post)
      #                            mean = as.vector(Phi_post),
      #                            sigma = kronecker(Sigma, Omega_post)) 
      
      # Phi_vec has length (mp+d) x m
      # Phi <- matrix(Phi_vec, ncol=m) # we are saving only vector Phi_vec
      
      # fast way, Carriero, p. 54
      V <- matrix(rnorm((m*p+d)*m), ncol = m) # [(mp+d) x m] matrix of standard normal
      Phi <- Phi_post + chol_Omega_post %*% V %*% t(chol(Sigma))

      Phi_vec <- as.vector(Phi)
      Sigma_vec <- as.vector(Sigma) # length = m x m
      
      answer[i,] <- c(Phi_vec, Sigma_vec)
    }
    
    # save as mcmc object to make some good functions available
    answer <- coda::as.mcmc(answer)
    
    # set prior attributes:
    attr(answer, "params")  <- data.frame(k=k,m=m,p=p,d=d, 
                                          T_in=T_in,T=T,T_dummy=T_dummy,
                                          constant=constant,
                                          keep=keep)
    
    attr(answer, "data") <- list(Y_in=Y_in, Z_in=Z_in, 
                                 X_dummy=priors$X_dummy, Y_dummy=priors$Y_dummy,
                                 X_wo_dummy=X_wo_dummy, Y_wo_dummy=Y_wo_dummy)
    
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
#' If NULL, then Y_in supplied for estimation will be used. Only last p values of Y_in are used.
#' @param Z_f future values of exogeneous variables
#' @param type ("prediction" by default) type of interval: "prediction" incorporates uncertainty about
#' future shocks; "credible" deals only with parameter uncertainty.
#' @param output (default "long") --- long or wide table for mean/quantiles of forecasts
#' @param out_of_sample logical, default is TRUE, whether forecasts are out of sample or not.
#' If forecasts are not out of sample, then parameter h is ignored
#' @param include (default is c("mean", "median", "sd")) what type of summary to provide
#' If include is NULL and level is NULL then the function will return raw mcmc predictions
#' @export
#' @return forecast results
#' @examples 
#' data(Yraw)
#' priors <- Carriero_priors(Yraw, p = 4)
#' model <- bvar_conjugate0(priors = priors)
#' forecast_conjugate(model, h=2, output="wide")
forecast_conjugate <- function(model, 
                               Y_in=NULL, 
                               Z_f=NULL,
                               output=c("long","wide"),
                               h=1, level=c(80,95),
                               type=c("prediction","credible"),
                               out_of_sample=TRUE,
                               include=c("mean","median","sd")) {

  
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
  
  # if Y_in is not supplied take Y_in from estimation
  if (is.null(Y_in)) Y_in <- attr(model, "data")$Y_in
  
  
  # in case of in-sample forecast h is set to T-T_dummy
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
    Phi <- matrix(model[i,1:(k*m)], nrow=k)
    Phi_transp <- t(Phi) # precalculate to do less operations in case h>1
    
    Sigma <- matrix(model[i,(k*m+1):(m*k + m*m)],nrow=m) # Sigma [m x m]
    # find square root of draw from Sigma (code is part of mvtnorm function)
    ev <- eigen( Sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("Omega_post is numerically not positive definite")
    }
    # precalculate R to do less operations in case h>1
    R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    
    
    for (j in 1:h) {
      
      
      if (out_of_sample) { # out-of-sample forecast
        # fill exogeneous values in x_t
        x_t[(m * p + 1):(m * p + d)] <- Z_f[j,]
        
        # fill endogeneous values in x_t (first out-of-sample forecast)
        if (j == 1)
          x_t[1:(m * p)] <- as.vector(t(Y_in)[,p:1])
        
        # fill endogeneous values recursively (second+ out-of-sample forecast)
        if (j > 1) {
          x_t[(m + 1):(m * p)] <- x_t[1:(m * (p - 1))]
          x_t[1:m] <- y_t
        }
      } else { # in-sample forecast
        x_t <- attr(model,"data")$X_wo_dummy[j,]
      }
      
      if (type=="prediction") e_t <- R %*% rnorm(m) 
      # e_t is 0 for bayesian credible intervals
      
      y_t <- Phi_transp %*% x_t + e_t
      forecast_raw[i, (m*(j-1)+1):(m*j)] <- y_t
    }
  }
  # save as mcmc object for standartisation
  forecast_raw <- coda::as.mcmc(forecast_raw)
  
  # we have m endogeneous variables and h forecasts for each
  varnames <- data.frame(y=rep(1:m, h), h=rep(1:h, each=m))
  
  forecast_summary <- NULL
  
  # calculate mean
  if ("mean" %in% include) {
    what <- rep("mean", h*m)
    value <- apply(forecast_raw, 2, mean)
    block <- cbind(varnames, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }

  # calculate median
  if ("median" %in% include) {
    what <- rep("median", h*m)
    value <- apply(forecast_raw, 2, median)
    block <- cbind(varnames, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }
  
  # sd
  if ("sd" %in% include) {
    what <- rep("sd", h*m)
    value <- apply(forecast_raw, 2, sd)
    block <- cbind(varnames, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }

  # calculate quantiles
  for (lev in level) {
    # lower
    what <- rep(paste0("lower_",lev), h*m)
    value <- apply(forecast_raw, 2, function(x) quantile(x, probs=(1-lev/100)/2))
    block <- cbind(varnames, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
    
    # upper
    what <- rep(paste0("upper_",lev), h*m)
    value <- apply(forecast_raw, 2, function(x) quantile(x, probs=(1+lev/100)/2))
    block <- cbind(varnames, what, value) # block of information
    forecast_summary <- rbind(forecast_summary, block)
  }
  
  rownames(forecast_summary)  <- NULL
  
  
  if (output=="wide") { # transform to wide format if requested
    forecast_summary <- reshape2::dcast(forecast_summary, y+h~what)
  }
  
  

  if (is.null(include) & is.null(level)) {
    # report only raw forecasts
    forecast_summary <- forecast_raw
  } else {
    # save raw forecasts for further analysis
    attr(forecast_summary, "forecast_raw") <- forecast_raw
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
#' priors <- Carriero_priors(Yraw, p = 4, lambdas = c(1,0.2,1,1))
#' model <- bvar_conjugate0(priors = priors)
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
  
  message("Number of lags, p = ", p)
  message("Number of endogeneos variables, m = ",m)
  message("Number of exogeneos variables (including constant), d = ",d)
  message("Number of parameters, k = mp + d = ",k)
  message("Initial number of observations, T_in = ",T_in)
  message("Number of dummy observations, T_dummy = ", T_dummy )
  message("Number of observations available for regression, T = T_in + T_dummy - p = ",T)
  
  post_mean <- apply(model, 2, mean)
  post_sd <- apply(model, 2, sd)
  
  message("Posterior mean of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
  print(matrix(head(post_mean, k*m), nrow=k))
  
  message("Posterior mean of Sigma (noise covariance) [m = ",m," x m = ",m,"]:")
  print(matrix(tail(post_mean, m*m), nrow=m))
    
  message("Posterior sd of Phi (VAR coefficients) [k = ",k," x m = ",m,"]:")
  print(matrix(head(post_sd, k*m), nrow=k))
  
  message("Posterior sd of Sigma (noise covariance) [m = ",m," x m = ",m,"]:")
  print(matrix(tail(post_sd, m*m), nrow=m))
  
}


