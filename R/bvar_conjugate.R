#' Set conjugate N-IW priors from lambdas
#' 
#' Set conjugate N-IW priors from lambdas
#' 
#' Set conjugate N-IW priors from lambdas
#' Based on Carriero p. 53
#'
#' @param p number of lags
#' @param Y multivariate time series
#' @param lambdas vector = (l0, l1, l2, l3, l4)
#' @param d (default is 1) number of exogeneous variables
#' @param VAR_in (either "levels" or "growth rates")
#' @return priors list containing Phi_prior [k x m], Omega_prior [k x k], S_prior [m x m], v_prior [1x1],
#' where k = mp+d
#' @export
#' @examples 
#' data(Yraw)
#' priors <- lambda2priors(Yraw)
lambda2priors <- function(Y, p=4, d=1, lambdas=c(1,0.2,1,1,1), 
                          VAR_in=c("levels","growth rates")) {
  l0 <- lambdas[1]
  l1 <- lambdas[2]
  l2 <- lambdas[3]
  l3 <- lambdas[4]
  l4 <- lambdas[5]
  
  m <- ncol(Y)
  k <- m*p+d
  
  VAR_in <- match.arg(VAR_in)
  
  if (!l2==1) warning("Conjugate N-IW is impossible for lambda_2 <> 1")
  
  # Litterman takes 6 lags in AR(p)
  
  # estimate sigma^2 from univariate AR(p) processes
  sigmas_sq <- rep(NA, m)
  for (j in 1:m) {
    y_uni <- Y[,j] # univariate time series
    AR_p <- forecast::Arima(y_uni, order = c(p,0,0)) # AR(p) model
    sigmas_sq[j] <- AR_p$sigma2
  }
  
  # set Phi_prior
  if (VAR_in=="levels") Phi_1 <- diag(m)
  if (VAR_in=="growth rates") Phi_1 <- matrix(0, m,m)
  Phi_prior <- t( cbind(Phi_1, matrix(0, nrow=k-m, ncol=m)) )
  
  S_prior <- diag(sigmas_sq)
  v_prior <- m+2
  
  # set Omega_prior
  
  s2i_ <- matrix(sigmas_sq, nrow = m, ncol = m)
  s2_j <- matrix(sigmas_sq, nrow = m, ncol = m, byrow = TRUE)
  
  # first we calculate prior variance of each element of Phi
  # they are located in [k x m] matrix like Phi itself
  Phi_vars <- NULL
  for (b in 1:p) {
    var_block <- l1*l2*s2i_/s2_j/b
    Phi_vars <- cbind(Phi_vars, var_block)
  }
  Phi_vars <- cbind(Phi_vars, l0*s2i_ )
  Phi_vars <- t(Phi_vars)
  
  # we vectorize Phi_vars
  Omega_diagonal <- as.vector(Phi_vars)
  # and set zero prior covariances
  Omega_prior <- diag(Omega_diagonal)
  
  
  
  
  
  priors <- list(v_prior=v_prior, S_prior=S_prior, 
                 Phi_prior=Phi_prior, Omega_prior=Omega_prior, Y_dummy=Y_dummy, X_dummy=X_dummy)
  
  return(priors)
}






#' Estimate Normal-Inverse-Wishart bayesian VAR model
#'
#' Estimate Normal-Inverse-Wishart bayesian VAR model
#'  
#' Estimate Normal-Inverse-Wishart bayesian VAR model
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
#' model <- bvar_conjugate0(Y)
bvar_conjugate0 <-
  function(Y_in, Z_in=NULL, constant=TRUE, p=4, keep=10000, verbose=FALSE,
           priors=list(Phi_prior=NULL, Omega_prior=NULL, S_prior=NULL, v_prior=NULL, 
                       Y_dummy=NULL, X_dummy=NULL) ) {

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
    
    
    
    # dimension of deterministic regressors
    d <- ncol(Z)
    
    # here we add dummy observations
    T_dummy <- 0 
    if (!is.null(priors$Y_dummy)) {
      if (!nrow(priors$Y_dummy)==nrow(priors$X_dummy)) stop("X_dummy and Y_dummy should have the same number of rows")
      T_dummy <- nrow(priors$Y_dummy)
    }
    Y <- rbind(priors$Y_dummy, Y)
    X <- rbind(priors$X_dummy, X)
        
    # get dimensions
    T <- nrow(Y)
    m <- ncol(Y)
    k <- m*p + d
    
    if (verbose) { 
      message("Number of lags, p =", p)
      message("Number of endogeneos variables, m = ",m)
      message("Number of exogeneos variables (including constant), d = ",d)
      message("Number of parameters, k = mp + d =",k)
      message("Initial number of observations, T_in = ",T_in)
      message("Number of dummy observations, T_dummy = ", T_dummy )
      message("Number of observations available for regression, T = T_in + T_dummy - p = ",T)
    }
    
    # extract priors from list for simplier notation
    v_prior <- priors$v_prior
    S_prior <- priors$S_prior
    Omega_prior <- priors$Omega_prior
    Phi_prior <- priors$Phi_prior
    
    # set some bad priors for lazy guys if not supplied 
    if (is.null(v_prior)) v_prior <- m + 1
    if (is.null(S_prior)) S_prior <- diag(m)
    if (is.null(Omega_prior)) Omega_prior <- 10*diag(k) 
    if (is.null(Phi_prior)) Phi_prior <- matrix(0, nrow=k, ncol=m)
    

    # Phi|Sigma ~ MN(Phi_prior,Sigma o Omega_prior)
    # Sigma ~IW(S_prior, v_prior)
    
    # convinient short-cuts
    XtX <- t(X) %*% X
    XtX_inv <- solve(XtX)
    
    # calculate posterior hyperparameters
    v_post <- v_prior + T
    Omega_post <- solve(solve(Omega_prior)+XtX)
    Phi_post <- Omega_post %*% (solve(Omega_prior) %*% Phi_prior)
    
    Phi_hat <- XtX_inv %*% t(X) %*% Y 
    E_hat <- Y - X %*% Phi_hat
    
    # Karlsson, p 15
    S_post <- S_prior + t(E_hat) %*% E_hat +
      t(Phi_prior - Phi_hat) %*% 
                solve(Omega_prior + XtX_inv) %*% 
                         (Phi_prior - Phi_hat)
    
    
    # Carriero, p 51 (mistake should be S_post=S_0+...)
    # S_post <- S_prior + t(E_hat) %*% E_hat +
    #             t(Phi_hat) %*% XtX %*% Phi_hat + 
    #               t(Phi_prior) %*% solve(Omega_prior) %*% Phi_prior -
    #                  t(Phi_post) %*% solve(Omega_post) %*% Phi_post
    
    # reserve space for Gibbs sampling replications
    answer <- matrix(0, nrow=keep, ncol = m*(m*p+d) + m*m)
    
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
    attr(answer, "prior")       <- list(type="conjugate",
                                   Phi_prior=Phi_prior,
                                   v_prior=v_prior,
                                   Omega_prior=Omega_prior,
                                   S_prior=S_prior)

    return(answer)
}

