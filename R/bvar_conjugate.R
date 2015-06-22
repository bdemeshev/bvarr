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
#' @param priors the list containing Phi_prior [k x m], Omega_prior [k x k], S_prior [m x m], v_prior [1x1],
#' where k = mp+d
#' @return the list containing all results of bayesian VAR estimation
#' @export
#' @examples
#' model <- bvar_conjugate0(Y)
bvar_conjugate0 <-
  function(Y_in, Z_in=NULL, constant=TRUE, p=2, keep=10000, verbose=FALSE,
           priors=list(Phi_prior=NULL, Omega_prior=NULL, S_prior=NULL, v_prior=NULL) ) {

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
    
    # here we add dummy observations...?
    # ...
        
    # get dimensions
    T <- nrow(Y)
    m <- ncol(Y)
    k <- m*p + d
    
    
    # extract priors from list
    v_prior <- priors$v_prior
    S_prior <- priors$S_prior
    Omega_prior <- priors$Omega_prior
    Phi_prior <- priors$Phi_prior
    
    # set some bad priors for lazy guys if not supplied 
    if (is.null(v_prior)) v_prior <- m + 2
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

