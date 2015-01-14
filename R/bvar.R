#' bvar
#'
#' @name bvar
#' @docType package
#' @author Gary Koop, Dimitris Korobilis, Boris Demeshev 
#' @import MCMCpack mvtnorm reshape2 ggplot2
NULL



#' US inflation, employement and interest rate data.frame
#'
#' Quarterly US data on inflation, unemployment and interest rate, 
#' 1953:Q1 - 2006:Q3 The variables are as follows:
#'
#' \itemize{
#' \item inflation inflation
#' \item employment employment
#' \item interest_rate interest rate
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Yraw
#' @usage data(Yraw)
#' @format A data frame with ... rows and 3 variables
NULL



#' Estimate six types of bayesian VAR models
#'
#' Estimate six types of bayesian VAR models
#'
#' @param Yraw the matrix or data.frame with endogeneous VAR variables
#' @param prior the type of prior: "diffuse","minnesota","conjugate", "independent","ssvs-Wishart","ssvs-ssvs"
#' @param W the matrix of exogenous variables
#' @param p the number of lags for endogeneous variables
#' @param constant (TRUE/FALSE) indicator, whether the constant is included
#' @param nsave the length of mcmc chain
#' @param nburn the length of burn-in of mcmc chain
#' @param it_print which iteration are printed
#' @param impulses whether to calculate irfs
#' @param ihor the number of lags for irfs
#' @param forecasting whether to calculate forecasts
#' @param repfor number of times to obtain a draw from the predictive
#' @param h number of forecast periods
#' @param a_i ssvs-ssvs hyperparameter, diagonal elements of SIGMA (Gamma density)
#' @param b_i ssvs-ssvs hyperparameter, diagonal elements of SIGMA (Gamma density)
#' @param kappa_0 ssvs-ssvs hyperparameter, variances for non-diagonal elements of SIGMA
#' @param kappa_1 ssvs-ssvs hyperparameter, variances for non-diagonal elements of SIGMA
#' @param p_i ssvs-ssvs hyperparameter, for Gamma ~ BERNOULLI(m,p_i), see eq. (14)
#' @param q_ij ssvs-ssvs  hyperparameter, for Omega_[j] ~ BERNOULLI(j,q_ij), see eq. (17)
#' @param a_bar minnesota hyperparameter vector (own lag, other lag, exogeneuos)
#' @return the list containing all results of bayesian VAR estimation
#' @export
#' @examples
#' bvar(Yraw)
bvar <-
function(Yraw, prior = "ssvs-ssvs", W = NULL, p = 4, constant = TRUE,
                 nsave = 10000, nburn = 2000, it_print = 2000, # gibbs-related
                 impulses = TRUE, ihor = 24,  # impulses-related
                 forecasting = TRUE, repfor = 50, h = 1, # forecasting-related
                 p_i = 0.5, q_ij = 0.5, # S-S hyperparameters 
                 kappa_0 = 0.1, kappa_1 = 6,   # S-S hyperparameters         
                 a_i = 0.01, b_i = 0.01,  # S-S hyperparameters 
                 a_bar = c(0.5,0.5,10^2) # minnesota hyperparameters
         ) {



# BVAR_FULL.m 
# This code replicates the results from the 1st empirical illustration 
# in Koop and Korobilis (2009).
#
# You can chose 6 different priors. For some priors, analytical results are
# available, so Monte Carlo Integration is used. For other priors, you need
# to use the Gibbs sampler. For Gibbs sampler models I take a number of
# 'burn-in draws', so that I keep only the draws which have converged.
#
# The specification of the prior hyperparmeters are in the file
# prior_hyper.m. See there for details.
#
# The convention used here is that ALPHA is the K x M matrix of VAR coefficients,
# alpha is the KM x 1 column vector of vectorized VAR coefficients, i.e.
# alpha = vec(ALPHA), and SIGMA is the M x M VAR covariance matrix.
#--------------------------------------------------------------------------
# Bayesian estimation, prediction and impulse response analysis in VAR
# models using posterior simulation. Dependent on your choice of forecasting,
# the VAR model is:
#
# In this code we provide direct (as opposed to iterated) forecasts
# Direct h-step ahead foreacsts:
#     Y(t+h) = A0 + Y(t) x A1 + ... + Y(t-p+1) x Ap + e(t+h)
#
# so that in this case there are also p lags of Y (from 0 to p-1).
#
# In any of the two cases, the model is written as:
#
#                   Y(t) = X(t) x A + e(t)
#
# where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
# also use the vector a which is defined as a=vec(A).
#--------------------------------------------------------------------------
# NOTES: The code sacrifices efficiency for clarity. It follows the
#        theoretical equations in the monograph and the manual.
#
# AUTHORS: Gary Koop and Dimitris Korobilis
# CONTACT: dikorombilis@yahoo.gr
#--------------------------------------------------------------------------
  


# [Yraw] = bvardgp(); # BB: translate to R!!!!

# Note that 'Yraw' is a matrix with T rows by M columns,
# where T is the number of time series observations (usually months or
# quarters), while M is the number of VAR dependent macro variables.


# p - vector of probabilities for 0
bernoullirnd <- function(p) {
  len <- length(p)
  ans <- ifelse(runif(len)<p,0,1)
  return(ans)
}


impulse <- function(By,smat,nstep) {
  # function response=impulse(By,smat,nstep), C. Sims' code.
  # smat is a square matrix of initial shock vectors.  To produce "orthogonalized
  # impulse responses" it should have the property that smat'*smat=sigma, where sigma
  # is the Var(u(t)) matrix and u(t) is the residual vector.  One way to get such a smat
  # is to set smat=chol(sigma).  To get the smat corresponding to a different ordering,
  # use smat=chol(P*Sigma*P')*P, where P is a permutation matrix.
  # By is a neq x nvar x nlags matrix.  neq=nvar, of course, but the first index runs over 
  # equations. In response, the first index runs over variables, the second over 
  # shocks (in effect, equations).
  
  neq <- dim(By)[1]
  nvar <- dim(By)[2]
  nlag <- dim(By)[3]

  resp <- array(0, dim=c(nvar,neq,nstep)) 
  resp[,,1] <- t(smat) # need lower triangular, last innovation untransformed
  for (it in 2:nstep) {
    for (ilag in 1:min(nlag,it-1)) {
      resp[,,it] <- resp[,,it]+By[,,ilag] %*% resp[,,it-ilag]
    }
  }
  return(resp)
}




# used for mlag2 function
makeblock <- function(X,i,p) {
  Xblock <- X[(p+1-i):(nrow(X)-i),] # get useful lines of X
  Xblock <- as.matrix(Xblock) # assure X is a matrix, not a vector
  Xblock <- rbind(matrix(0,nrow=p,ncol=ncol(X)),Xblock)  # append p zero lines at top
  return(Xblock)
}

mlag2 <- function(X,p) {
  X <- as.matrix(X)
  
  # we need to bind horizontally p blocks
  Xlag <- matrix(nrow=nrow(X),ncol=0) # create empty matrix with correct number of raws
  for (i in 1:p) Xlag <- cbind(Xlag,makeblock(X,i,p)) # bind blocks horizontally
  return(Xlag)  
}


# even if the prior is not ssvs we'll have these objects in output
gammas_mat <- NULL
omega_mat <- NULL
gamma_draws <- NULL
omega_draws <- NULL


# -------- BB: GO-GO-GO

# For models using analytical results, there are no convergence issues 
# (You are not adviced to change the next 3 lines)
if (prior %in% c("diffuse","minnesota","conjugate") ) nburn <- 0

ntot <- nsave + nburn  # Total number of draws
 
Yraw <- as.matrix(Yraw)
# Get initial dimensions of dependent variable
Traw <- nrow(Yraw)
M <- ncol(Yraw)


if (forecasting & h<=0) 
  stop("You have set forecasting, but the forecast horizon h<=0")
possible.priors <- c("diffuse","minnesota","conjugate",
                     "independent","ssvs-wishart","ssvs-ssvs")
prior <- tolower(prior)
if (!prior %in% possible.priors) stop("Possible priors are: ", paste(possible.priors,collapse=", ")) 

# The model specification is different when implementing direct forecasts,
# compared to the specification when computing iterated forecasts.
if (forecasting) {
  # Now create VAR specification according to forecast method

  Y1 <- Yraw[(h+1):Traw,]
  Y2 <- Yraw[2:(Traw-h),]
  Traw <- Traw - h - 1
} else {
  Y1 <- Yraw
  Y2 <- Yraw
}

# Generate lagged Y matrix. This will be part of the X matrix
Ylag <- mlag2(Y2,p) # Y is [T x M]. ylag is [T x (Mp)]

# print(dimnames(Ylag))

# Now define matrix X which has all the R.H.S. variables (constant, lags of
# the dependent variable and exogenous regressors/dummies).
# Note that in this example I do not include exogenous variables (other macro
# variables, dummies, or trends). You can load a file with exogenous
# variables, call them, say W, and then extend variable X1 in line 133, as:
#            X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:) W(p+1:Traw,:)];
# and line 135 as:
#            X1 = [Ylag(p+1:Traw,:)  W(p+1:Traw,:)];

X1 <- Ylag[(p+1):Traw,]
if (!is.null(W)) X1 <- cbind(X1,W[(p+1):Traw,]) # append exogenous variables 
if (constant) { 
  X1 <- cbind(1, X1) # add column of ones to the left of Ylag 
  dimnames(X1)[[2]][1] <- "constant" # name to understand
}



# Get size of final matrix X
Traw3 <- nrow(X1)
K <- ncol(X1)



# Create the block diagonal matrix Z
Z1 <- kronecker(diag(M),X1)

# print(dimnames(Z1))

# Form Y matrix accordingly
# Delete first "LAGS" rows to match the dimensions of X matrix
Y1 <- Y1[(p+1):Traw,] # This is the final Y matrix used for the VAR

# Traw was the dimesnion of the initial data. T is the number of actual 
# time series observations of Y and X
T <- Traw - p

#========= FORECASTING SET-UP:
# Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if (forecasting) {
  Y_pred <- matrix(0,nsave*repfor,M) # Matrix to save prediction draws
  colnames(Y_pred) <- colnames(Yraw)
  
  PL <- rep(0,nsave)                 # BB: in R PL is vector! 
  # BB: (PL was nsave x 1 matrix in matlab)

  # Direct forecasts, we only need to keep the last observation for evaluation
  Y <- head(Y1,-1)                # Y is equal to Y1 without 1 last observation
  X <- head(X1,-1)
  Z <- kronecker(diag(M),X)  
  T <- T - 1
} else { # if no prediction is present, keep all observations
  Y <- Y1
  X <- X1
  Z <- Z1
}

#========= IMPULSE RESPONSES SET-UP:
# Create matrices to store forecasts

if (impulses) {   # BB: REDO? Cleaner separation of estimation/prediction/irf calculation?
  
  # BB: save all impulse responses in standard way
  all_responses <- array(0,c(nsave,M,M,ihor))
  dimnames(all_responses)[[2]] <- colnames(Yraw)
  dimnames(all_responses)[[3]] <- colnames(Yraw)
  # dimnames(all_responses)[[4]] <- paste("lag",1:ihor)
  
  # print(dimnames(all_responses))
  
  bigj <- matrix(0,M,M*p)  
  bigj[1:M,1:M] <- diag(M)
  rownames(bigj) <- colnames(Yraw)
}

#-----------------------------PRELIMINARIES--------------------------------
#First get ML estimators
tXX <- t(X) %*% X # save time!!!
inv_tXX <- solve(tXX)
tXY <- t(X) %*% Y

A_OLS <- inv_tXX %*% tXY # This is the matrix of regression coefficients
a_OLS <- as.vector(A_OLS)         # This is the vector of parameters, i.e. it holds
# that a_OLS = vec(A_OLS)
SSE <- t(Y - X%*%A_OLS)%*%(Y - X%*%A_OLS)   # Sum of squared errors (BB: a matrix!)
SIGMA_OLS <- SSE/(T-K+1)

# Initialize Bayesian posterior parameters using OLS values
alpha <- a_OLS     # This is the single draw from the posterior of alpha
ALPHA <- A_OLS     # This is the single draw from the posterior of ALPHA
SSE_Gibbs <- SSE   # This is the SSE based on each draw of ALPHA
SIGMA <- SIGMA_OLS # This is the single draw from the posterior of SIGMA
IXY <-  kronecker(diag(M),tXY)

# Storage space for posterior draws
alpha_draws <- matrix(0,nsave,K*M)   # save draws of alpha
ALPHA_draws <- array(0,c(nsave,K,M))   # save draws of ALPHA
dimnames(ALPHA_draws)[[3]] <- colnames(Yraw)
SIGMA_draws <- array(0,c(nsave,M,M))   # save draws of SIGMA
dimnames(SIGMA_draws)[[2]] <- colnames(Yraw)
dimnames(SIGMA_draws)[[3]] <- colnames(Yraw)


# input:
# p ---  Number of lags on dependent variables
# M --- dimension of y_t
# output
# priors :)



#if (prior == "diffuse") { # Diffuse (1)
  # I guess there is nothing to specify in this case!
  # Posteriors depend on OLS quantities
#}

if (prior == "minnesota") { # minnesota-Whishart (2)
  # Prior mean on VAR regression coefficients
  # !!!!!!! АХТУНГ!!!!!! ПОХОЖЕ БЕЗ КОНСТАНТЫ НАДО УБРАТЬ rep(0,M) !!!!!!!!!!!!
  A_prior <- rbind(rep(0,M),0.9*diag(M),matrix(0,(p-1)*M,M)) #<---- prior mean of ALPHA (parameter matrix) 
  a_prior <- as.vector(A_prior)               #<---- prior mean of alpha (parameter vector)
  
  # minnesota Variance on VAR regression coefficients
  # First define the hyperparameters 'a_bar_i'
  # BB: moved into function definition # a_bar <- c(0.5,0.5,10^2)
  
  # Now get residual variances of univariate p_MIN-lag autoregressions. Here
  # we just run the AR(p) model on each equation, ignoring the constant
  # and exogenous variables (if they have been specified for the original
  # VAR model)
  p_MIN <- 6
  sigma_sq <- rep(0,M) # vector to store residual variances
  for (i in 1:M) {
    # Create lags of dependent variables   
    Ylag_i <- mlag2(Yraw[,i],p)
    Ylag_i <- Ylag_i[(p_MIN+1):Traw,]
    X_i = cbind(rep(1,Traw-p_MIN),Ylag_i)
    Y_i <- Yraw[(p_MIN+1):Traw,i]
    # OLS estimates of i-th equation
    alpha_i <- solve(t(X_i)%*%X_i)%*%(t(X_i)%*%Y_i)
    sigma_sq[i] <- 1/(Traw-p_MIN)*sum((Y_i - X_i%*%alpha_i)^2)
  }
  # Now define prior hyperparameters.
  # Create an array of dimensions K x M, which will contain the K diagonal
  # elements of the covariance matrix, in each of the M equations.
  V_i <- matrix(0,K,M)
  
  # index in each equation which are the own lags
  ind <- matrix(0,M,p)
  for (i in 1:M)
    ind[i,] <- constant+seq(from=i,by=M,len=p) # matlab, i:M:K 
  # but i:M:K sometime produces sequence longer than p
  
  for (i in 1:M) { # for each i-th equation
    for (j in 1:K) {   # for each j-th RHS variable
      if (constant) {
        if (j==1) V_i[j,i] <- a_bar[3]*sigma_sq[i] # variance on constant                
        if (j %in% ind[i,]) 
          V_i[j,i] <- a_bar[1]/(ceiling((j-1)/M)^2) # variance on own lags  
        if ((j>1)& (!(j %in% ind[i,]))) {         
          ll <- (j-1) %% M # remainder of division by M
          if (ll==0) ll <- M # replace remainder==0 by M
          
          V_i[j,i] <- (a_bar[2]*sigma_sq[i])/((ceiling((j-1)/M)^2)*sigma_sq[ll]) 
        }
      }
      if (!constant) {
        if (j %in% ind[i,])
          V_i[j,i] <- a_bar[1]/(ceiling(j/M)^2) # variance on own lags
        # !!!! was (j-1)/M in the bvar_analyt but j in bvar_full aux script
        if (!(j %in% ind[i,])) {
          
          ll <- j %% M # remainder of division by M
          if (ll==0) ll <- M # replace remainder==0 by M
          
          V_i[j,i] <- (a_bar[2]*sigma_sq[i])/((ceiling(j/M)^2)*sigma_sq[ll])            
          # !!!! was (j-1)/M in the bvar_analyt but j in bvar_full aux script
          
        }
      }
    }} # i,j cycle
  
  # Now V is a diagonal matrix with diagonal elements the V_i
  V_prior <- diag(as.vector(V_i))  # this is the prior variance of the vector alpha
  
  # NOTE: No prior for SIGMA. SIGMA is simply a diagonal matrix with each
  # diagonal element equal to sigma_sq(i). See Kadiyala and Karlsson (1997)
  SIGMA <- diag(sigma_sq)
}    
if (prior == "conjugate") { # Normal-Whishart
  # Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
  A_prior = matrix(0,K,M)          #<---- prior mean of ALPHA (parameter matrix)
  a_prior <- as.vector(A_prior)    #<---- prior mean of alpha (parameter vector)
  V_prior <- 10*diag(K)            #<---- prior variance of alpha
  
  # Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
  v_prior <- M+1                 #<---- prior Degrees of Freedom (DoF) of SIGMA
  S_prior <- diag(M)      #<---- prior scale of SIGMA
  inv_S_prior <- solve(S_prior)     # BB: funny way to obtain identity matrix :)
}    
if (prior == "independent")  { # Independent Normal-Wishart
  n <- K*M # Total number of parameters (size of vector alpha)
  a_prior <- rep(0,n)    #<---- prior mean of alpha (parameter vector)
  V_prior <- 10*diag(n)   #<---- prior variance of alpha
  
  # Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
  v_prior <- M+1             #<---- prior Degrees of Freedom (DoF) of SIGMA
  S_prior <- diag(M)         #<---- prior scale of SIGMA
  inv_S_prior <- solve(S_prior) # BB: funny way to obtain identity matrix :)
}    
if (prior == "ssvs-wishart" | prior == "ssvs-ssvs") { # ssvs on alpha, Wishart or ssvs on SIGMA    
  n <- K*M # Total number of parameters (size of vector alpha)
  # mean of alpha
  a_prior <- rep(0,n) 
  
  # This is the std of the OLS estimate ofalpha. You can use this to 
  # scale tau_0 and tau_1 (see below) if you want.
  sigma_alpha <- sqrt(diag(kronecker(SIGMA,inv_tXX)))
  # otherwise, set ' sigma_alpha <- rep(1,n) '
  
  # ssvs variances for alpha
  tau_0 <- 0.1*sigma_alpha   # Set tau_[0i], tau_[1i]
  tau_1 <- 10*sigma_alpha
  

  if (prior == "ssvs-wishart")  { # Wishart on SIGMA
    # Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
    v_prior <- M+1             #<---- prior Degrees of Freedom (DoF) of SIGMA
    S_prior <- diag(M)         #<---- prior scale of SIGMA  
    inv_S_prior <- solve(S_prior) # BB: inverse of identity :)
  }
  
  
  
  # Initialize Gamma and Omega vectors
  gammas <- rep(1,n)    # vector of Gamma
  omega <- as.list(rep(NA,M-1)) # omega=cell(1,M-1);  # ?????????????
  for (kk_1 in 1:(M-1)) {
    omega[[kk_1]] <- rep(1,kk_1) # omega{kk_1} = ones(kk_1,1);  # Omega_j # ??????????? make a list?
  }
  
  # Set space in memory for some vectors that we are using in ssvs
  gamma_draws <- matrix(0,nsave,n) # vector of gamma draws
  omega_draws <- matrix(0,nsave,0.5*M*(M-1)) # vector of omega draws    
}

#-------------------- Prior specification ends here

# BB: reserve space for vectors and matrices 
# BB: better to use appropriate size
hh <- as.list(rep(NA,M-1)) # cell(1,M-1);
DD_j <- as.list(rep(NA,M-1)) # cell(1,M-1);
eta <- as.list(rep(NA,M-1)) # eta = cell(1,M-1);


psi_ii_sq <- rep(0,M) 


#========================== Start Sampling ================================
#==========================================================================

for (irep in 1:ntot)  { # Start the Gibbs "loop"
  if (irep %% it_print == 0) # print every it_print iterations
     message("Iteration ",irep," out of ",ntot)

  #--------- Draw ALPHA and SIGMA with Diffuse Prior
  if (prior == "diffuse") {
  # Posterior of alpha|SIGMA,Data ~ Normal
    V_post <- kronecker(SIGMA,inv_tXX)
    alpha <- as.vector( a_OLS + t(chol(V_post))%*%rnorm(K*M) ) # Draw alpha
    ALPHA <- matrix(alpha,K,M) # Create draw of ALPHA       
                        
    # Posterior of SIGMA|Data ~ iW(SSE_Gibbs,T-K) 
    SIGMA <- riwish(T-K,SSE_Gibbs) # matlab: inv(wish(inv(SSE_Gibbs),T-K)) # Draw SIGMA                        
  } 
  
  #--------- Draw ALPHA and SIGMA with minnesota Prior
  if (prior == "minnesota") {
    # Draw ALPHA
    for (i in 1:M) {
      V_block_inv <- solve( V_prior[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] )
      V_post <- solve( V_block_inv + tXX / SIGMA[i,i] )
      a_post <- V_post %*% (V_block_inv %*% a_prior[((i-1)*K+1):(i*K)] +  t(X) %*% Y[,i] / SIGMA[i,i] )
      alpha[((i-1)*K+1):(i*K)] <- a_post + t(chol(V_post))%*%rnorm(K) # Draw alpha
    }
    ALPHA <- matrix(alpha,K,M) # Create draw in terms of ALPHA
  
    # SIGMA in this case is a known matrix, whose form is decided in
    # the prior (see prior_hyper.m)
  }
  #--------- Draw ALPHA and SIGMA with Normal-Wishart Prior
  if (prior == "conjugate") {
    # ******Get all the required quantities for the posteriors       
    V_post <- solve( solve(V_prior) + tXX )
    A_post <- V_post %*% (solve(V_prior)%*%A_prior + tXX%*%A_OLS)
    a_post <- as.vector(A_post)
                 
    S_post <- SSE + S_prior + t(A_OLS)%*% tXX %*%A_OLS + t(A_prior)%*%solve(V_prior)%*%A_prior - t(A_post)%*%(solve(V_prior) + tXX)%*%A_post
    v_post <- T + v_prior
    
    # This is the covariance for the posterior density of alpha
    COV <- kronecker(SIGMA,V_post)
    
    # Posterior of alpha|SIGMA,Data ~ Normal
    alpha <- as.vector( a_post + t(chol(COV))%*%rnorm(K*M) ) # Draw alpha
    ALPHA <- matrix(alpha,nrow=K,ncol=M) # Draw of ALPHA

    # Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
    SIGMA <- riwish(v_post,S_post) # inv(wish(inv(S_post),v_post));% Draw SIGMA
  }

  if (prior == "independent") {
    VARIANCE <- kronecker(solve(SIGMA),diag(T))
    V_post <- solve(V_prior + t(Z) %*% VARIANCE %*% Z)
    a_post <- V_post %*% (V_prior %*% a_prior + t(Z) %*% VARIANCE %*% as.vector(Y) )
    alpha <- as.vector( a_post + t(chol(V_post)) %*% rnorm(n) ) # Draw of alpha
        
    ALPHA <- matrix(alpha,nrow=K,ncol=M) # Draw of ALPHA
        
    # Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
    v_post <- T + v_prior
    S_post <- S_prior + t(Y - X %*% ALPHA) %*% (Y - X %*% ALPHA)
    SIGMA <- riwish(v_post,S_post) # Draw SIGMA        
  }
  
  # --------- Draw ALPHA and SIGMA using ssvs prior 
  if (prior == "ssvs-wishart" | prior == "ssvs-ssvs") {
    # Draw SIGMA
    if (prior == "ssvs-wishart") { # Wishart
      # Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
      v_post <- T + v_prior
      S_post <- inv_S_prior + t(Y - X %*% ALPHA) %*% (Y - X %*% ALPHA)
      SIGMA <- riwish(v_post,S_post) # Draw SIGMA
    }
    
    
  if (prior == "ssvs-ssvs") { # ssvs
    # Draw psi|alpha,gamma,omega,DATA from the GAMMA dist.
    
    for (kk_5 in 1:(M-1)) {
      
      het <- ifelse(omega[[kk_5]]==0,kappa_0,kappa_1)
      hh[[kk_5]] <- as.vector(het)
      
      
      # BB: diag(6) gives 6*6 matrix in R but 1*1 matrix in matlab
      # D_j[[kk_5]] <- diag(hh[[kk_5]],nrow=length(hh[[kk_5]]))
      # Now create covariance matrix D_[j]*R_[j]*D_[j], see eq. (15)
      # DD_j[[kk_5]] <- D_j[[kk_5]] %*% D_j[[kk_5]]
      DD_j[[kk_5]] <- diag(hh[[kk_5]]^2,nrow=length(hh[[kk_5]]))
    }
    
    # Create B_[i] matrix
    B_i <- rep(0,M)
    
    # BB: rr=1 is a special case
    B_i[1] <- b_i + 0.5*SSE[1,1]
    psi_ii_sq[1] <- rgamma(1,shape=a_i + 0.5*T,rate=B_i[1])
    
    # BB: rr=2..M
    for (rr in 2:M) {
      s_i <- SSE_Gibbs[1:(rr - 1),rr]
      S_i <- SSE_Gibbs[1:(rr-1),1:(rr-1)]
      
      DiDi <- DD_j[[rr-1]]
      
      Delta_j <- solve(S_i + solve(DiDi))
      
      B_i[rr] <- b_i + 0.5*(SSE_Gibbs[rr,rr] - t(s_i) %*% Delta_j %*% s_i)
      psi_ii_sq[rr] <- rgamma(1,shape=a_i + 0.5*T,rate=B_i[rr])
      
      miu_j = - sqrt(psi_ii_sq[rr])*(Delta_j %*% s_i)
      eta[[rr-1]] <- miu_j + t(chol(Delta_j)) %*% rnorm(rr-1)
      
    }
 
    # Draw omega|eta,psi,phi,gamma,omega,DATA from BERNOULLI dist.
    omega_vec <- NULL # temporary vector to store draws of omega
    for (kk_9 in 1:(M-1)) {
      # BB: time waste: omeg_g <- omega[[kk_9]]
      eta_g <- eta[[kk_9]]

      u_ij1 <- (1/kappa_0)*exp(-0.5*(eta_g^2)/(kappa_0^2))*q_ij
      u_ij2 <- (1/kappa_1)*exp(-0.5*(eta_g^2)/(kappa_1^2))*(1-q_ij)
      ost <- u_ij1/(u_ij1 + u_ij2)
      omeg_g <- bernoullirnd(ost)
      omega_vec <- c(omega_vec, omeg_g) #ok<AGROW>
            
      # BB: size(omega(kk_9))=kk_9 # u_[ij1], u_[ij2], see eqs. (32 - 33)                          
      
      omega[[kk_9]] <- omeg_g #ok<AGROW>
    }

    # Create PSI matrix from individual elements of "psi_ii_sq" and "eta"
     
    PSI_ALL <- diag(sqrt(psi_ii_sq))
    
    for (nn_2 in 1:(M-1)) { # Now non-diagonal elements
      PSI_ALL [1:nrow(eta[[nn_2]]),nn_2+1] <- eta[[nn_2]]
    }
    # Create SIGMA
    SIGMA <- solve(PSI_ALL %*% t(PSI_ALL))        
  } # END DRAWING SIGMA, if prior==ssvs

  # Draw alpha              
  # Hyperparameters for alpha|gamma ~ N_[m](0,D*D)

  # h_i is tau_0 if gamma=0 and tau_1 if gamma=1
  # BB: gammas may be only 0 or 1
  h_i <- ifelse(gammas==0,tau_0,tau_1)
  
  D <- diag(h_i) # Create D. Here D=diag(h_i) will also do
  # BB: the simplest way is the TRUE one!
  
  DD <- D %*% D # Prior covariance matrix for Phi_m
  isig <- solve(SIGMA)
  psi_xx <- kronecker(isig,tXX)
  inv_DD <- solve(DD)
  V_post <- solve(psi_xx + inv_DD)
  # a_post = V_post %*% ((psi_xx) %*% a_OLS + (solve(DD)) %*% a_prior)

  visig <- as.vector(isig)
  a_post <- as.vector( V_post %*% (IXY %*% visig + inv_DD %*% a_prior) )
        
  alpha <- as.vector( a_post + t(chol(V_post)) %*% rnorm(n) ) # Draw alpha

  ALPHA <- matrix(alpha, nrow=K,ncol=M) # Draw of ALPHA

  # Draw gamma|phi,psi,eta,omega,DATA from BERNOULLI dist.    
  # BB: in matlab u_i1, u_i2, gst are scalars, in R code --- vectors
  u_i1 <- (1/tau_0)*exp(-0.5*(alpha/tau_0)^2)*p_i          
  u_i2 <- (1/tau_1)*exp(-0.5*(alpha/tau_1)^2)*(1-p_i)
  gst <- u_i1/(u_i1 + u_i2)
  gammas <- bernoullirnd(gst) # ok<AGROW>
  

  # Save new Sum of Squared Errors (SSE) based on draw of ALPHA  
  SSE_Gibbs <- t(Y - X %*% ALPHA) %*% (Y - X %*% ALPHA)
    
  } # end of common part for ssvs and ssvs-w priors
  # =============Estimation ends here

  # ****************************|Predictions, Responses, etc|***************************
    if (irep > nburn) {
      # =========FORECASTING:
      if (forecasting) {
  
        Y_temp <- matrix(0,nrow=repfor,ncol=M)
        # compute 'repfor' predictions for each draw of ALPHA and SIGMA
        for (ii in 1:repfor) {
          X_fore = c(Y[T,], X[T,2:(M*(p-1)+1)] )
          if (constant) X_fore <- c(1,X_fore)
          # Forecast of T+1 conditional on data at time T
          Y_temp[ii,] <- X_fore %*% ALPHA + rnorm(M) %*% chol(SIGMA)
        }
        # Matrix of predictions
        Y_pred[((irep-nburn-1)*repfor+1):((irep-nburn)*repfor),] <- Y_temp
        # Predictive likelihood
        PL[irep-nburn] <- dmvnorm (Y1[T+1,],X[T,] %*% ALPHA,SIGMA)
        #if (PL[irep-nburn] == 0) { # BB: ACHTUNG!!!! VERY STRANGE
          # BB: maybe the reason was to avoid log(0)
          # BB: in R there is no such problem, so these lines were commented out
        #  PL[irep-nburn] <- 1
        # }
  
  } # end forecasting
  #=========Forecasting ends here
  
  #=========IMPULSE RESPONSES:
    if (impulses) {
      # ------------Identification code I:
      Bv <- array(0,dim=c(M,M,p))
      for (i_1 in 1:p) {

        alpha_index <- (1+(i_1-1)*M):(i_1*M) # BB: correction for constant=FALSE
        if (constant) alpha_index <- alpha_index + 1
        Bv[,,i_1] <- ALPHA[alpha_index,]              
      }
  
      # st dev matrix for structural VAR
      shock <- t(chol(SIGMA))
      d <- diag(diag(shock)) # BB: ??????
      shock <- solve(d) %*% shock

      
  # responses <- impulse(Bv,shock,ihor)
  
  # BB: we need to save all responses in standard way
  all_responses[irep-nburn,,,] <- impulse(Bv,shock,ihor) # responses
    
    }
  
  #----- Save draws of the parameters
  alpha_draws[irep-nburn,] <- alpha
  ALPHA_draws[irep-nburn,,] <- ALPHA
  SIGMA_draws[irep-nburn,,] <- SIGMA
  if (prior %in% c("ssvs-wishart","ssvs-ssvs")) { 
    gamma_draws[irep-nburn,] <- gammas #ok<AGROW>
    if (prior == "ssvs-ssvs")
      omega_draws[irep-nburn,] <- omega_vec #ok<AGROW>
  }
  
  } # end saving results
  
  
} # end main Gibbs loop

#====================== End Sampling Posteriors ===========================
  


#Posterior mean of parameters:
ALPHA_mean <- apply(ALPHA_draws,c(2,3),mean) # squeeze(mean(ALPHA_draws,1)) #posterior mean of ALPHA
SIGMA_mean <- apply(SIGMA_draws,c(2,3),mean) # squeeze(mean(SIGMA_draws,1)) #posterior mean of SIGMA

#Posterior standard deviations of parameters:
ALPHA_std <- apply(ALPHA_draws,c(2,3),sd) # squeeze(std(ALPHA_draws,1)) #posterior std of ALPHA
SIGMA_std <- apply(SIGMA_draws,c(2,3),sd) # squeeze(std(SIGMA_draws,1)) #posterior std of SIGMA

#or you can use 'ALPHA_COV = cov(alpha_draws,1);' to get the full
#covariance matrix of the posterior of alpha (of dimensions [KM x KM] )

if (prior %in% c("ssvs-wishart","ssvs-ssvs") ) {   
  # Find average of restriction indices Gamma
  gammas <- apply(gamma_draws,2,mean) # BB: 2 is correct! matlab: mean(gamma_draws,1);
  gammas_mat <- matrix(gammas,nrow=K,ncol=M)
  if (prior == "ssvs-ssvs") {
    # Find average of restriction indices Omega
    omega <- apply(omega_draws,2,mean) # BB: 2 is correct, matlab: mean(omega_draws,1)';
    omega_mat <- matrix(0,nrow=M,ncol=M)
    for (nn_5 in 1:(M-1)) {
      ggg <- omega[((nn_5-1)*(nn_5)/2 + 1):(nn_5*(nn_5+1)/2)]
      omega_mat[1:length(ggg),nn_5+1] <- ggg
    }
  }
}




# prediction part
# mean prediction and log predictive likelihood
  Y_pred_mean <- apply(Y_pred,2,mean) # mean prediction
  Y_pred_std <- apply(Y_pred,2,sd) # std prediction
  log_PL <- mean((log(PL)),1)
  
  #This are the true values of Y at T+h:
  true_value <- Y1[T+1,]
  
  Y_pred_melt <- melt(Y_pred)
  



  answer <- list(all_responses = all_responses,
                 Y_pred_melt = Y_pred_melt,
                 Y_pred_std = Y_pred_std,
                 Y_pred_mean = Y_pred_mean,
                 log_PL = log_PL,
                 true_value = true_value,
                 gammas_mat = gammas_mat, # ssvs-specific
                 omega_mat = omega_mat,   # ssvs-specific
                 ALPHA_std = ALPHA_std,
                 SIGMA_std = SIGMA_std,
                 ALPHA_mean = ALPHA_mean,
                 SIGMA_mean = SIGMA_mean,
                 alpha_draws = alpha_draws,
                 ALPHA_draws = ALPHA_draws,
                 SIGMA_draws = SIGMA_draws,
                 gamma_draws = gamma_draws, # ssvs-specific
                 omega_draws = omega_draws  # ssvs-specific
                 )
  return(answer)
}






#' Print summary for bayesian VAR model
#'
#' Print summary for bayesian VAR model
#'
#' @param bvar.model the list containing all results of bayesian VAR estimation
#' @export
#' @examples
#' bvar(Yraw)
bvar.summary <- function(bvar.model) {
  # Print some directions to the user
  message('Please find the means and variances of the VAR parameters in the vectors')
  message('ALPHA_mean and ALPHA_std for the VAR regression coefficients, and ')
  message('SIGMA_mean and SIGMA_std for the VAR covariance matrix. The predictive')
  message('mean and standard deviation are in Y_pred_mean and Y_pred_std, respectively.')
  message('The log Predictive Likelihood is given by variable log_PL. The true value')
  message('of y(t+h) is given in the variable true_value. For example the mean squared')
  message('forecast error can be obtained using the command')
  message('                MSFE = (Y_pred_mean - true_value).^2')
  message('If you are using the ssvs prior, you can get the averages of the restriction')
  message('indices $\\gamma$ and $\\omega$. These are in the variables gammas_mat and omega_mat') 
  
  
}

#' Draw irfs of bayesian VAR model
#'
#' Draw irfs of bayesian VAR model
#'
#' @param bvar.model the list containing all results of bayesian VAR estimation
#' @param qus the vector of quantiles for irfs
#' @return the ggplot2 graph of IRFs
#' @export
#' @examples
#' bvar(Yraw)
bvar.imp.plot <- function(bvar.model, qus = c(0.1, 0.5, 0.90)) {
  
  imp_responses <- apply(bvar.model$all_responses,c(2,3,4),quantile,probs=qus)
  
  imp_resp_melt <- melt(imp_responses)
  colnames(imp_resp_melt) <- c("probability","from","to","lag","impulse")
  
  
  p <- ggplot(data=imp_resp_melt,aes(x=lag,y=impulse)) +
    geom_line(aes(col=probability)) + facet_wrap(from~to,scales = "free")
  # print(p)
  return(p)
}

#' Histograms of forecasts
#'
#' This function draws histograms of forecasts
#'
#' @param bvar.model the list containing all results of bayesian VAR estimation
#' @return the ggplot graph of forecast histograms
#' @export
#' @examples
#' bvar(Yraw)
bvar.pred.plot <- function(bvar.model) {
  p <- ggplot(data=bvar.model$Y_pred_melt,aes(x=value)) + geom_histogram() + 
    facet_wrap(~Var2)
  # print(p)
  return(p)
}




