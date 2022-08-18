####################################################
## Purpose: Generate longitudinal data determined ##
##          by covariance pattern models          ##
## By: Kevin Josey                                ##
####################################################

# Dependencies ------------------------------------------------------------

library(nlme) # needed for gls()
set.seed(06261992) # for replication purposes


# Functions for generating different covariance patterns ------------------

# Exchangeable covariance pattern
xch <- function(sig2, rho, p){
  
  if(length(sig2) != 1)
    stop("length(sig2) != 1")
  
  if(length(rho) != 1)
    stop("length(rho) != 1")
  
  R <- matrix(rho, p, p)
  diag(R) <- 1
  D <- diag(sqrt(sig2), p, p)
  V <- D %*% R %*% D
  return(V)
  
}

# AR(1) covariance pattern
ar1 <- function(sig2, phi, p){
  
  if(length(sig2) != 1)
    stop("length(sig2) != 1")
  
  if(length(phi) != 1)
    stop("length(phi) != 1")
  
  times <- 1:p
  H <- abs(outer(times, times, "-"))
  V <- sig2 * phi^H
  return(V)
  
}

# Unstructured covariance pattern
un <- function(sig2, R){
  
  if (!isSymmetric(R))
    stop("R must be a square/symmetric matrix")
  
  if( any(eigen(R)$values <= 0) )
    stop("R must be positive definite")
  
  if (any(diag(R) != 1))
    stop("R must have 1 along the diagonals")
  
  if (any(abs(R) > 1))
    stop("off diagonal entrie of R must be between (-1,1)")
  
  if (length(sig2) == 1)
    sig2 <- rep(sig2, times = nrow(R))
  else if (length(sig2) != nrow(R))
    stop("length(sig2) != nrow(R)")
  
  if (any(sig2 <= 0))
    stop("sig2 must be strictly positive")
  
  p <- nrow(R)
  D <- diag(sqrt(sig2), p, p)
  V <- D %*% R %*% D
  
  return(V)
  
}


# Simulate covariance pattern models --------------------------------------

# data dimensions
n1 <- 100 # sample size for group 1
n2 <- 100 # sample size for group 2
n <- n1 + n2 # total sample size
p <- 4 # number of repeated measurements

# covariance matrices
V_xch <- xch(sig2 = 4, rho = 0.3, p = p) # exchangeable
V_ar1 <- ar1(sig2 = 4, phi = 0.3, p = p) # AR(1)
R <- matrix(c(1, 0.1, -0.2, -0.3, 
              0.1, 1, 0.4, 0.8,
              -0.2, 0.4, 1, 0.6, 
              -0.3, 0.8, 0.6, 1),
            nrow = 4, ncol = 4) 
V_un <- un(sig2 = c(1,2,3,4), R = R) # unstructured-ish

# design matrix
grp <- factor( rep(1:2, times = c(n1*p,n2*p)) ) # identify treatment assignment
time <- rep(0:(p-1), times = n) # we can also factor this variable
X <- model.matrix(~ grp*time) # reference cell coding
id <- rep(1:n, each = p) # id indexing clusters

# fixed effect coefs
beta <- c(1, -1, 2, 0.5) 
# c(intercept, grp2, time, grp2*time)

E_xch <- E_ar1 <- E_un <- matrix(NA, nrow = n, ncol = p) # error matrix
# vectorizing the error matrix only works when the design consists of
# an equal number of planned measurements per independent sampling unit

# generate errors for each independent sampling unit (by rows)
for (i in 1:n) {
  
  E_xch[i,] <- t(chol(V_xch)) %*% rnorm(p, mean = 0, sd = 1)
  E_ar1[i,] <- t(chol(V_ar1)) %*% rnorm(p, mean = 0, sd = 1)
  E_un[i,] <- t(chol(V_un)) %*% rnorm(p, mean = 0, sd = 1)
  # want L not U in LU Decomposition
  # we can also use mvrnorm() from MASS package (this uses SVD)
  
}

# vectorize error matrix by rows (independent sampling units)
e_xch <- as.vector(t(E_xch)) 
e_ar1 <- as.vector(t(E_ar1))  
e_un <- as.vector(t(E_un)) 

# generate responses
y_xch <- c( X %*% beta + e_xch )
y_ar1 <- c( X %*% beta + e_ar1 )
y_un <- c( X %*% beta + e_un )

# fit mixed models
xch_fit <- gls(y_xch ~ grp*time, correlation = corCompSymm(form = ~1|id)) # exchangeable model
ar1_fit <- gls(y_ar1 ~ grp*time, correlation = corAR1(form = ~1|id)) # AR(1) model
un_fit <- gls(y_un ~ grp*time, correlation = corSymm(form = ~1|id),
              weights = varIdent(form = ~1|time)) # unstructured model


# Feasible Generalized Least Squares (without nlme) -----------------------

y_un <- y_un[order(id, time)] # make sure y is ordered appropriately

# status variables
tol <- 1e-3 # tolerance for convergence
max_iter <- 100 # maximum number of iterations allowed
iter <- 0 # initialize iteration count
converged <- FALSE # convergence indicator

b_hat <- c( solve(t(X) %*% X) %*% t(X) %*% y_un ) # initialize beta with ordinary least squares
Sigma_hat <- diag(1, nrow = p, ncol = p)

while(!converged) {
  
  iter <- iter + 1 # track how many iterations algo takes
  b_0 <- b_hat # b.hat will be updated
  Sigma_0 <- Sigma_hat
  rss <- matrix(0, p, p) # initialize RSS matrix
  
  # find RSS
  for(j in 1:n){
    
    y_tmp <- y_un[id == j]
    X_tmp <- X[id == j,]
    rss <- rss + (y_tmp - X_tmp %*% b_0) %*% t(y_tmp - X_tmp %*% b_0)
    
  }
  
  Sigma_hat <- rss/(n - ncol(X)) # RMSE matrix
  Sigma_inv <- solve(Sigma_hat) # invert this now (more effecient computationally)
  Omega <- kronecker(diag(1, n, n), Sigma_inv) # diagonalize inverse of the RMSE
  
  b_hat <- c( solve(t(X) %*% Omega %*% X) %*% t(X) %*% Omega %*% y_un ) # gls estimate
  
  if ( max(abs(b_hat - b_0)) < tol & max(abs(Sigma_hat - Sigma_0)) < tol ) # check convergence
  
  if (iter == max_iter) # prevents algorithm from running forever
    break
  
}

converged # did the model converge?
Sigma_hat # estimate of Sigma
getVarCov(un_fit)
b_hat # estimate of beta
coef(un_fit)
