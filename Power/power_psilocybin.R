library(geepack)
library(swdpwr)
library(MASS)
library(tidyr)

## Functions for covariance matrices

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
  
  if (length(p) == 1)
    times <- 1:p
  else
    times <- p
  
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

## parameter values

alpha <- 0.05
n_seq <- seq(10, 100, by = 10)
days_seq <- c(0,1,7,28,56,84)
trt_seq <- c("Low + Ctrl", "High + Ctrl", "Low + Trt", "High + Trt")
rho <- 0.7
n.iter <- 100
result <- data.frame()

# start simulation

for (i in 1:length(n_seq)) {
  
  n <- n_seq[i]
  
  trt <- rep(trt_seq, each = n)
  R <- ar1(sig2 = 1, phi = rho, p = 0:5)
  Sigma <- un(sig2 = c(6*(1 + 0.25*(0:(length(days_seq) - 1))))^2, R = R)
  
  test.stat <- vector(mode = "numeric", length = n.iter)
  
  for (j in 1:n.iter) {
  
    y_wide <- matrix(35 - rep(c(0,0.5,1), times = c(4*n, 4*n, 16*n))*
                       (12*as.numeric(trt == "High + Trt" | trt == "High + Ctrl") +
                          6*as.numeric(trt == "Low + Trt" | trt == "Low + Ctrl")), nrow = 4*n) + 
      mvrnorm(4*n, rep(0, 6), Sigma)
    
    colnames(y_wide) <- 0:5
    
    data <- data.frame(y_wide) %>% 
      pivot_longer(
        cols = 1:6, 
        names_to = "time",
        values_to = "y",
      )
    
    data$time <- factor(data$time, label = days_seq)
    id <- factor(rep(1:nrow(y_wide), each = 6))
    trt_long <- factor(rep(trt, each = 6), 
                       levels = c("Low + Ctrl", "High + Ctrl",
                                  "Low + Trt", "High + Trt")) 
    
    data <- data.frame(data, id = id, trt = trt_long)
    model <- geeglm(y ~ trt*time, data = data, id = id, waves = time, corstr = "unstructured")
    
    # 8-week Difference Estimator
    # contr <- rbind(c(0, 0, -1, 1, rep(0, 15), -1, 1, rep(0, 3)),
    #                   c(0, -1, 0, 1, rep(0, 14), -1, 0, 1, rep(0, 3)))
    
    # DiD estimator
    contr <- rbind(c(rep(0, 19), -1, 1, rep(0, 3)),
                   c(rep(0, 18), -1, 0, 1, rep(0, 3)))
    
    vcov_gee = model$geese$vbeta
    
    contrast_est = coef(model) %*% t(contr)
    contrast_se = sqrt(contr %*% vcov_gee %*% t(contr))
    
    output = data.frame(Estimate = contrast_est[1,1],
                        SE = diag(contrast_se)[1]) %>% 
      mutate(LCI = Estimate - qnorm(1 - alpha/(2))*SE,
             UCI = Estimate + qnorm(1 - alpha/(2))*SE) # Bonferroni baked in
             
    test.stat[j] <- as.numeric(all(output$UCI < 0))
    
  }
  
  result <- rbind(result, data.frame(n = n, pwr = mean(test.stat, na.rm = T)))
  
}
