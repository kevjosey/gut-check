###############################################
## PURPOSE: Simulations for Linear Models    ##
## BY:      Kevin Josey                      ##
###############################################

# Dependencies ------------------------------------------------------------

library(nlme) # for gls()
set.seed(06041989) # for replication purposes


# Generating Random Variables ---------------------------------------------

n <- 10000

## Continuous Random Variables

# Exponential
exp_dist <- rexp(n, rate = 1/3)
hist(exp_dist, main = "Exponential Distribution, lambda = 3")

par(mfrow = c(1,3)) # 3 panel plots

# Gamma
gamma_dist <- rgamma(n, shape = 20/2, rate = 1/2)
hist(gamma_dist, main = "Gamma with Rate Parameter")
gamma_dist_alt <- rgamma(n, shape = 20/2, scale = 2)
hist(gamma_dist_alt, main = "Gamma with Scale Parameter")

# Chi-square: special case of gamma
chisq_dist <- rchisq(n, df = 20)
hist(chisq_dist, main = "Chi-Square Equivalent")

# Normal
par(mfrow = c(1,1))
norm_dist <- rnorm(n, mean = 0, sd = 4)
hist(norm_dist, main = "Normal Distribution")

## Discrete Random Variables

# Poisson
pois_dist <- rpois(n, lambda = 4)
hist(pois_dist, main = "Poisson Distribution")

# Bernoulli
bern_dist <- rbinom(n, size = 1, prob = 1/3)
table(bern_dist)


# Probability Integral Transformation -------------------------------------

U <- runif(n)

# Exponential
lambda <- 3
Y <- -lambda*log(1-U)
mean(Y) # mean = lambda
var(Y) # var = lambda^2

# f(x) = 1/2 cos(x); F(X) = 1/2 sin(X) + 1/2; -pi/2 < X < pi/2
X <- asin(2*U - 1)
mean(X) # mean = 0
var(X) # var = 1/4*(pi^2-8)


# General Linear Model ----------------------------------------------------

## One-way ANOVA

rm(list = ls())

# Design Matrix
n1 <- 15 # sample size for group 1
n2 <- 20 # sample size for group 2
n3 <- 25 # sample size for group 3
n <- n1 + n2 + n3

grp <- factor( rep(1:3, times = c(n1,n2,n3)) )
X <- model.matrix(~ 0 + grp) 
# cell means coding (only time I will use this)

# Random Reponses
e <- rnorm(n, sd = 1.5) # mse = 2.25
mu <- c(1,1,2) # parameters

y <- X %*% mu + e

# Fitting the Model
one_way_fit <- lm(y ~ 0 + grp) # note that we remove the intercept
summary(one_way_fit)

# ANOVA by Hand
C <- rbind(c(1, -1, 0),
           c(1, 0, -1)) # contrast matrix
mu_hat <- solve(t(X) %*% X) %*% t(X) %*% y # OLS estimator
rmse <- t(y - X %*% mu_hat) %*% (y - X %*% mu_hat) / (n - 3) # residual MSE

I <- diag(1, nrow = n, ncol = n)
proj_mat <- X %*% solve(t(X) %*% X) %*% t(X) # projection matrix
rmse_alt <- t(y) %*% (I - proj_mat) %*% y / (n - 3) # alternative RSME

mu_0 <- c(1,1,1) # null hypothesis
mu_var <- solve(t(X) %*% X) * as.vector(rmse) # variance of mu.hat

wald_stat <- t(C %*% (mu_hat - mu_0)) %*% solve(C %*% mu_var %*% t(C)) %*% (C %*% (mu_hat - mu_0)) / 2 # Wald statistic

pval <- pf(wald_stat, 2, 57, lower.tail = FALSE) # p-value

anova(lm(y ~ grp)) # check

## Two-Way ANOVA w/ interaction

rm(list = ls())

# Design Matrix (factor 1: 3 lvls, factor 2: 2 levels; balanced design)
n1 <- 30
n2 <- 30
n3 <- 30

m1 <- m2 <- 45

a1 <- factor( rep(1:3, times = c(n1,n2,n3)) )
a2 <- factor( rep(1:2, times = (m1 + m2)/2) )

X <- model.matrix(~ a1*a2) # reference cell coding, full rank matrix

# Response
e <- rnorm(nrow(X), sd = 1.5) # MSE = 2.25; independence
beta <- c(10, 1, -1, 2, 0.5, -0.5) # c(int, factor1.2, factor1.3, factor2.2, factor1.2:factor2.2, factor1.3:factor2.2)

y <- X %*% beta + e

# Fitting the Model
two_way_fit <- lm(y ~ a1*a2)
anova(two_way_fit)
summary(two_way_fit)


# Repeated Measures ANOVA ----------------------------------------------

rm(list = ls())

p <- 2 # repeated measurements

# Covariance Matrix
R <- matrix(1, nrow = p, ncol = p)
G <- diag(2, nrow = p, ncol = p)

V <- R + G # W/I subject covariance

# Design Matrix
n1 <- 100 # sample size for group 1
n2 <- 100 # sample size for group 2
n3 <- 100 # sample size for group 3
n <- n1 + n2 + n3

grp <- factor( rep(1:3, times = c(n1*p,n2*p,n3*p)) )
time <- factor( rep(1:p, times = n) )
X <- model.matrix(~ grp*time) # reference cell coding
id <- rep(1:n, each = p) # id indexing clusters

# Responses
beta <- c(10, 1, 2, 0.2, -0.1, 0.1) 
# c(int, grp2, grp3, time2, grp2*time2, grp3*time2)

E <- matrix(NA, nrow = n, ncol = p)
# vectorizing the error matrix only works when each
# independent sampling unit has an equal number of planned 
# repeated measurements.

for (i in 1:n) {
  
  E[i,] <- t(chol(V)) %*% rnorm(p, mean = 0, sd = 1)
  # want L not U in LU Decomposition
  
}

e <- as.vector(t(E)) # vectorize by rows of the error matrix 
y <- X %*% beta + e

# Fit Mixed Model
cs_fit <- gls(y ~ grp*time, correlation = corCompSymm(form = ~1|id)) # generalized least squares (nlme)
summary(cs_fit) # ICC = 1/(1+2) = 0.3333 = rho
anova(cs_fit)

# Difference in Differences Hypothesis by Hand ( test interactions )

L <- rbind(c(0,0,0,0,1,0),
           c(0,0,0,0,0,1)) # contrast matrix

nu_1 <- 2 # numerator degrees of freedom
nu_2 <- n - 3 # residual degrees of freedom

beta_hat <- coef(cs_fit) # estimate
theta_hat <- L %*% beta_hat # contrast estimates
theta_0 <- as.matrix( c(0, 0) ) # null hypothesis

# if you can find a better way of extract the covariance estimate, please let me know
V_tmp <- matrix(coef(cs_fit$modelStruct$corStruct, unconstrained = FALSE), nrow = p, ncol = p) 
diag(V_tmp) <- 1
V_hat <- V_tmp*cs_fit$sigma^2 # within subject variance estimate
Sig_hat <- kronecker(diag(1, nrow = n, ncol = n), V_hat) # expand V

wald_var <- solve(L %*% solve(t(X) %*% solve(Sig_hat) %*% X) %*% t(L)) # variance of theta_hat
wald_stat <- t(theta_hat - theta_0) %*% wald_var %*% (theta_hat - theta_0) / nu_1 # wald statistic

pval <- pf(wald_stat, nu_1, nu_2, lower.tail = FALSE) # p-value


# Power and Sample Size ---------------------------------------------------

## Power

# Monte Carlo Simulation for Power
# Fix: sample size, nuisance parameters, covariates
# For each true beta (effect size more or less):
# (1) Generate 10,000 samples. 
# (2) Run the desired hypothesis test for each sample.
# (3) Save whether you rejected the null or not.
# (4) Proportion of samples that reject is the empirical power

rm(list = ls())

iter <- 10000 # no. of iterations
n1 <- 15 # sample size group 1
n2 <- 20 # sample size group 2
n3 <- 25 # sample size group 3

grp <- factor( rep(1:3, times = c(n1,n2,n3)) )
X <- model.matrix(~ grp) # Reference Cell Coding

sd <- 3 # MSE = 9, independent
scen <- seq(0, 3, by = 0.5)

out <- vector(mode = "numeric", length = length(scen)) # initialize output

for(j in 1:length(scen)) { # nested for loops are bad, see if you can't improve this with the apply class of functions
  
  beta <- c(1, 0, scen[j]) # concentrated effects
  n_reject <- 0 # initialize rejection count

  for (i in 1:iter) {
  
    e <- rnorm(nrow(X), sd = sd) # residual error
    y <- X %*% beta + e
  
    fit <- lm(y ~ grp) # fitted model
    a <- anova(fit)
  
    n_reject <- n_reject + ( a$`F value`[1] > qf(0.95, 2, 57) ) 
    # ndf = rank(C) = 2, ddf = N - rank(X) = 60 - 3 = 57
  
  }
  
  out[j] <- n_reject/iter
    
}

pow <- cbind(scen, out)
plot(out ~ scen, data = pow, type = "l", lwd = 2,
     ylab = "Power", xlab = "Concentrated Effect Size", 
     main = "Power with Respect to Concentrated Effects")
abline(h = 0.8, col = "red", lty = 3, lwd = 2)
segments(x0 = 2.5, y0 = -0.1, y1 = 0.8, lty = 3, col = "red", lwd = 2)

## Sample Size

# Consider an equal allocation study (group sample sizes are equal)
# Fix coefficients and nuisance parameters
# Similar procedure as sample size with a few modifications

rm(list = ls())

iter <- 10000 # no. of iterations
sd <- 3 # mse = 9
try_n <- seq(5, 50, by = 5) # test different sample sizes
beta <- c(1, 0, 2) # note this is our effect size

out <- vector(mode = "numeric", length = length(try_n)) # initialize output

for(j in 1:length(try_n)) { # nested for loops are bad, see if you can't improve this with the apply class of functions
  
  n <- try_n[j]
  grp <- factor( rep(1:3, each = n) )
  X <- model.matrix(~ grp) # Reference Cell Coding
  n_reject <- 0 # initialize rejection count
  
  for (i in 1:iter) {
    
    e <- rnorm(nrow(X), sd = sd) 
    y <- X %*% beta + e
    
    fit <- lm(y ~ grp) # fitted model
    a <- anova(fit)
    
    n_reject <- n_reject + ( a$`F value`[1] > qf(0.95, 2, 57) ) 
    # ndf = rank(C) = 2, ddf = N - rank(X) = 60 - 3 = 57
    
  }
  
  out[j] <- n_reject/iter
  
}

samp_size <- cbind(try_n, out)
