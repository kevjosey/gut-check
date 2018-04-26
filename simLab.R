###############################################
## PURPOSE: Demonstration of Simulations for ##
##          Linear Models                    ##
## BY:      Kevin Josey                      ##
###############################################

# Dependencies ------------------------------------------------------------

library(nlme) # for gls()
library(lme4) # for lmer()
set.seed(06041987)


# Generating Random Variables ---------------------------------------------

n <- 10000

## Continuous Random Variables

# Exponential
exp_dist <- rexp(n, rate = 1/3)
hist(exp_dist, main = "Exponential Distribution, lambda = 3")

# Gamma
par(mfrow = c(1,3))
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
table(pois_dist)

# Bernoulli
bern_dist <- rbinom(n, size = 1, prob = 1/3)
table(bern_dist)


# Probability Integral Transformation -------------------------------------

U <- runif(n)

# Exponential
lambda <- 3
Y <- -lambda*log(1-U)
mean(Y) # lambda
var(Y) # lambda^2

# f(x) = 1/2 cos(x); F(X) = 1/2 sin(X) + 1/2; -pi/2 < X < pi/2

X <- asin(2*U - 1)
mean(X) # mu = 0
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
e <- rnorm(nrow(X), sd = 1.5) # mse = 2.25
mu <- c(1,1,2) # parameters

y <- X %*% mu + e

# Fitting the Model
one.way.fit <- lm(y ~ 0 + grp) # note that we remove the intercept
summary(one.way.fit)

# ANOVA by Hand
C <- rbind(c(1, -1, 0),
           c(1, 0, -1)) # contrast matrix
mu.hat <- solve(t(X) %*% X) %*% t(X) %*% y # OLS estimator
rmse <- t(y - X %*% mu.hat) %*% (y - X %*% mu.hat) / (n - 3) # residual mse

I <- diag(1, nrow = n, ncol = n)
proj.mat <- X %*% solve(t(X) %*% X) %*% t(X) # projection matrix
rmse.alt <- t(y) %*% (I - proj.mat) %*% y / (n - 3) # alternative rmse

mu.0 <- c(1,1,1) # null hypothesis
mu.var <- solve(t(X) %*% X) * as.vector(rmse) # variance of mu.hat

wald.stat <- t(C %*% (mu.hat - mu.0)) %*% solve(C %*% mu.var %*% t(C)) %*% (C %*% (mu.hat - mu.0)) / 2 # Wald statistic

pval <- pf(wald.stat, 2, 57, lower.tail = FALSE) # p-value

anova(lm(y ~ grp)) # check

## Two-Way ANOVA w/ interaction

rm(list = ls())

# Design Matrix (factor 1: 3 lvls, factor 2: 2 levels, balanced design)
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
two.way.fit <- lm(y ~ a1*a2)
anova(two.way.fit)
summary(two.way.fit)

## ANCOVA

# Design Matrix
n1 <- 50
n2 <- 50

grp <- factor( rep(c(1,2), times = c(n1,n2)) )
cv <- rnorm(n1 + n2, mean = 2, sd = 1)
X <- model.matrix(~ grp + cv) # reference cell coding

# Response
e <- rnorm(nrow(X), sd = 2) # mse = 4
beta <- c(10, 0.8, 0.1) 

y <- X %*% beta + e

# Fitting the Model
ancova.fit <- lm(y ~ grp + cv)
summary(ancova.fit)


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
  # want L not U in LU Decomposition; use unit variance
  
}

e <- as.vector(t(E)) # vectorize by rows of the error matrix 
y <- X %*% beta + e

# Fit Mixed Model
cs.fit <- gls(y ~ grp*time, correlation = corCompSymm(form = ~1|id)) # generalized least squares (nlme)
summary(cs.fit) # ICC = 1/(1+2) = 0.3333 = rho
anova(cs.fit)

# Difference in Differences Hypothesis by Hand ( test interactions )

L <- rbind(c(0,0,0,0,1,0),
           c(0,0,0,0,0,1)) # contrast matrix

nu.e = n - 3 # residual degrees of freedom
nu.1 = 2 # numerator degrees of freedom
nu.2 = nu.e - 1 + 1 # denominator defrees of freedom

beta.hat <- coef(cs.fit) # estimate
theta.hat <- L %*% beta.hat # contrast estimates
theta.0 <- as.matrix( c(0, 0) ) # null hypothesis

# if you can find a better way of extract the covariance estimate, please let me know
V.tmp <- matrix(coef(cs.fit$modelStruct$corStruct, unconstrained = FALSE), nrow = p, ncol = p) 
diag(V.tmp) <- 1
V.hat <- V.tmp*cs.fit$sigma^2 # within subject variance estimate
Sig.hat = kronecker(diag(1, nrow = n, ncol = n), V.hat) # expand V

wald.var <- solve(L %*% solve(t(X) %*% solve(Sig.hat) %*% X) %*% t(L)) # variance of theta.hat
wald.stat <- t(theta.hat - theta.0) %*% wald.var %*% (theta.hat - theta.0) / nu.1 # wald statistic

pval <- pf(wald.stat, nu.1, nu.2, lower.tail = FALSE) # p-value


# Linear Mixed Effects Model ----------------------------------------------

# Using everything from before but fitting a random intercept model 
# instead of a "compound symmetric pattern" model

d <- rep(rnorm(n, 0, sqrt(1)), each = p)
e <- rnorm(n*p, 0, sqrt(2)) 
# Same ICC as the covariance pattern model,
# errors are independent,
# between cluster heterogeneity is expressed in the random effect

y_n <- X %*% beta + d + e # Z_{ij} = 1 for all i,j

# Fit Model
rand.int.fit <- lmer(y_n ~ grp*time + (1|id)) # using lme4 package
summary(rand.int.fit) # Same distribution as cs.fit response variable


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

# One-Way ANOVA

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
  n.reject <- 0 # initialize rejection count

  for (i in 1:iter) {
  
    e <- rnorm(nrow(X), sd = sd) # residual error
    y <- X %*% beta + e
  
    fit <- lm(y ~ grp) # fitted model
    a <- anova(fit)
  
    n.reject <- n.reject + ( a$`F value`[1] > qf(0.95, 2, 57) ) 
    # ndf = rank(C) = 2, ddf = N - rank(X) = 60 - 3 = 57
  
  }
  
  out[j] <- n.reject/iter
    
}

pow <- cbind(scen, out)
plot(out ~ scen, data = pow, type = "l", lwd = 2,
     ylab = "Power", xlab = "Concentrated Effect Size", 
     main = "Power with Respect to Concentrated Effects")
abline(h = 0.8, col = "red", lty = 3, lwd = 2)
segments(x0 = 2.5, y0 = -0.1, y1 = 0.8, lty = 3, col = "red", lwd = 2)

## Sample Size

# Consider an equal allocation study (group sample sizes are equal)
# Fix: effect parameters, nuisance parameters, covariates
# Same procedure as sample size

rm(list = ls())

# One-Way ANOVA

iter <- 10000 # no. of iterations
sd <- 3 # mse = 9
try.n <- seq(5, 50, by = 5) # test different sample sizes
beta <- c(1, 0, 2) # note this is our "effect size"

out <- vector(mode = "numeric", length = length(try.n)) # initialize output

for(j in 1:length(try.n)) { # nested for loops are bad, see if you can't improve this with the apply class of functions
  
  n <- try.n[j]
  grp <- factor( rep(1:3, each = n) )
  X <- model.matrix(~ grp) # Reference Cell Coding
  n.reject <- 0 # initialize rejection count
  
  for (i in 1:iter) {
    
    e <- rnorm(nrow(X), sd = sd) 
    y <- X %*% beta + e
    
    fit <- lm(y ~ grp) # fitted model
    a <- anova(fit)
    
    n.reject <- n.reject + ( a$`F value`[1] > qf(0.95, 2, 57) ) 
    # ndf = rank(C) = 2, ddf = N - rank(X) = 60 - 3 = 57
    
  }
  
  out[j] <- n.reject/iter
  
}

samp.size <- cbind(try.n, out)
