library(clusterPower)
library(swdpwr)
library(interactionR)
library(geepack)
library(dplyr)

# data dimensions
n <- 120 # number of subjects
m <- 5 # number of repeated measurements (unequal spacing)
n.iter <- 1000 # number of simulations

# data
pcp_id <- 1:m # PCP ID variable
site_id <- factor(sample(1:l, size = m, replace = T)) # Random Site Allocation
therapy <- rep(c(0,1), times = c(n/2, n/2)) # treatment


pcp0 <- sample(pcp_id[1:(m/2)], size = n/2, replace = T)
pcp1 <- sample(pcp_id[(m/2 + 1):m], size = n/2, replace = T)
pcp <- c(pcp0, pcp1)
pcp <- pcp[order(pcp)]
site <- factor(rep(site_id, times = table(pcp)))

covars <- data.frame(pcp = pcp, treat = treat, site = site)

# parameters
p0 <- 0.25
p1 <- 0.31
icc <- 0.0

beta <- c(log(p0/(1 - p0)), 0,
          rnorm(l - 1, 0, 0.5)) # random effects for site modeled as fixed effects

sig2_seq <- seq(0.1, 2, by = 0.05)

# calculate icc
icc_seq <- sapply(sig2_seq, function(z, ...) {
  BinICC(link = "logit", meanresponse_start = p0, tau2 = z)$ICC
})

sig2 <- sig2_seq[which.min(abs(icc_seq - icc))]

# GEE simulation
test_gee <- vector(mode = "numeric", length = n.iter)

for (i in 1:n.iter) {
  
  print(i)
  
  X <- model.matrix(~ treat + site, data = covars)
  alpha <- rnorm(m, 0, sqrt(sig2))
  beta[3:(l+1)] <- rnorm(l - 1, 0, 0.5)
  mu <- plogis(alpha[pcp] + c(X %*% beta))
  y <- rbinom(nrow(X), size = 1, prob = mu)
  
  fit <- geeglm(y ~ treat + site, family = binomial(link = "logit"), id = pcp, 
                data = data.frame(y = y, treat = treat, site = site, pcp = pcp),
                corstr = "exchangeable")
  
  cc <- coef(summary(fit))
  cut <- with(as.data.frame(cc), Estimate[2] + qnorm(0.95)*Std.err[2])
  test_gee[i] <- as.numeric(cut < (qlogis(p1) - qlogis(p0)))
  
}

mean(test_gee) # power

## with longpower.R
