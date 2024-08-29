library(clusterPower)
library(swdpwr)
library(geepack)
library(dplyr)

# data dimensions
n <- 3000
m <- 500 # number of PCP Providers
l <- 5 # number of sites
n.iter <- 1000 # number of simulations

# data
pcp_id <- 1:m # PCP ID variable
site_id <- factor(sample(1:l, size = m, replace = T)) # Random Site Allocation
site <- factor(rep(site_id, each = n/m))
treat <- rep(c(0,1), times = c(2*n/3, n/3)) # treatment
pcp <- factor(rep(pcp_id, each = n/m))
covars <- data.frame(pcp = pcp, treat = treat, site = site)

# parameters
p0 <- 0.25
p1 <- 0.31
icc <- 0.0

beta <- c(qlogis(p0), qlogis(p1) - qlogis(p0),
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
