
library(clusterPower)
library(parallel)
library(swdpwr)
library(geepack)
library(sandwich)
library(dplyr)
options(dplyr.summarise.inform = FALSE)

# data dimensions
n.seq <- seq(240, 600, by = 20) # number of schools
n.iter <- 1000 # number of simulations

# parameters
p0 <- 0.2
p1 <- 0.75*p0
icc <- 0

sig2.seq <- seq(0.1, 2.3, by = 0.1)

# calculate icc
icc.seq <- sapply(sig2.seq, function(z, ...) {
  BinICC(link = "logit", meanresponse_start = p0, tau2 = z)$ICC
})

output <- data.frame()

out_list <- mclapply(icc.seq, function(icc, ...) {
  
  sig2 <- sig2.seq[which.min(abs(icc.seq - icc))]
  output <- data.frame()
  
  for (n in n.seq) {

    school_id <- 1:n
    school_size <- sample(200:1000, size = n)
    area_id <- rep(1:4, times = n/4)
    treat_id <- rep(c(0,1), times = c(n/2, n/2)) # treatment
    
    dat_id <- data.frame(school = school_id, area = area_id, treat = treat_id)
    dat <- dat_id[rep(school_id, times = school_size),]
    
    beta <- c(qlogis(p0), qlogis(p1) - qlogis(p0), rnorm(2)) # random effects for site modeled as fixed effects
    
    test_gee <- sapply(1:n.iter, function(i, ...) {
      
      print(i)
      
      X <- model.matrix(~ treat*area, data = dat)
      alpha <- rnorm(n, 0, sqrt(sig2))
      beta[3:4] <- rnorm(2)
      mu <- plogis(alpha[dat$school] + c(X %*% beta))
      dat$y <- rbinom(nrow(X), size = 1, prob = mu)
      
      dat_reduce <- dat %>% group_by(school, area, treat) %>% summarize(y = sum(y), size = n())
      dat_reduce$ybar <- with(dat_reduce, y/size)
      
      fit <- glm(ybar ~ treat*area, family = quasipoisson(link = "log"), data = dat_reduce, weight = size)
      
      est <- coef(fit)[2]
      se <- sqrt(vcov(fit)[2,2])
      cut <- est + qnorm(0.95)*se
      return(mean(cut < 0))
      
    })
    
    output <- rbind(output, data.frame(n = n, pwr = mean(test_gee))) # power
    
  } 
  
  return(output) 
  
}, mc.cores = 8)


output <- cbind(icc = rep(round(icc.seq, 3), each = length(n.seq)), do.call(rbind, out_list))
