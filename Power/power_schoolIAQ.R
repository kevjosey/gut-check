
library(parallel)
library(swdpwr)
library(geepack)
library(sandwich)
library(dplyr)
options(dplyr.summarise.inform = FALSE)

# data dimensions
n.seq <- seq(24, 504, by = 12) # number of schools
n.iter <- 1000 # number of simulations

# parameters
p0 <- 1.6
p1 <- 0.75*p0
p0_fall <- 1.2

sig2.seq <- seq(0.01, 1, by = 0.01)

# calculate icc
icc.seq <- c(0, 0.05, 0.1, 0.15, 0.2)
icc.seq.all <- sapply(sig2.seq, function(z, ...) {
  iccCounts:::r_Pois(log(p1), log(sqrt(z)))
})

out_list <- mclapply(n.seq, function(n, ...) {

    school_id <- 1:n
    school_size <- sample(200:1000, size = n, replace = TRUE)
    area_id <- rep(1:4, times = n/4)
    treat_id <- rep(c(0,1), times = c(n/2, n/2)) # treatment
    
    dat_id <- data.frame(school = school_id, area = area_id, treat = treat_id)
    dat <- dat_id[rep(school_id, times = school_size),]
    dat$area <- factor(dat$area)
    
    beta <- c(log(p0), log(p1) - log(p0), rnorm((2*(nlevels(dat$area) - 1)))) # random effects for site modeled as fixed effects
    output <- data.frame()
    
    for (icc in icc.seq) {
    
      sig2 <- sig2.seq[which.min(abs(icc.seq.all - icc))]
    
      test_gee <- sapply(1:n.iter, function(i, ...) {
        
        print(i)
        
        X <- model.matrix(~ treat*area, data = dat)
        alpha <- rnorm(n, 0, sqrt(sig2))
        beta[3:(2*(nlevels(dat$area) - 1) + 2)] <- rnorm((2*(nlevels(dat$area) - 1)), 0, sd = 0.1)
        mu <- exp(alpha[dat$school] + c(X %*% beta))
        dat$y <- rpois(nrow(X), lambda = mu)
        
        dat_reduce <- dat %>% group_by(school, area, treat) %>% summarize(y = sum(y), size = n())
        dat_reduce$ybar <- with(dat_reduce, y/size)
        
        fit <- glm(ybar ~ treat, family = quasipoisson(link = "log"), data = dat_reduce, weight = size)
        
        est <- coef(fit)[2]
        se <- sqrt(vcovHC(fit)[2,2])
        cut <- est + qnorm(0.95)*se
        return(mean(cut < 0))
        
      })
    
      output <- rbind(output, data.frame(icc = round(icc, 3), pwr = mean(test_gee))) # power
    
  } 
  
  return(output) 
  
}, mc.cores = 30)

output <- cbind(n = rep(n.seq, each = length(icc.seq)), do.call(rbind, out_list))

write.csv(output, file = "~/Documents/breathe_power_results.csv")
