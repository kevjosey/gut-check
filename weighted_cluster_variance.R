library(dplyr)

n <- 10000 # sample size
m <- 400 # clusters
n.iter <- 10000

mu <- vector(mode = "numeric", length = n.iter)
tau <- vector(mode = "numeric", length = n.iter)
ci <- matrix(NA, nrow = n.iter, ncol = 2)

for (i in 1:n.iter) {
  
  x <- rnorm(n, 0, 1)
  z <- rnorm(m, -1, 1)
  w <- sample(1:100, n, replace = T, prob = (100:1)/100) #offset
  id <- sample(1:m, n, replace = TRUE)
  zx <- z[id]
  alpha <- rnorm(m, -2, 1)

  y <- rbinom(n, size = w, exp(alpha[id])/(1 + exp(alpha[id])))
  ybar <- y/w
  data <- data.frame(id = id, ybar = ybar, m = m)
  mu[i] <- weighted.mean(exp(alpha[id])/(1 + exp(alpha[id])), w = w)
  
  tau[i] <- weighted.mean(ybar, w = w)
  dat <- data.frame(id = id, ybar = ybar, w = w)
  
  tau.dat1 <- dat %>% group_by(id) %>% 
    summarise(diff = weighted.mean(ybar, w = w) - tau[i],
              w = sum(w))
  tau.dat2 <- dat %>% group_by(id) %>% 
    mutate(diff = ybar - weighted.mean(ybar, w = w))
  
  # tau_var <- (sum(w^2*(ybar - tau[i])^2)/sum(w))/sum(w)
  
  tau_var <- with(tau.dat1, sum(w*(diff)^2)/sum(w))/(nrow(tau.dat1) - 1) +
    with(tau.dat2, sum(w*(diff)^2)/sum(w))/(nrow(tau.dat2) - nrow(tau.dat1))
  
  ci[i,] <- c(tau[i] - 1.96*sqrt(tau_var), 
              tau[i] + 1.96*sqrt(tau_var))  
  
}

mean(ci[,1] < mean(mu) & ci[,2] > mean(mu))
