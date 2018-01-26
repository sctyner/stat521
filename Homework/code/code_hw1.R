# write function to calculate the HT variance estimator and the HT variance of the simple sample 
library(tidyverse)
# data, inclusion probs
y <- c(5,3,7)
pi_y <- c(.8,.7,.6)

# samples, pr(A)s
samps <- data_frame(A = list(c(1,2), c(1,3), c(2,3), c(1,2,3)), 
                    PrA = c(.4, .3, .2, .1))
# pairs, 2-way inclusion probs
pi_kl <- data.frame(k = c(1,1,2,1,2,3, 2,3,3), l = c(2,3,3,1,2,3, 1,1,2), pi_kl = c(.5, .4, .3, pi_y,.5, .4, .3))

# matrix form 
my_PI <- pi_kl %>% spread(l, pi_kl) %>% select(-k) %>% as.matrix()


# Compute V(t_ht) for U = {1,2,3}
Delta <- matrix(0, nrow = 3, ncol = 3)
for(i in 1:3){
  for(j in 1:3){
      Delta[i,j] <- my_PI[i,j] - (my_PI[i,i] * my_PI[j,j])
    }
  }
Check_y <- (y/pi_y) %*% t((y/pi_y))
v_ht <- sum(colSums(Delta * Check_y))
v_ht

# Compute estimator for each 
v_ht_hat <- function(samp, dat = y, pikl = my_PI) {
  n <- length(samp)
  Check_delta <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      k <- samp[i]
      l <- samp[j]
      Check_delta[i,j] <- (pikl[k,l] - pikl[k,k]* pikl[l,l]) / pikl[k,l]
    }
  }
  y_s <- dat[samp]
  piy_s <- diag(pikl)[samp] 
  check_y <- (y_s/piy_s) %*% t((y_s/piy_s))
  est_v_ht <- sum(colSums(Check_delta * check_y))
  return(est_v_ht)
}

samps2 <- samps %>% mutate(V_ht_est = map_dbl(A, v_ht_hat), 
                           EV = PrA * V_ht_est)
samps2$EV %>% sum
# hm not getting it. try this package
#install.packages("samplingbook")
library(samplingbook)
htestimate(y = y, N = 3, PI = my_PI, method = "ht")


