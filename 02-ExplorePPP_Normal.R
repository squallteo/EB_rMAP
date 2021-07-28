rm(list=ls())

package2load <- c("RBesT", "tidyverse", "doParallel", "ggplot2")
lapply(package2load, require, character.only = TRUE)

set.seed(712)


ppp_cut <- 0.1 #c for ppp

#derive and approximate MAP prior
dt <- crohn
sigma <- 88
dt$se_yh <- sigma/sqrt(dt$n)
#meta analysis of historical data
# with(dt, meta::metamean(n = n, mean = y, sd = rep(sigma,length(study))))
#point est is -50. Decide to consider true current mean from -60 to -40
map_mcmc <- gMAP(cbind(y, se_yh) ~ 1 | study, 
                  #weight = n,
                  data=dt,
                  family=gaussian,
                  beta.prior=cbind(-50, 88),
                  tau.dist="HalfNormal",tau.prior=cbind(0,11))
#approximate the MAP
map_hat <- automixfit(map_mcmc)


# map_mcmc1 <- gMAP(cbind(y, se_yh) ~ 1 | study, 
#                  weight = n,
#                  data=dt,
#                  family=gaussian,
#                  beta.prior=cbind(0, 100),
#                  tau.dist="HalfNormal",tau.prior=cbind(0,44))
# #approximate the MAP
# map_hat1 <- automixfit(map_mcmc1)



sigma(map_hat) <- 88
n_c <- 30
y_cvec <- seq(-55, -40, 0.1)

for(y in 1:length(y_cvec)){
  y_c <- y_cvec[y]
  #calculate ppp and pick the optimal w_V
  w <- seq(0, 1, 0.01)
  ppp <- rep(NA, length(w))
  for(i in 1:length(w)){
    rmap <- robustify(map_hat, weight=w[i], mean=-50, sigma=50)
    rmap_pred <- preddist(rmap, n=n_c)
    p_lower <- pmix(rmap_pred, y_c)
    ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
  }
  # View(tibble(w, ppp))
  pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(ppp==max(ppp) & pass)
  w_eb <- ifelse(nrow(pppdt) == 0, 1, max(pppdt$w))
  ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))
  tt <- tibble(y_c, w_eb, ppp_eb)
  if(y==1){wdt <- tt}
  else {wdt <- rbind(wdt, tt)}
} 

ggplot(wdt, aes(x=y_c, y = w_eb)) + geom_line()

