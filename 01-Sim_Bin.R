rm(list=ls())

package2load <- c("RBesT", "tidyverse", "doParallel", "ggplot2")
lapply(package2load, require, character.only = TRUE)

set.seed(712)

#meta analysis of historical data
# with(dt, meta::metaprop(event = r, n = n, method = "Inverse"))
#point est is 0.25 (0.202, 0.312). Decide to consider true current rate from 0.2 to 0.32

#current control
n_c <- 50
pvec_c <- seq(0.20, 0.32, 0.01)
w_vec <- c(0.2, 0.5, 0.8) #w_v in rMAP
a_c <- 1; b_c <- 1 #vague beta prior for p(c)
#EB rMAP parameters
ppp_cut <- 0.1
#current treatment
n_t <- 100
effsize <- 0.2
prior_t <- mixbeta(c(1, 1, 1), param = "ab")

#decision rule to claim trial success: pr(p_t - p_c > Qcut) > Pcut
Qcut <- 0
Pcut <- 0.95
# success_rule = decision2S(pc = 0.95, qc = 0, lower.tail = F, link = "identity")
npostdist <- 20000
nsim <- 5000
#################################
#################################
#derive and approximate MAP prior
map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | study,
                 data=AS,
                 tau.dist="HalfNormal",
                 tau.prior=1,
                 beta.prior=2,
                 family=binomial)
map_hat <- automixfit(map_mcmc)

#################################
#################################
#for binary outcome, enumerate all possible y_c and calculate w_eb at each y_c
for(y in 0:n_c){
  y_c <- y
  #calculate ppp and pick the optimal w_V
  w <- seq(0, 1, 0.01)
  ppp <- rep(NA, length(w))
  for(i in 1:length(w)){
    rmap <- robustify(map_hat, weight=w[i], mean = a_c/(a_c+b_c), n = a_c+b_c- 1)
    rmap_pred <- preddist(rmap, n=n_c)
    p_lower <- pmix(rmap_pred, y_c)
    ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
  }
  # View(tibble(w, ppp))
  pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(ppp==max(ppp) & pass)
  w_eb <- ifelse(nrow(pppdt) == 0, 1, max(pppdt$w))
  ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))
  tt <- tibble(y_c, w_eb, ppp_eb) %>% rename(y=y_c)
  if(y==0){wdt <- tt}
  else {wdt <- rbind(wdt, tt)}
}  
# ggplot(wdt, aes(x=y, y = w_eb)) + geom_line()

#################################
#################################
ncores <- min(parallel::detectCores(), 40)
cl = makeCluster(ncores-1)
registerDoParallel(cl)

for(p in 1:length(pvec_c)){
  p_c <- pvec_c[p]
  
  #parallel computing at simulation level
  results <-
    foreach(s = 1:nsim, .combine = rbind, .packages = c("RBesT","tidyverse"), .errorhandling = "remove") %dopar% {
    set.seed(s+712)
    #current treatment arm
    y_t <- rbinom(1, n_t, p_c + effsize)
    postmix_t <- postmix(prior_t,r = y_t, n = n_t)
    post_t <- rmix(mix = postmix_t, n = npostdist)
    
    #current control arm
    y_c <- rbinom(1, n_c, p_c)
    #extract w_eb corresponding to observed data
    w_eb <- as.numeric(wdt %>% filter(y==y_c) %>% select(w_eb))
    
    w_rmap <- c(w_eb, w_vec) #first element is w_eb
    
    dvec <- estvec <- rep(NA, length(w_rmap))
    for(w in 1:length(w_rmap)){
      w_v <- w_rmap[w]
      #robustification, note that RBesT uses the mean/sample size-1 parametrization, see help of "robustify" for details
      rmap_c <- robustify(map_hat, weight = w_v, mean = a_c/(a_c+b_c), n = a_c+b_c- 1)
      postmix_rmap_c <- postmix(rmap_c,r = y_c, n = n_c)
      post_rmap_c <- rmix(mix = postmix_rmap_c, n = npostdist)
      estvec[w] <- median(post_rmap_c)
      
      post_diff <- post_t - post_rmap_c
      dvec[w] <- (mean(post_diff> Qcut) > Pcut)
    }
    
    c(dvec, estvec)
    
  }
  
  tt <- tibble(Rate=p_c, Method = c("EB", w_vec),
               PoS = colMeans(results[,1:(length(w_vec)+1)], na.rm = T),
               Est = colMeans(results[,-(1:(length(w_vec)+1))], na.rm = T), 
               Bias = Est-Rate, 
               Var = apply(results[,-(1:(length(w_vec)+1))], 2, var, na.rm=T),
               MSE = Bias^2 + Var)
  if(p==1) {outdt <- tt}
  else {outdt <- rbind(outdt, tt)}
  
}


ggplot(outdt, aes(x=Rate, y=PoS, color = Method)) + geom_line()
ggplot(outdt, aes(x=Rate, y=Bias, color = Method)) + geom_line()
ggplot(outdt, aes(x=Rate, y=MSE, color = Method)) + geom_line()

save.image("BinSim.RData")

