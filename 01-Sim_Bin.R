rm(list=ls())

package2load <- c("RBesT", "tidyverse", "doParallel", "ggplot2")
lapply(package2load, require, character.only = TRUE)

set.seed(712)

#current control
n_c <- 50
pvec_c <- seq(0.10, 0.50, 0.02)
w_vec <- c(0, 0.5, 1) #w_v in rMAP
a_c <- 1; b_c <- 1 #vague beta prior for p(c)
#EB rMAP parameters
ppp_cut <- 0.8

#current treatment
n_t <- 100
es <- 0
# es <- 0.2
#prior of treatment
prior_t <- mixbeta(c(1, 1, 1))

#decision rule to claim trial success: pr(p_t - p_c > Qcut) > Pcut
Qcut <- 0
Pcut <- 0.95
# success_rule = decision2S(pc = 0.95, qc = 0, lower.tail = F, link = "identity")
npostdist <- 20000
nsim <- 5000
#################################
#################################
#meta analysis of historical data
# with(dt, meta::metaprop(event = r, n = n, method = "Inverse"))
#point est is 0.25 (0.202, 0.312). Decide to consider true current rate from 0.2 to 0.32
#derive and approximate MAP prior
map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | study,
                 data=AS,
                 tau.dist="HalfNormal",
                 tau.prior=1,
                 beta.prior=2,
                 family=binomial)
map_hat <- automixfit(map_mcmc)
ess(map_hat)

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
  # pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(ppp==max(ppp) & pass)
  # w_eb <- ifelse(nrow(pppdt) == 0, 1, max(pppdt$w))
  pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
  w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
  ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))
  tt <- tibble(y_c, w_eb, ppp_eb) %>% rename(y=y_c)
  if(y==0){wdt <- tt}
  else {wdt <- rbind(wdt, tt)}
}  
ggplot(wdt, aes(x=y, y = w_eb)) + geom_line()

#################################
#################################
ncores <- min(parallel::detectCores(), 40)
cl <- makeCluster(ncores-1)
registerDoParallel(cl)

for(p in 1:length(pvec_c)){
  p_c <- pvec_c[p]
  
  #parallel computing at simulation level
  results <-
    foreach(s = 1:nsim, .combine = rbind, .packages = c("RBesT","tidyverse"), .errorhandling = "remove") %dopar% {
      set.seed(s+712)
      #current control arm
      y_c <- rbinom(1, n_c, p_c)
      #current treatment arm
      y_t <- rbinom(1, n_t, p_c+es)
      postmix_t <- postmix(prior_t, r = y_t, n = n_t)
      post_t <- rmix(mix = postmix_t, n = npostdist)
      
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
        
        tt <- post_t - post_rmap_c
        estvec[w] <- median(tt)
        dvec[w] <- (mean(tt > Qcut) > Pcut)
      }
      
      c(dvec, estvec)
      # c(y_c, w_rmap[1], dvec, estvec)
      
    }
  
  tt <- tibble(Rate=p_c, Method = c("EB", w_vec),
               PoS = colMeans(results[,1:(length(w_vec)+1)], na.rm = T),
               Est = colMeans(results[,-(1:(length(w_vec)+1))], na.rm = T), 
               Bias = Est-es, 
               Var = apply(results[,-(1:(length(w_vec)+1))], 2, var, na.rm=T),
               MSE = Bias^2 + Var,
               Nsim = nrow(results))
  if(p==1) {outdt <- tt}
  else {outdt <- rbind(outdt, tt)}
  
}


plotdt <- outdt %>% filter(!(Method %in% c("0.25", "0.75")))

plotlst <- list()
library(ggpubr)

if(es==0){
  plotlst[[1]] <- 
    ggplot(plotdt, aes(x=Rate, y=PoS, color = Method)) + geom_line(size=1) + geom_vline(xintercept=0.25, linetype="dashed") +
    xlab("True Current Rate") + ylab("Error Rate") + theme_bw() + 
    scale_color_discrete(name="Mixture\nWeight") +
    theme(legend.position = "none")
  
  plotlst[[2]] <- 
    ggplot(plotdt, aes(x=Rate, y=abs(Bias), color = Method)) + geom_line(size=1) + geom_vline(xintercept=0.25, linetype="dashed") +
    xlab("True Current Rate") + ylab("Absolute Bias") + theme_bw() + 
    scale_color_discrete(name="Mixture\nWeight") +
    theme(legend.position = "none")
  
  plotlst[[3]] <- 
    ggplot(plotdt, aes(x=Rate, y=MSE, color = Method)) + geom_line(size=1) + geom_vline(xintercept=0.25, linetype="dashed") +
    xlab("True Current Rate") + ylab("Mean Square Error") + theme_bw() + 
    scale_color_discrete(name="Mixture\nWeight")
  
  png("Sim_Bin_Null.png", width = 2700, height = 1000, res = 300)
  ggarrange(plotlst[[1]], plotlst[[2]], plotlst[[3]],
            nrow = 1, ncol = 3)
  dev.off()
  
  # save.image("BinSim_Null.RData")
} else{
  plotlst[[1]] <- 
    ggplot(plotdt, aes(x=Rate, y=PoS, color = Method)) + geom_line(size=1) + geom_vline(xintercept=0.25, linetype="dashed") +
    xlab("True Current Rate") + ylab("Probability of Success") + theme_bw() + 
    scale_color_discrete(name="Mixture\nWeight") +
    theme(legend.position = "none")
  
  plotlst[[2]] <- 
    ggplot(plotdt, aes(x=Rate, y=abs(Bias), color = Method)) + geom_line(size=1) + geom_vline(xintercept=0.25, linetype="dashed") +
    xlab("True Current Rate") + ylab("Absolute Bias") + theme_bw() + 
    scale_color_discrete(name="Mixture\nWeight") +
    theme(legend.position = "none")
  
  plotlst[[3]] <- 
    ggplot(plotdt, aes(x=Rate, y=MSE, color = Method)) + geom_line(size=1) + geom_vline(xintercept=0.25, linetype="dashed") +
    xlab("True Current Rate") + ylab("Mean Square Error") + theme_bw() + 
    scale_color_discrete(name="Mixture\nWeight")
  
  png("Sim_Bin_Alt.png", width = 2700, height = 1000, res = 300)
  ggarrange(plotlst[[1]], plotlst[[2]], plotlst[[3]],
            nrow = 1, ncol = 3)
  dev.off()
  
  # save.image("BinSim_Alt.RData")
}



