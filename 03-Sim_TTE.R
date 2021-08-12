rm(list=ls())

package2load <- c("RBesT", "tidyverse", "doParallel", "ggplot2")
lapply(package2load, require, character.only = TRUE)

hist_exp <- c(5, 10, 15, 20)
hist_haz <- 0.4
curr_exp <- 30
curr_haz <- c(0.3, 0.4, 0.5, 0.6, 0.7)

Qcut <- 0.5
Pcut <- 0.9

ppp_cut <- 0.75
w_vec <- c(0, 0.5, 1)

npostdist <- 20000
nsim <- 1000

#################################
#################################
ncores <- min(parallel::detectCores(), 40)
cl <- makeCluster(ncores-1)
registerDoParallel(cl)

for(h in 1:length(curr_haz)){
  print(h)
  results <- 
    foreach(s = 1:nsim, .combine = rbind, .packages = c("RBesT","tidyverse"), .errorhandling = "remove") %dopar% {
      set.seed(s+712)
      hist_r <- rpois(rep(1,4), hist_haz*hist_exp)
      histdt <- tibble(study = 1:4, r = hist_r, exp = hist_exp)
      #MAP prior
      map_mcmc <- gMAP(r ~ 1 + offset(log(exp)) | study, data = histdt, family = poisson, 
                       tau.dist = "HalfNormal", tau.prior = cbind(0, 0.5),
                       beta.prior=cbind(0, 10))
      map_hat <- automixfit(map_mcmc)
      
      curr_r <- rpois(1, curr_haz[h]*curr_exp)
      
      w <- seq(0, 1, 0.01)
      ppp <- rep(NA, length(w))
      for(i in 1:length(w)){
        rmap <- robustify(map_hat, weight=w[i], mean=summary(map_hat)[4], n=1)
        rmap_pred <- preddist(rmap, n=curr_r)
        p_lower <- pmix(rmap_pred, curr_r/curr_exp)
        ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
      }
      pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
      w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
      
      w_rmap <- c(w_eb, w_vec) #first element is w_eb
      
      dvec <- estvec <- rep(NA, length(w_rmap))
      for(w in 1:length(w_rmap)){
        w_v <- w_rmap[w]
        rmap <- robustify(map_hat, weight=w_v, mean=summary(map_hat)[4], n=1)
        postmix_rmap <- postmix(rmap, n = curr_exp , m = curr_r/curr_exp)
        post_rmap_c <- rmix(mix = postmix_rmap, n = npostdist)
        estvec[w] <- median(post_rmap_c)
        dvec[w] <- (mean(post_rmap_c < Qcut) > Pcut)
      }
      
      c(dvec, estvec)
    }
  
  tt <- tibble(Hazard=curr_haz[h], Method = c("EB", w_vec),
               PoS = colMeans(results[,1:(length(w_vec)+1)], na.rm = T),
               Est = colMeans(results[,-(1:(length(w_vec)+1))], na.rm = T), 
               Bias = Est-Hazard, 
               Var = apply(results[,-(1:(length(w_vec)+1))], 2, var, na.rm=T),
               MSE = Bias^2 + Var)
  if(h==1) {outdt <- tt}
  else {outdt <- rbind(outdt, tt)}
}

# save.image("TTESim.RData")


# for(i in 1:3){
#   tt <- outdt %>% filter(Hazard == curr_haz[i]) %>% mutate(Bias=round(abs(Bias*1000),1), MSE=round(MSE*1000,1)) %>% select(PoS, Bias, MSE)
#   if(i==1){p1 <- tt}
#   else{p1 <- cbind(p1, tt)}
# }
# write.csv(p1, "p1.csv")
# 
# for(i in 4:5){
#   tt <- outdt %>% filter(Hazard == curr_haz[i]) %>% mutate(Bias=round(abs(Bias*1000),1), MSE=round(MSE*1000,1)) %>% select(PoS, Bias, MSE)
#   if(i==4){p2 <- tt}
#   else{p2 <- cbind(p2, tt)}
# }
# write.csv(p2, "p2.csv")
