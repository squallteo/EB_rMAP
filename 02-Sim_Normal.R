rm(list=ls())

package2load <- c("RBesT", "tidyverse", "doParallel", "ggplot2")
lapply(package2load, require, character.only = TRUE)

set.seed(712)

#current control
n_c <- 50
yvec_c <- seq(-60, -40, 0.1) #a fine grid of y_c for w_eb
muvec_c <- seq(-60, -40, 1) #a grid of mu_c for simulations
# w_vec <- c(0, 0.25, 0.5, 0.75, 1) #w_v in rMAP
w_vec <- c(0, 0.25, 0.5, 0.75, 1) #w_v in rMAP
m_v <- -50; sd_v <- 50 #mean and sd of vague beta prior for mu_c
#EB rMAP parameters
ppp_cut <- 0.9

#decision rule to claim trial success: pr(p_c < Qcut) > Pcut
Qcut <- -40
Pcut <- 0.975

npostdist <- 20000
nsim <- 5000
#################################
#################################
dt <- crohn[-5,]
sigma <- 40
dt$se_yh <- sigma/sqrt(dt$n)
#meta analysis of historical data
# with(dt, meta::metamean(n = n, mean = y, sd = rep(sigma,length(study))))
#point est is -50. Decide to consider true current mean from -60 to -40
#derive and approximate MAP prior
map_mcmc <- gMAP(cbind(y, se_yh) ~ 1 | study, 
                 #weight = n,
                 data=dt,
                 family=gaussian,
                 beta.prior=cbind(-50, sigma),
                 tau.dist="HalfNormal",tau.prior=cbind(0,5))
#approximate the MAP
map_hat <- automixfit(map_mcmc)
sigma(map_hat) <- sigma
ess(map_hat)

#################################
#################################
#for normal outcome, create a find grid for y_c and map the realized random variable to the nearest grid value to determine w_eb
for(y in 1:length(yvec_c)){
  y_c <- yvec_c[y]
  #calculate ppp and pick the optimal w_V
  w <- seq(0, 1, 0.01)
  ppp <- rep(NA, length(w))
  for(i in 1:length(w)){
    rmap <- robustify(map_hat, weight=w[i], mean=m_v, n=1, sigma=sigma)
    rmap_pred <- preddist(rmap, n=n_c, sigma=sigma)
    p_lower <- pmix(rmap_pred, y_c)
    ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
  }
  # View(tibble(w, ppp))
  # pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(ppp==max(ppp) & pass)
  # w_eb <- ifelse(nrow(pppdt) == 0, 1, max(pppdt$w))
  pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
  w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
  ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))
  tt <- tibble(y_c, w_eb, ppp_eb)
  if(y==1){wdt <- tt}
  else {wdt <- rbind(wdt, tt)}
} 
ggplot(wdt, aes(x=y_c, y = w_eb)) + geom_line()

#################################
#################################
ncores <- min(parallel::detectCores(), 40)
cl = makeCluster(ncores-1)
registerDoParallel(cl)

for(p in 1:length(muvec_c)){
  mu_c <- muvec_c[p]
  
  #parallel computing at simulation level
  results <-
    foreach(s = 1:nsim, .combine = rbind, .packages = c("RBesT","tidyverse"), .errorhandling = "remove") %dopar% {
      set.seed(s+712)
      #current control arm
      yall_c <- rnorm(n_c, mu_c, sigma)
      y_c <- mean(yall_c)
      #extract w_eb corresponding to observed data
      #if y_c is out of the grid, set w_eb to 1 (way off)
      if(y_c < min(yvec_c) | y_c > max(yvec_c)){
        w_eb <- 1 
      }
      else{
        tt <- abs(yvec_c - y_c)
        tt1 <- which(tt==min(tt))
        w_eb <- as.numeric(wdt %>% filter(y_c==yvec_c[tt1]) %>% select(w_eb))
      }
      w_rmap <- c(w_eb, w_vec) #first element is w_eb

      dvec <- estvec <- rep(NA, length(w_rmap))
      for(w in 1:length(w_rmap)){
        w_v <- w_rmap[w]
        rmap_c <- robustify(map_hat, weight = w_v, mean = m_v, n=1, sigma = sigma)
        postmix_rmap_c <- postmix(rmap_c, n = n_c, m = y_c)
        post_rmap_c <- rmix(mix = postmix_rmap_c, n = npostdist)
        estvec[w] <- median(post_rmap_c)
        dvec[w] <- (mean(post_rmap_c < Qcut) > Pcut)
      }
      
      c(dvec, estvec)
      
    }
  
  tt <- tibble(Mean=mu_c, Method = c("EB", w_vec),
               PoS = colMeans(results[,1:(length(w_vec)+1)], na.rm = T),
               Est = colMeans(results[,-(1:(length(w_vec)+1))], na.rm = T), 
               Bias = Est-Mean, 
               Var = apply(results[,-(1:(length(w_vec)+1))], 2, var, na.rm=T),
               MSE = Bias^2 + Var)
  if(p==1) {outdt <- tt}
  else {outdt <- rbind(outdt, tt)}
  
}

plotdt <- outdt %>% filter(!(Method %in% c("0.25", "0.75")))
ggplot(plotdt, aes(x=Mean, y=Bias, color = Method)) + geom_line()
ggplot(plotdt, aes(x=Mean, y=MSE, color = Method)) + geom_line()
ggplot(plotdt, aes(x=Mean, y=PoS, color = Method)) + geom_line()

# save.image("NormalSim.RData")

