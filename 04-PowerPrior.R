rm(list=ls())

library(BayesPPD)
library(RBesT)
library(tidyverse)
library(doParallel)

set.seed(712)

#current control
n_c <- 50
pvec_c <- seq(0.10, 0.50, 0.02)
a0 <- 0.5

#current treatment
n_t <- 100
es <- 0
# es <- 0.2
#prior of treatment
prior_t <- c(1,1)

#decision rule to claim trial success: pr(p_t - p_c > Qcut) > Pcut
Qcut <- 0
Pcut <- 0.95
# success_rule = decision2S(pc = 0.95, qc = 0, lower.tail = F, link = "identity")
npostdist <- 20000
nsim <- 5000

data(AS)
ppdt <- AS %>% select(r, n) %>% mutate(a0 = a0)

ncores <- min(parallel::detectCores(), 40)
cl <- makeCluster(ncores-1)
registerDoParallel(cl)

for(p in 1:length(pvec_c)){
  p_c <- pvec_c[p]
  
  results <-
    foreach(s = 1:nsim, .combine = rbind, .packages = c("BayesPPD"), 
            .errorhandling = "remove") %dopar% {
              #current control arm
              y_c <- rbinom(1, n_c, p_c)
              ppfit <- 
                two.grp.fixed.a0(data.type = "Bernoulli", y.c = y_c, n.c = n_c,
                                 historical = as.matrix(ppdt))
              post_pp_c <- rbeta(npostdist, ppfit[1], ppfit[2])
              #current treatment arm
              y_t <- rbinom(1, n_t, p_c+es)
              post_t <- rbeta(npostdist, prior_t[1] + y_t, prior_t[2] + n_t - y_t)
              
              tt <- post_t - post_pp_c
              
              c(median(tt), (mean(tt > Qcut) > Pcut))
            } 
  tt <- tibble(Rate=p_c, Method = "PP",
               PoS = mean(results[,2], na.rm = T),
               Est = mean(results[,1], na.rm = T), 
               Bias = Est-es, 
               Var = var(results[,1], na.rm=T),
               MSE = Bias^2 + Var,
               Nsim = nrow(results))
  if(p==1) {outdt <- tt}
  else {outdt <- rbind(outdt, tt)}
}


plotdt <- outdt %>% filter(Rate %in% c(0.1, 0.2, 0.26, 0.36, 0.5))
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
  
  png("Sim_Bin_PP_Null.png", width = 2700, height = 1000, res = 300)
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





