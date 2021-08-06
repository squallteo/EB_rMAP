library(MASS)
library(coda)
library(R2jags)
library(survival)
library(RBesT)
library(ggplot2)
library(tidyverse)

rm(list=ls())
set.seed(712)

source("03-BUGSModel.R")
source("03-Func.R")
##FIOCCO data set from Roychoudhury and Neuenschwander (2020) paper in Statistics in Medicine
FIOCCO.n.events <- c(1,  3,  3,  4,  3,  0,  0,  2,  0,  6,  0,  0,  9,  1,  0, 10,  6,  6,  5,  9,
                     9,  3,  0,  0,  1,  3,  5,  7,  9,  4,  5, 10,  0,  0,  3,  7,  1,  2,  2,  4, 
                     3,  1,  3,  0,  0,  0,  0,  0,  5,  3,  6,  2,  3,  3,  0,  2,  1,  1,  1,  0, 
                     0,  6,  3, 12,  8,  2, 3,  2, 11,  1,  0, 10,  2,  2,  5,  3,  3,  3,  2,  3, 
                     3,  0,  0,  0,  0,  1,  3,  4, 1,  1,  4,  1,  6,  0,  0,  0,  2,  0,  3,  1, 
                     4,  0,  1,  1,  0,  0,  0,  0,  1,  5, 17, 0,  2,  7,  8,  4,  0,  6,  2,  0)
FIOCCO.exp.time <- c(  9.4,  8.8,  7.9,  7.0,  6.1,  5.8,  5.8,  7.3,  8.8,  7.6,  6.2, 10.0, 21.1, 
                       19.9, 19.8, 18.5, 16.5, 15.0, 13.6, 15.7, 16.2, 13.6, 12.5, 20.1, 21.9, 21.4, 
                       20.4, 18.9, 16.9, 15.2, 14.1, 16.2, 18.5, 18.3, 17.0, 24.5,  5.6,  5.2,  4.8, 
                       4.0,  3.1,  2.6,  2.1,  2.3,  2.9,   2.9,  2.9,  4.7,  6.4,  5.4,  4.2,  3.2, 
                       2.6,  1.9,  1.5,  1.7, 1.5,  1.0,  0.6,  0.7, 17.8, 17.0, 15.9, 14.0, 11.5, 
                       10.2,  9.6, 11.9, 12.4,  9.9,  9.4, 12.1,  8.0,  7.5,  6.6,  5.6,  4.9,  4.1, 
                       3.5,  3.8,  3.6,  2.9,  2.9,  4.7,  9.2,  9.1,  8.6,  7.8,  7.1,  6.9,  6.2,  
                       7.4,  8.0,  6.7,  6.6, 10.7,  5.2,  5.0,  4.6,  4.1,  3.5,  3.0,  2.9,  3.5, 
                       4.2,  4.2,  4.1,  6.7, 23.4, 22.6, 19.9, 17.8, 17.5, 16.4, 14.5, 17.2, 21.0, 
                       19.7, 17.4, 27.5)




FIOCCO.MAP.Prior <- MAC.Surv.anal_jags(Nobs              = 108,
                                       study             = sort(rep(1:9,12)),
                                       Nstudies          = 9,
                                       Nint              = 12,
                                       Ncov              = 1,
                                       Nout              = 1,
                                       int.low           = rep(1:12,9),
                                       int.high          = rep(1:12,9),
                                       int.length        = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                             0.33, 0.42, 0.42, 0.41, 0.67),
                                       n.events          = FIOCCO.n.events[1:108],
                                       exp.time          = FIOCCO.exp.time[1:108],
                                       X                 = matrix(0,108,1),
                                       Prior.mu.mean.ex  = c(0, 10),
                                       Prior.rho.ex      = c(0,10),
                                       prior.mu.mean.nex = matrix(rep(0,120), nrow=10, ncol=12),
                                       prior.mu.sd.nex   = matrix(rep(2,120), nrow=10, ncol=12),
                                       p.exch            = matrix(rep(1,120), nrow=10),
                                       Prior.beta        = matrix(c(0,10), nrow=2) ,
                                       beta.cutoffs      = 0,
                                       Prior.tau.study   = c(0, 0.5),
                                       Prior.tau.time    = c(-1.386294, 0.707293),
                                       MAP.prior         = TRUE,
                                       pars              = c("tau.study","log.hazard.pred"),
                                       R.seed            = 10,
                                       bugs.seed         = 12
)


map_lst <- NULL
ppp_lst <- NULL
ppp_cut <- 0.8
w_eb <- rep(-1, 12)

for(t in 1:12){
  # curr_r <- rdata[18+t]
  # curr_E <- Edata[18+t]
  curr_r <- FIOCCO.n.events[108+t]
  curr_E <- FIOCCO.exp.time[108+t]
  loghaz <- FIOCCO.MAP.Prior$BUGSoutput$sims.list$log.hazard.pred[,,t]
  map_hat <- RBesT::automixfit(loghaz, type="norm")
  sigma(map_hat) <- sd(loghaz)
  map_lst[[t]] <- map_hat
  
  w <- seq(0, 1, 0.01)
  ppp <- rep(NA, length(w))
  for(i in 1:length(w)){
    rmap <- robustify(map_hat, weight=w[i], mean=0, sigma=2)
    rmap_pred <- preddist(rmap, n=curr_r)
    p_lower <- pmix(rmap_pred, log(curr_r/curr_E))
    ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
  }
  pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
  ppp_lst[[t]] <- tibble(w, ppp) %>% ggplot(aes(x=w,y=ppp)) + geom_line() + ggtitle(t)
  w_eb[t] <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
}

#analysis with EM-rMAP, MAP prior and vague prior
analysis_EB <- MAC.Surv.anal_jags(Nobs              = 120,
                                  study             = sort(rep(1:10,12)),
                                  Nstudies          = 10,
                                  Nint              = 12,
                                  Ncov              = 1,
                                  Nout              = 1,
                                  int.low           = rep(1:12,10),
                                  int.high          = rep(1:12,10),
                                  int.length        = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                       0.33, 0.42, 0.42, 0.41, 0.67),
                                  n.events          = FIOCCO.n.events[1:120],
                                  exp.time          = FIOCCO.exp.time[1:120],
                                  X                 = matrix(0,120,1),
                                  Prior.mu.mean.ex  = c(0, 10),
                                  Prior.rho.ex      = c(0,10),
                                  prior.mu.mean.nex = matrix(rep(0,120), nrow=10, ncol=12),
                                  prior.mu.sd.nex   = matrix(rep(2,120), nrow=10, ncol=12),
                                  p.exch            = 1 - rep(1,10) %o% w_eb, #EB weights
                                  Prior.beta        = matrix(c(0,10), nrow=2) ,
                                  beta.cutoffs      = 0,
                                  Prior.tau.study   = c(0, 0.5),
                                  Prior.tau.time    = c(-1.386294, 0.707293),
                                  MAP.prior         = FALSE,
                                  pars              = c("tau.study","log.hazard.pred"),
                                  R.seed            = 10,
                                  bugs.seed         = 12
)

analysis_MAP <- MAC.Surv.anal_jags(Nobs              = 120,
                                  study             = sort(rep(1:10,12)),
                                  Nstudies          = 10,
                                  Nint              = 12,
                                  Ncov              = 1,
                                  Nout              = 1,
                                  int.low           = rep(1:12,10),
                                  int.high          = rep(1:12,10),
                                  int.length        = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                        0.33, 0.42, 0.42, 0.41, 0.67),
                                  n.events          = FIOCCO.n.events[1:120],
                                  exp.time          = FIOCCO.exp.time[1:120],
                                  X                 = matrix(0,120,1),
                                  Prior.mu.mean.ex  = c(0, 10),
                                  Prior.rho.ex      = c(0,10),
                                  prior.mu.mean.nex = matrix(rep(0,120), nrow=10, ncol=12),
                                  prior.mu.sd.nex   = matrix(rep(2,120), nrow=10, ncol=12),
                                  p.exch            = matrix(rep(1,120), nrow=10),
                                  Prior.beta        = matrix(c(0,10), nrow=2) ,
                                  beta.cutoffs      = 0,
                                  Prior.tau.study   = c(0, 0.5),
                                  Prior.tau.time    = c(-1.386294, 0.707293),
                                  MAP.prior         = FALSE,
                                  pars              = c("tau.study","log.hazard.pred"),
                                  R.seed            = 10,
                                  bugs.seed         = 12
)

analysis_vague <- MAC.Surv.anal_jags(Nobs              = 120,
                                   study             = sort(rep(1:10,12)),
                                   Nstudies          = 10,
                                   Nint              = 12,
                                   Ncov              = 1,
                                   Nout              = 1,
                                   int.low           = rep(1:12,10),
                                   int.high          = rep(1:12,10),
                                   int.length        = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                         0.33, 0.42, 0.42, 0.41, 0.67),
                                   n.events          = FIOCCO.n.events[1:120],
                                   exp.time          = FIOCCO.exp.time[1:120],
                                   X                 = matrix(0,120,1),
                                   Prior.mu.mean.ex  = c(0, 10),
                                   Prior.rho.ex      = c(0,10),
                                   prior.mu.mean.nex = matrix(rep(0,120), nrow=10, ncol=12),
                                   prior.mu.sd.nex   = matrix(rep(2,120), nrow=10, ncol=12),
                                   p.exch            = matrix(rep(0,120), nrow=10),
                                   Prior.beta        = matrix(c(0,10), nrow=2) ,
                                   beta.cutoffs      = 0,
                                   Prior.tau.study   = c(0, 0.5),
                                   Prior.tau.time    = c(-1.386294, 0.707293),
                                   MAP.prior         = FALSE,
                                   pars              = c("tau.study","log.hazard.pred"),
                                   R.seed            = 10,
                                   bugs.seed         = 12
)

# save.image("TTEAnalysis.RData")

for(i in 1:12){
  #EB
  loghaz <- analysis_EB$BUGSoutput$sims.list$log.hazard.pred[,,i]
  tt <- round(exp(c(median(loghaz), quantile(loghaz, c(0.025, 0.975)))), 3)
  res1 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
  #MAP
  loghaz <- analysis_MAP$BUGSoutput$sims.list$log.hazard.pred[,,i]
  tt <- round(exp(c(median(loghaz), quantile(loghaz, c(0.025, 0.975)))), 3)
  res2 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
  #Vague
  loghaz <- analysis_vague$BUGSoutput$sims.list$log.hazard.pred[,,i]
  tt <- round(exp(c(median(loghaz), quantile(loghaz, c(0.025, 0.975)))), 3)
  res3 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
  
  out_int <- tibble(Interval = i, w_eb = w_eb[i], EB = res1, MAP = res2, Vague = res3)
  
  if(i==1){resdt <- out_int}
  else(resdt <- rbind(resdt, out_int))
}

resdt



########################################################
########################################################
########################################################
########################################################
#consolidate data
# int1 <- 1:4
# int2 <- 5:12
# 
# tt <- matrix(FIOCCO.n.events, nrow = 12, ncol = 10, byrow = F)
# rdata <- c(rbind(colSums(tt[int1,]), colSums(tt[int2,])))
# 
# tt <- matrix(FIOCCO.exp.time, nrow = 12, ncol = 10, byrow = F)
# Edata <- c(rbind(colSums(tt[int1,]), colSums(tt[int2,])))

# FIOCCO.MAP.Prior <- MAC.Surv.anal_jags(Nobs              = 18,
#                                        study             = sort(rep(1:9,2)),
#                                        Nstudies          = 9,
#                                        Nint              = 2,
#                                        Ncov              = 1,
#                                        Nout              = 1,
#                                        int.low           = rep(1:2,9),
#                                        int.high          = rep(1:2,9),
#                                        int.length        = c(1,3),
#                                        n.events          = rdata[1:18],
#                                        exp.time          = Edata[1:18],
#                                        X                 = matrix(0,18,1),
#                                        Prior.mu.mean.ex  = c(0,10), 
#                                        Prior.rho.ex      = c(0,10),
#                                        # prior.mu.mean.nex = matrix(rep(0,18), nrow=9, ncol=2),
#                                        # prior.mu.sd.nex   = matrix(rep(1,18), nrow=9, ncol=2),
#                                        # p.exch            = matrix(rep(1,18), nrow=9),
#                                        prior.mu.mean.nex = matrix(rep(0,20), nrow=10, ncol=2),
#                                        prior.mu.sd.nex   = matrix(rep(2,20), nrow=10, ncol=2),
#                                        p.exch            = matrix(rep(1,20), nrow=10),
#                                        Prior.beta        = matrix(c(0,10), nrow=2), 
#                                        beta.cutoffs      = 0,
#                                        Prior.tau.study   = c(0, 0.5), 
#                                        Prior.tau.time    = c(-1.386294, 0.707293),
#                                        MAP.prior         = TRUE,
#                                        pars              = c("tau.study","log.hazard.pred"),
#                                        R.seed            = 10,
#                                        bugs.seed         = 12
# )


# #Poisson-Gamma model?
# haz <- exp(loghaz)
# poisson_mix <- RBesT::automixfit(haz, type="gamma")
# likelihood(poisson_mix) <- "poisson"
# preddist(poisson_mix, n=1)
# 
# mixgamma(c(1,1,1), param =  ,likelihood = "poisson")
# 
# sample2 <- rmix(poisson_mix,10000)
# tb2 <- tibble(Mix="Gamma", x=sample2)
# tb1 <- tibble(Mix="Normal", x=sample1)
# 
# plotdt <- rbind(tb1, tb2)
# plotdt %>%
#   ggplot(aes(x=x, color=Mix)) + geom_density()
