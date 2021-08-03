library(MASS)
library(coda)
library(R2jags)
library(survival)
library(RBesT)
library(ggplot2)
library(tidyverse)

source("03-BUGSModel.R")
source("03-Func.R")
##FIOCCO data set 
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
                                       prior.mu.sd.nex   = matrix(rep(1,120), nrow=10, ncol=12),
                                       p.exch            = matrix(rep(1,120), nrow=10),
                                       Prior.beta        = matrix(c(0,10), nrow=2) , 
                                       beta.cutoffs      = 0,
                                       Prior.tau.study   = c(0, 0.5), 
                                       Prior.tau.time    = c(-1.386294, 0.707293),
                                       # MAP.prior         = TRUE,
                                       pars              = c("tau.study","log.hazard.pred"),
                                       R.seed            = 10,
                                       bugs.seed         = 12
)

l <- 2
loghaz <- FIOCCO.MAP.Prior$BUGSoutput$sims.list$log.hazard.pred[,,l]
normal_mix <- RBesT::automixfit(loghaz, type="norm")
sigma(normal_mix) <- 1
ess(normal_mix)
sample1 <- exp(rmix(normal_mix,10000))


#Poisson-Gamma model?
haz <- exp(loghaz)
poisson_mix <- RBesT::automixfit(haz, type="gamma")
likelihood(poisson_mix) <- "poisson"
preddist(poisson_mix, n=1)

mixgamma(c(1,1,1), param =  ,likelihood = "poisson")

sample2 <- rmix(poisson_mix,10000)
tb2 <- tibble(Mix="Gamma", x=sample2)
tb1 <- tibble(Mix="Normal", x=sample1)

plotdt <- rbind(tb1, tb2)
plotdt %>%
  ggplot(aes(x=x, color=Mix)) + geom_density()
