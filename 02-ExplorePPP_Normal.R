rm(list=ls())

package2load <- c("RBesT", "tidyverse", "ggplot2")
lapply(package2load, require, character.only = TRUE)

set.seed(712)

ppp_cutvec <- c(0.8, 0.85, 0.9) #gamma for ppp
n_cvec <- c(25, 50, 100)
y_cvec <- seq(-60, -40, 0.1)
m_v <- -50; sd_v <- 50 #mean and sd of vague beta prior for mu_c
#derive and approximate MAP prior
dt <- crohn[-5,]
sigma <- 40
dt$se_yh <- sigma/sqrt(dt$n)
#meta analysis of historical data
# with(dt, meta::metamean(n = n, mean = y, sd = rep(sigma,length(study))))
#point est is -46.8. Decide to consider true current mean from -60 to -40
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

f_v <- mixnorm(c(1, -50, 1), sigma=sigma, param="mn")
ess(f_v)

for(k in 1:length(ppp_cutvec)){
  ppp_cut <- ppp_cutvec[k]
  for(j in 1:length(n_cvec)){
    n_c <- n_cvec[j]
    for(y in 1:length(y_cvec)){
      y_c <- y_cvec[y]
      #calculate ppp and pick the optimal w_V
      w <- seq(0, 1, 0.01)
      ppp <- rep(NA, length(w))
      for(i in 1:length(w)){
        rmap <- robustify(map_hat, weight=w[i], mean=m_v, n=1, sigma=sigma)
        rmap_pred <- preddist(rmap, n=n_c, sigma=sigma)
        p_lower <- pmix(rmap_pred, y_c)
        ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
      }
      pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
      w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
      ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))
      tt <- tibble(y_c, w_eb, ppp_eb)
      if(y==1){wdt <- tt}
      else {wdt <- rbind(wdt, tt)}
    }
    tt <- wdt %>% mutate(Gamma = ppp_cut, SS = n_c)
    if(k==1 & j==1){outdt <- tt}
    else{outdt <- rbind(outdt, tt)}
  }
}

plotlst <- list()
plotlst[[1]] <-
outdt %>% filter(Gamma==0.9) %>%
  ggplot(aes(x=y_c, y = w_eb, color=factor(SS))) + geom_line(size=1) + geom_vline(xintercept=-46.8, linetype="dashed") +
  xlab("Observed Mean Response") + ylab("EB-rMAP Weight") + theme_bw() +
  scale_color_discrete(name="Sample\nSize")

plotlst[[2]] <-
outdt %>% filter(SS==50) %>%
  ggplot(aes(x=y_c, y = w_eb, color=factor(Gamma))) + geom_line(size=1) + geom_vline(xintercept=-46.8, linetype="dashed") +
  xlab("Observed Mean Response") + ylab("EB-rMAP Weight") + theme_bw() +
  scale_color_discrete(name="Gamma")

# save.image("CompareWeights_Normal.RData")

# library(ggpubr)
# png("EBweights.png", width = 2700, height = 1000, res = 300)
# ggarrange(plotlst[[1]], plotlst[[2]],
#           nrow = 1, ncol = 2)
# dev.off()
