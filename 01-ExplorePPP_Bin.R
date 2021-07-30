rm(list=ls())

package2load <- c("RBesT", "tidyverse", "ggplot2")
lapply(package2load, require, character.only = TRUE)

set.seed(712)

ppp_cutvec <- c(0.8, 0.85, 0.9) #gamma for ppp
n_cvec <- c(25, 50, 100)
a_c <- b_c <- 1

# with(AS, meta::metaprop(event = r, n = n, method = "Inverse"))
#historical rate around 0.25 by meta analysis
#derive and approximate MAP prior
map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | study,
                 data=AS,
                 tau.dist="HalfNormal",
                 tau.prior=1,
                 beta.prior=2,
                 family=binomial)
map_hat <- automixfit(map_mcmc)#approximate the MAP
ess(map_hat)


for(k in 1:length(ppp_cutvec)){
  ppp_cut <- ppp_cutvec[k]
  for(j in 1:length(n_cvec)){
    n_c <- n_cvec[j]
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
      pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
      w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
      ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))
      tt <- tibble(y_c, w_eb, ppp_eb)
      if(y==0){wdt <- tt}
      else {wdt <- rbind(wdt, tt)}
    }  
    
    tt <- wdt %>% mutate(Gamma = ppp_cut, SS = n_c)
    if(k==1 & j==1){outdt <- tt}
    else{outdt <- rbind(outdt, tt)}
  }
}


plotlst <- list()
plotlst[[1]] <-
  outdt %>% filter(Gamma==0.85) %>%
  ggplot(aes(x=y_c, y = w_eb, color=factor(SS))) + geom_line(size=1) + geom_vline(xintercept=0.25*n_c, linetype="dashed") +
  xlab("Observed Mean Response") + ylab("EB-rMAP Weight") + theme_bw() + ggtitle("Gamma: 0.85") + 
  scale_color_discrete(name="Sample\nSize")

plotlst[[2]] <-
  outdt %>% filter(SS==50) %>%
  ggplot(aes(x=y_c, y = w_eb, color=factor(Gamma))) + geom_line(size=1) + geom_vline(xintercept=0.25*n_c, linetype="dashed") +
  xlab("Observed Mean Response") + ylab("EB-rMAP Weight") + theme_bw() + ggtitle("Current Sample Size: 50") + 
  scale_color_discrete(name="Gamma")

# save.image("CompareWeights_Normal.RData")

# library(ggpubr)
# png("EBweights.png", width = 2700, height = 1000, res = 300)
# ggarrange(plotlst[[1]], plotlst[[2]],
#           nrow = 1, ncol = 2)
# dev.off()
