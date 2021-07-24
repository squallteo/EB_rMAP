rm(list=ls())

package2load <- c("RBesT", "tidyverse", "doParallel")
lapply(package2load, require, character.only = TRUE)

set.seed(712)

n_c <- 60
p_cvec <- seq(0.15, 0.35, 0.01)
es <- 0.2
ppp_cut <- 0 #c for ppp
nsim <- 500
# with(AS, meta::metaprop(event = r, n = n, method = "Inverse"))
#historical rate around 0.25 by meta analysis
#derive and approximate MAP prior
map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | study,
                 data=AS,
                 tau.dist="HalfNormal",
                 tau.prior=1,
                 beta.prior=2,
                 family=binomial)
map_hat <- automixfit(map_mcmc)

for(p in 1:length(p_cvec)){
  p_c <- p_cvec[p]
  
  
  #parallel computing at simulation level
  for(s in 1:nsim){
    y_c <- rbinom(1, n_c, p_c)
    #calculate ppp and pick the optimal w_V
    w <- seq(0, 1, 0.01)
    ppp <- rep(NA, length(w))
    for(i in 1:length(w)){
      rmap <- robustify(map_hat, weight=w[i], mean=1/2)
      rmap_pred <- preddist(rmap, n=n_c)
      p_lower <- pmix(rmap_pred, y_c)
      ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
    }
    # View(tibble(w, ppp))
    pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp > ppp_cut)) %>% filter(ppp==max(ppp) & pass)
    w_eb <- ifelse(nrow(pppdt) == 0, 1, max(pppdt$w))
    ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))
    
    if(s==1){w_ebvec <- w_eb; ppp_ebvec <- ppp_eb; y_cvec <- y_c}
    else {w_ebvec <- c(w_ebvec, w_eb); ppp_ebvec <- c(ppp_ebvec, ppp_eb); y_cvec <- c(y_cvec, y_c)}
    
  }
  
  tibble(rate = p_c, y_c = y_cvec, w_eb = w_ebvec, ppp_eb = ppp_ebvec)

  
  
}

