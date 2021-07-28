rm(list=ls())

package2load <- c("RBesT", "tidyverse", "doParallel", "ggplot2")
lapply(package2load, require, character.only = TRUE)

set.seed(712)

n_c <- 20
p_cvec <- seq(0.15, 0.35, 0.01)
a_c <- b_c <- 1
es <- 0.2
ppp_cut <- 0.1 #c for ppp
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
  tt <- tibble(y_c, w_eb, ppp_eb)
  if(y==0){wdt <- tt}
  else {wdt <- rbind(wdt, tt)}
}  
  
ggplot(wdt, aes(x=y_c, y = w_eb)) + geom_line()