library(RBesT)
library(ggplot2)
library(tidyverse)
set.seed(712)

# with(AS, meta::metaprop(event = r, n = n, method = "Inverse"))
#historical rate around 0.25 by meta analysis
map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | study,
                 data=AS,
                 tau.dist="HalfNormal",
                 tau.prior=1,
                 beta.prior=2,
                 family=binomial)
map_hat <- automixfit(map_mcmc)

w <- seq(0, 1, 0.01)
ppp <- rep(NA, length(w))
for(i in 1:length(w)){
  rmap <- robustify(map_hat, weight=w[i], mean=1/2)
  rmap_pred <- preddist(rmap, n=80)
  p_lower <- pmix(rmap_pred, 36)
  ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
}

plotdt <- tibble(w, ppp)
plotdt %>% filter(ppp==max(ppp))
ggplot(plotdt, aes(x=w, y=ppp)) + geom_line()

