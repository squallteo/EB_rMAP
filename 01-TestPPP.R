library(RBesT)
set.seed(712)
map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | study,
                 data=AS,
                 tau.dist="HalfNormal",
                 tau.prior=1,
                 beta.prior=2,
                 family=binomial)
map_hat <- automixfit(map_mcmc)
rmap <- robustify(map_hat, weight=0.2, mean=1/2)
rmap_pred <- preddist(rmap, n=10)
summary(rmap_pred)
plot(rmap_pred)

pmix(rmap_pred, 4)

