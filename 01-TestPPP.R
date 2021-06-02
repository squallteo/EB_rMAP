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
plot(rmap_pred)

pmix(rmap_pred, 4)


############
# Example 2: 2-comp Beta mixture

bm <- mixbeta( inf=c(0.8,15,50),rob=c(0.2,1,1))
plot(bm)
bmPred <- preddist(bm,n=10)
summary(bmPred)
plot(bmPred)
mdn <- qmix(bmPred,0.5)
mdn
d <- dmix(bmPred,x=0:10)

n.sim <- 100000
r <-  rmix(bmPred,n.sim)
d
table(r)/n.sim


# Example 3: 3-comp Normal mixture

m3 <- mixnorm( c(0.50,-0.2,0.1),c(0.25,0,0.2), c(0.25,0,0.5), sigma=10)
print(m3)
summary(m3)
plot(m3)
predm3 <- preddist(m3,n=2)
plot(predm3)
print(predm3)
summary(predm3)
