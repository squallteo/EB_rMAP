rm(list=ls())

package2load <- c("RBesT", "tidyverse", "ggplot2")
lapply(package2load, require, character.only = TRUE)
set.seed(712)

ppp_cut <- 0.75

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

rmat <- matrix(FIOCCO.n.events, nrow = 12, ncol = 10, byrow = F)
Emat <- matrix(FIOCCO.exp.time, nrow = 12, ncol = 10, byrow = F)

map_lst <- NULL
for(l in 1:12){
  #meta analysis of historical data
  histdt <- tibble(study = 1:9, r = rmat[l, 1:9], exp = Emat[l, 1:9])
  #MAP prior
  map_mcmc <- gMAP(r ~ 1 + offset(log(exp)) | study, data = histdt, family = poisson, 
                   tau.dist = "HalfNormal", tau.prior = cbind(0, 0.5),
                   beta.prior=cbind(0, 10))
  map_hat <- automixfit(map_mcmc)
  map_lst[[l]] <- map_hat
  #vague prior
  # f_v <- mixgamma(c(1,hrate,1), param = "mn", likelihood = "poisson")
}

# save.image("TTEAnalysis_GP.RData")

w_eb <- rep(-1, 12)
for(l in 1:12){
  histdt <- tibble(study = 1:9, r = rmat[l, 1:9], exp = Emat[l, 1:9])
  metafit <- meta::metarate(r, exp, data=histdt)
  hrate <- round(exp(metafit$TE.fixed),3)
  crate <- round(rmat[l, 10]/Emat[l, 10],3)
  map_hat <- map_lst[[l]]
  #calculate ppp and pick the optimal w_V
  w <- seq(0, 1, 0.01)
  ppp <- rep(NA, length(w))
  for(i in 1:length(w)){
    rmap <- robustify(map_hat, weight=w[i], mean=summary(map_hat)[4], n=1)
    rmap_pred <- preddist(rmap, n=rmat[l, 10])
    p_lower <- pmix(rmap_pred, crate)
    ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
  }
  pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
  w_eb[l] <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
  
  #EB
  rmap <- robustify(map_hat, weight=w_eb[l], mean=summary(map_hat)[4], n=1)
  postmix_rmap <- postmix(rmap, n = Emat[l, 10] , m = rmat[l, 10]/Emat[l, 10])
  post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
  tt <- round(quantile(post_rmap_c, c(0.5, 0.025, 0.975)), 3)
  res1 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
  #MAP
  rmap <- robustify(map_hat, weight=0, mean=summary(map_hat)[4], n=1)
  postmix_rmap <- postmix(rmap, n = Emat[l, 10] , m = rmat[l, 10]/Emat[l, 10])
  post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
  tt <- round(quantile(post_rmap_c, c(0.5, 0.025, 0.975)), 3)
  res2 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
  #Vague
  rmap <- robustify(map_hat, weight=1, mean=summary(map_hat)[4], n=1)
  postmix_rmap <- postmix(rmap, n = Emat[l, 10] , m = rmat[l, 10]/Emat[l, 10])
  post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
  tt <- round(quantile(post_rmap_c, c(0.5, 0.025, 0.975)), 3)
  res3 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
  
  # out_int <- tibble(Interval = l, MetaMean = hrate, HistMean = summary(map_hat)[1], HistMedian = summary(map_hat)[4], HistSD = summary(map_hat)[2],
  #                   Curr = crate, w_eb = w_eb[l], EB = res1, MAP = res2, Vague = res3)
  out_int <- tibble(Interval = l, HistMedian = round(summary(map_hat)[4], 3),
                    Curr = crate, w_eb = w_eb[l], EB = res1, MAP = res2, Vague = res3)
  if(l==1){resdt <- out_int}
  else(resdt <- rbind(resdt, out_int))
}
resdt

kableExtra::kbl(resdt, format="latex")

#########################################
#########################################
#########################################
#plot
l <- 6
map_hat <- map_lst[[l]]

#EB
rmap <- robustify(map_hat, weight=w_eb[l], mean=summary(map_hat)[4], n=1)
postmix_rmap <- postmix(rmap, n = Emat[l, 10] , m = rmat[l, 10]/Emat[l, 10])
post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
dt1 <- tibble(Method="EB-rMAP", Sample=post_rmap_c)
#MAP
rmap <- robustify(map_hat, weight=0, mean=summary(map_hat)[4], n=1)
postmix_rmap <- postmix(rmap, n = Emat[l, 10] , m = rmat[l, 10]/Emat[l, 10])
post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
dt2 <- tibble(Method="MAP (EX)", Sample=post_rmap_c)
#Vague
rmap <- robustify(map_hat, weight=1, mean=summary(map_hat)[4], n=1)
postmix_rmap <- postmix(rmap, n = Emat[l, 10] , m = rmat[l, 10]/Emat[l, 10])
post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
dt3 <- tibble(Method="Vague (NEX)", Sample=post_rmap_c)

plotdt <- rbind(dt1, dt2, dt3)


png("HazardDensity.png", width = 1200, height = 1000, res = 300)
plotdt %>% 
  ggplot(aes(x=Sample, color=Method)) + geom_density(size=1) +
  xlab("Hazard Rate") + ylab("Density") + theme_bw() +
  scale_x_continuous(labels = seq(0,1.5, 0.25), breaks = seq(0,1.5, 0.25)) +
  theme(legend.position = c(0.7, 0.7),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10))
dev.off()
