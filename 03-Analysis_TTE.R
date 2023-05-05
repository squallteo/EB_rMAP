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

#one interval analysis from year 0 to 1.5
rvec <- colSums(rmat[1:6,])
Evec <- colSums(Emat[1:6,])

histdt <- tibble(study = 1:9, r = rvec[1:9], exp = Evec[1:9])
meta::metarate(r, exp, study, data=histdt)

map_mcmc <- gMAP(r ~ 1 + offset(log(exp)) | study, data = histdt, family = poisson, 
                 tau.dist = "HalfNormal", tau.prior = cbind(0, 0.5),
                 beta.prior=cbind(0, 10))
map_hat <- automixfit(map_mcmc)
vague <- mixgamma(c(1, summary(map_hat)[1], 1), param = "mn")

sample_MAP <- tibble(x=rmix(map_hat,100000), Prior = "MAP")
sample_vague <- tibble(x=rmix(vague,100000), Prior = "Vague")

(p1 <-
rbind(sample_MAP, sample_vague) %>% 
  ggplot(aes(x=x, fill=Prior)) + geom_density(alpha=0.5) +
  scale_x_continuous(limits = c(0,5)) + xlab("Hazard Rate") +
  theme_bw() +
  theme(axis.title = element_text(face="bold",size=20),
        axis.text = element_text(size=20),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.position = c(0.5, 0.5)
  )
)


ppp_cutvec <- c(0.85, 0.9, 0.95)
y_cvec <- seq(0, 51, 1)
n_cvec <- c(Evec[10])

for(j in 1:length(n_cvec)){
  n_c <- n_cvec[j]
  for(k in 1:length(ppp_cutvec)){
    ppp_cut <- ppp_cutvec[k]
    for(y in 1:length(y_cvec)){
      y_c <- y_cvec[y]
      #calculate ppp and pick the optimal w_V
      w <- seq(0, 1, 0.01)
      ppp <- rep(NA, length(w))
      for(i in 1:length(w)){
        rmap <- robustify(map_hat, weight=w[i], mean=summary(map_hat)[1], n=1)
        rmap_pred <- preddist(rmap, n=n_c)
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


(p2 <-
outdt %>% mutate(hrate = y_c/SS) %>% filter(SS==117.6) %>%
  ggplot(aes(x=y_c, y = w_eb, color=factor(Gamma))) + geom_line(size=1) + geom_vline(xintercept=0.37*n_c, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 50, 5)) + 
  xlab("Observed Number of Events") + ylab("EB-rMAP Weight") + theme_bw() +  
  scale_color_discrete(name="Gamma") + 
  theme(axis.title = element_text(face="bold",size=20),
        axis.text = element_text(size=20),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.position = c(0.25, 0.25)
  )
    
)


# png(filename = "Analysis1.png", width = 1200, height = 600)
# gridExtra::grid.arrange(p1, p2, nrow = 1)
# dev.off()

outdt %>% filter(y_c==rvec[10])



#use gamma=0.9, w_eb = 0.54 when r = 32
w_eb <- 0.54
#EB
rmap <- robustify(map_hat, weight=w_eb, mean=summary(map_hat)[4], n=1)
postmix_rmap <- postmix(rmap, n = Evec[10], m = rvec[10]/Evec[10])
post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
tt <- round(quantile(post_rmap_c, c(0.5, 0.025, 0.975)), 3)
res1 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
dt1 <- tibble(Method="EB-rMAP", Sample=post_rmap_c)
#MAP
rmap <- robustify(map_hat, weight=0, mean=summary(map_hat)[4], n=1)
postmix_rmap <- postmix(rmap, n = Evec[10], m = rvec[10]/Evec[10])
post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
tt <- round(quantile(post_rmap_c, c(0.5, 0.025, 0.975)), 3)
res2 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
dt2 <- tibble(Method="MAP", Sample=post_rmap_c)
#Vague
rmap <- robustify(map_hat, weight=1, mean=summary(map_hat)[4], n=1)
postmix_rmap <- postmix(rmap, n = Evec[10], m = rvec[10]/Evec[10])
post_rmap_c <- rmix(mix = postmix_rmap, n = 20000)
tt <- round(quantile(post_rmap_c, c(0.5, 0.025, 0.975)), 3)
res3 <- paste(tt[1], " (", tt[2], ", ", tt[3], ")", sep="")
dt3 <- tibble(Method="Vague", Sample=post_rmap_c)

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

##############################################
##############################################
##############################################





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
