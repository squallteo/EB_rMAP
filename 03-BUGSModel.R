MAC.Surv.WB <- function()
{
  # tau.study: Between study variation (Half normal)
  # tau_k of equation (5) in paper
  Prior.tau.study.prec <- pow(Prior.tau.study[2],-2)
  for (t in 1:Nint) {
    tau.study[t] ~ dnorm(Prior.tau.study[1], Prior.tau.study.prec) %_% I(0, )
    tau.study.prec[t] <- pow(tau.study[t],-2)
  }
  #tau.time: correlation for piecewise exponential pieces (log-normal)
  Prior.tau.time.prec <- pow(Prior.tau.time[2],-2)
  tau.time ~ dlnorm(Prior.tau.time[1], Prior.tau.time.prec)
  tau.time.prec <- pow(tau.time,-2)
  # priors for regression parameter: h covariates
  for (h in 1:Ncov) {
    prior.beta.prec[h] <- pow(Prior.beta[2,h],-2)
    beta[h,1] ~ dnorm(Prior.beta[1,h],prior.beta.prec[h])
  }
  
  ###########################
  #EXNEX structure of mu (equation (7) in paper)
  ###########################
  #EX structure for mu is 1st order NDLM 
  mu.prec.ex <- pow(Prior.mu.mean.ex[2],-2)
  mu.mean.ex ~ dnorm(Prior.mu.mean.ex[1],mu.prec.ex)
  mu1.ex ~ dnorm(mu.mean.ex,tau.time.prec)
  mu.ex[1] <- mu1.ex
  prec.rho.ex <- pow(Prior.rho.ex[2],-2)
  
  #(equation (8) in paper)
  for(t in 2:Nint){
    mu.dlm.ex[t] <- mu.ex[t-1] + rho.ex[t-1]
    #variance of mu: discount factor X tau.time.prec 
    mu.time.prec.ex[t] <- tau.time.prec/w.ex
    mu.ex[t] ~ dnorm(mu.dlm.ex[t], mu.time.prec.ex[t])
    rho.ex[t-1] ~ dnorm(Prior.rho.ex[1], prec.rho.ex)
  }
  w.ex ~dunif(w1,w2)
  #NEX structure for mu: unrelated structure 
  for( s in 1:Nstudies){ 
    for(t in 1:Nint) {
      prior.mu.prec.nex[s,t] <- pow(prior.mu.sd.nex[s,t],-2)
      mu.nex[s,t] ~ dnorm(prior.mu.mean.nex[s,t],prior.mu.prec.nex[s,t])
    }
  }
  #baseline hazard lambda_jk, equation (5) in paper
  # hazard-base: for each study (+ population mean + prediction) and time-interval t, no covariates
  for ( s in 1:Nstudies) {
    for ( t in 1:Nint) {
      #random effect
      RE[s,t] ~ dnorm(0,tau.study.prec[t])
      #ex-nex indicator
      Z[s,t] ~ dbin(p.exch[s,t],1)
      
      #For ex
      #if only one study, then no RE
      log.hazard.base.ex[s,t] <- mu.ex[t] + step(Nstudies-1.5)*RE[s,t]
      
      #For Nexch
      log.hazard.base.nex[s,t] <- mu.nex[s,t]
      log.hazard.base[s,t] <- Z[s,t]*log.hazard.base.ex[s,t]+(1-Z[s,t])*log.hazard.base.nex[s,t]
      hazard.base[s,t] <- exp(log.hazard.base[s,t])
    }
  }
  # likelihood: pick hazards according to study and time-invervals (int.low to int.high) for each 
  #observation j
  # note: hazard is per unit time, not depending on length of interval t
  for(j in 1:Nobs) {
    for ( t in 1:Nint) {
      # log.hazard for all time intervals
      log(hazard.obs[j,t]) <- log.hazard.base[study[j],t] + inprod(X[j,1:Ncov],beta[1:Ncov,1])
    }
    
    
    #Poisson likelihood
    #alpha is Poisson rate
    #this inner product term doesn't really do anything... It is just lambda*E
    alpha[j] <- (inprod(hazard.obs[j,int.low[j]:int.high[j]],int.length[int.low[j]:int.high[j]])
                 /sum(int.length[int.low[j]:int.high[j]]))*exp.time[j]
    n.events[j] ~ dpois(alpha[j])
    n.events.pred[j] ~ dpois(alpha[j])
  }
  # outputs of interest (from covariate patterns in Xout)
  #Nout > 1 if there are covaraites
  #survival probabilities
  for ( h in 1:Nout) {
    for ( s in 1:Nstudies) {
      # mean (index = Nstudies)
      # surv: pattern x study x time
      surv1[ h,s,1] <- 1
      for ( t in 1:Nint) {
        #covariate effect is on the log-hazard rate scale
        #this is hazard, not hazard.base
        log(hazard[h,s,t]) <- log.hazard.base[s,t] + inprod(Xout[h,1:Ncov],beta[1:Ncov,1])
        log.hazard[h,s,t] <- log(hazard.base[s,t])
        surv1[h,s,t+1] <- surv1[h,s,t]*exp(-hazard[h,s,t]*int.length[t])
        surv[h,s,t] <- surv1[h,s,t+1]
      }
    }
  }
  ######### Prediction##########
  for ( h in 1:Nout) {
    for ( t in 1:Nint) {
      hazard.pred[h,t] <- hazard[h,Nstudies,t]
      log.hazard.pred[h,t] <- log(hazard[h,Nstudies,t]+pow(10,-6))
      surv.pred[h,t] <- surv[h,Nstudies,t]
    }
  }
  for (j in 1:Ncov) {
    beta.ge.cutoffs[j] <- step(beta[j,1]-beta.cutoffs[j,1])
  }
}  
