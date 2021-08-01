##Description of arguments of the function
#Nobs                       =  Total number of data points, capped at Nint*Nstudies, no need to enter 0/0 data point(?)
#study                       =  Study indicator
#Nstudies                 =  Number of studies
#Nint                         =  Number of intervals
#Ncov                        =  Number of covariates
#int.low, int.high     =  Interval indicator
#int.length                =  Length of interval
#n.events                  =  Number of events at each interval
#exp.time                  =  Exposure time for each interval
#X                                =  Important covariates
#Prior.mu.mean.ex  =  Mean and standard deviation for normal priors for exchangeability for the first interval
#Prior.rho.ex              =  Mean and standard deviation for normal priors for random gradient of 1st order NDLM
#w1, w2                       =  Upper and lower bound for uniform prior of the smoothing factor for 1st order NDLM model for exchangeability 
#prior.mu.mean.nex  =  Mean of normal prior for non-exchangeability
#prior.mu.sd.nex         =  Standard deviation for normal prior of non-exchangeability
#p.exch                          =  Prior probability of exchangeability
#Prior.beta                   =  Mean and standard deviation for normal priors of regression coefficients  
#beta.cutoffs               =  Cut-off for treatment effect
#Prior.tau.study          =  location and scale parameter of half-normal prior for between trial heterogeneity, identical for all time intervals
#Prior.tau.time            =  Mean and standard deviation for log-normal prior for variance compoment of 1st order NDLM
#MAP.Prior                   =  If TRUE: derives MAP prior 
#pars                             =  Parameters to keep in each MCMC run
#bugs.directory           =  Directory path where WinBUGS14.exe file resides
#R.seed                         =  Seed to generate initial value in R (requires for reprducibility) 
#bugs.seed                   =  WinBUGS seed (requires for reprducibility)

MAC.Surv.anal <- function(Nobs                = NULL,
                          study                = NULL,
                          Nstudies           = NULL,
                          Nint                   = NULL,
                          Ncov                  = 1,
                          Nout                  = 1,
                          int.low              = NULL,
                          int.high             = NULL,
                          int.length         = NULL,
                          n.events           = NULL,
                          exp.time          = NULL,
                          X                       = NULL,
                          Xout                  = matrix(0,1,1),
                          Prior.mu.mean.ex   = NULL, 
                          Prior.rho.ex              = NULL,
                          w1                              = 0,
                          w2                              = 1,
                          prior.mu.mean.nex = NULL,
                          prior.mu.sd.nex       = NULL,
                          p.exch                       = NULL,
                          Prior.beta                 = NULL, 
                          beta.cutoffs             = NULL,
                          Prior.tau.study        = NULL, 
                          Prior.tau.time          = NULL,
                          MAP.prior                = FALSE,
                          pars                           = c("tau.study","log.hazard","mu.ex", "hazard"),
                          bugs.directory          = "C:/Users/hongzhang/Apps/WinBUGS14",
                          R.seed                        = 10,
                          bugs.seed                  = 12
)
{  
  set.seed(R.seed)
  beta.cutoffs <- cbind(NULL,beta.cutoffs)
  if(MAP.prior){Nstudies <- Nstudies+1}
  #WinBUGS model
  model <- MAC.Surv.WB
  
  #Data for JAGS format
  data     = list("Nstudies", "Nint", "Ncov", "Nobs", "Nout","study","int.low","int.high","int.length",
                  "n.events","exp.time",
                  "X","Xout",
                  "Prior.mu.mean.ex","Prior.rho.ex",
                  "prior.mu.mean.nex","prior.mu.sd.nex",
                  "p.exch","w1","w2",
                  "Prior.beta",
                  "beta.cutoffs",
                  "Prior.tau.study","Prior.tau.time"
  )
  #Initial values
  hazard0 = (sum(n.events)+0.5)/sum(exp.time)
  initsfun = function(i)
    list(
      mu1.ex = rnorm(1,log(hazard0),0.25),
      tau.study = rgamma(12,1,1),
      tau.time  = rgamma(1,1,1),
      mu.mean.ex = rnorm(1,log(hazard0),0.1),
      rho.ex = rnorm(Nint-1,0,0.05),
      w.ex = runif(1,0,1),
      beta=cbind(NULL,rnorm(Ncov,0,1))
    )
  inits <- lapply(rep(1,3),initsfun)
  
  
  #WinBUGS run
  fit = bugs(
    data=data,
    inits=inits,
    par=pars,
    model=model,
    n.chains=3,n.burnin=8000,n.iter=16000,n.thin=1,
    bugs.directory= bugs.directory,
    bugs.seed= bugs.seed,
    DIC= TRUE,
    debug= F
  )
  fit$sims.matrix = NULL
  fit$sims.array = NULL
  #WinBUGS summary
  summary <- fit$summary 
  R2WB <- fit
  
  output <- list(summary=summary,R2WB =R2WB)
  
  return(output)
  
}


MAC.Surv.anal_jags <- function(Nobs                = NULL,
                          study                = NULL,
                          Nstudies           = NULL,
                          Nint                   = NULL,
                          Ncov                  = 1,
                          Nout                  = 1,
                          int.low              = NULL,
                          int.high             = NULL,
                          int.length         = NULL,
                          n.events           = NULL,
                          exp.time          = NULL,
                          X                       = NULL,
                          Xout                  = matrix(0,1,1),
                          Prior.mu.mean.ex   = NULL, 
                          Prior.rho.ex              = NULL,
                          w1                              = 0,
                          w2                              = 1,
                          prior.mu.mean.nex = NULL,
                          prior.mu.sd.nex       = NULL,
                          p.exch                       = NULL,
                          Prior.beta                 = NULL, 
                          beta.cutoffs             = NULL,
                          Prior.tau.study        = NULL, 
                          Prior.tau.time          = NULL,
                          MAP.prior                = FALSE,
                          pars                           = c("tau.study","log.hazard","mu.ex", "hazard"),
                          bugs.directory          = "C:/Program Files/WinBUGS14",
                          R.seed                        = 10,
                          bugs.seed                  = 12
)
{  
  # browser()
  
  set.seed(R.seed)
  beta.cutoffs <- cbind(NULL,beta.cutoffs)
  if(MAP.prior){Nstudies <- Nstudies+1}
  #WinBUGS model
  model <- MAC.Surv.WB
  
  #Data for JAGS format
  data     = list("Nstudies", "Nint", "Ncov", "Nobs", "Nout","study","int.low","int.high","int.length",
                  "n.events","exp.time",
                  "X","Xout",
                  "Prior.mu.mean.ex","Prior.rho.ex",
                  "prior.mu.mean.nex","prior.mu.sd.nex",
                  "p.exch","w1","w2",
                  "Prior.beta",
                  "beta.cutoffs",
                  "Prior.tau.study","Prior.tau.time"
  )
  #Initial values
  hazard0 = (sum(n.events)+0.5)/sum(exp.time)
  initsfun = function(i)
    list(
      mu1.ex = rnorm(1,log(hazard0),0.25),
      tau.study = rgamma(12,1,1),
      tau.time  = rgamma(1,1,1),
      mu.mean.ex = rnorm(1,log(hazard0),0.1),
      rho.ex = rnorm(Nint-1,0,0.05),
      w.ex = runif(1,0,1),
      beta=cbind(NULL,rnorm(Ncov,0,1))
    )
  inits <- lapply(rep(1,3),initsfun)
  
  
  #WinBUGS run
  # fit = bugs(
  #   data=data,
  #   inits=inits,
  #   par=pars,
  #   model=model,
  #   n.chains=3,n.burnin=8000,n.iter=16000,n.thin=1,
  #   bugs.directory= bugs.directory,
  #   bugs.seed= bugs.seed,
  #   DIC= TRUE,
  #   debug= FALSE
  # )
  #JAGS run
  fit <- jags(
    data=data,
    inits=inits,
    parameters.to.save = pars,
    model.file = model,
    n.chains=3,n.burnin=8000,n.iter=16000,n.thin=1,
    jags.seed= bugs.seed,
    DIC= TRUE
  )
  
  return(fit)
  # fit$sims.matrix = NULL
  # fit$sims.array = NULL
  #WinBUGS summary
  # summary <- fit$summary 
  # R2WB <- fit
  # 
  # output <- list(summary=summary,R2WB =R2WB)
  # 
  # return(output)
  
}
