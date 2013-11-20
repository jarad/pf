source("sir_functions.r")
source("pf_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# pf - function to run particle filter given n number of particles, filt = "BF", "APF", or "KD" for which filter to run, resamp = "multinomial", "residual", "stratified", or "systematic" for which resampling method to use, and prior = "normal" or "uniform" for which prior to use on unknown parameters
# Returns nothing; saves .rdata data file
pf <- function(n, filt, resamp, prior, mod.data, mod.fit, nonunif = "ess", thresh = 0.8, seed = rnorm(1), progress = TRUE)
{
  # Load simulated data
  load(paste(dpath,"sim-",mod.data,".rdata",sep=""))
  y = mysim$sim$y
  P = mysim$true.params$P
  b = mysim$true.params$b
  varsigma = mysim$true.params$varsigma
  sigma = mysim$true.params$sigma
  eta = mysim$true.params$eta
  nu = mysim$true.params$theta[3]

  # Create functions to draw fixed parameters and map theta to original scale
  if(mod.fit == "orig")
  {
    if(prior == "uniform")
    { 
      thetal = c(0.1400, 0.0900, 0.9500)
      thetau = c(0.5000, 0.1430, 1.3000)
#      rtheta = function() u2theta(runif(3,thetal,thetau),thetal,thetau)
      ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
      rtheta = function() u2theta(runif(2,thetal[1:2],thetau[1:2]),thetal[1:2],thetau[1:2])
    } else if(prior == "semi-uniform") {
      thetal = c(0.1400, 0.0900, 0.9500)
      thetau = c(0.5000, 0.1430, 1.3000)
      rtheta = function() runif(3,thetal,thetau)
      ftheta = function(theta,param=1) theta
    } else if(prior == "lognormal") {
      theta.mean = c(-1.3296, -2.1764, 0.1055)
      theta.sd = c(.3248, .1183, .0800)
#      rtheta = function() rnorm(3,theta.mean,theta.sd)
      ftheta = function(theta,param=1) exp(theta)
      rtheta = function() rnorm(2, theta.mean[1:2], theta.sd[1:2])
    } else if(prior == "semi-lognormal") {
      theta.mean = c(-1.3296, -2.1764, 0.1055)
      theta.sd = c(.3248, .1183, .0800)
      rtheta = function() exp(rnorm(3,theta.mean,theta.sd))
      ftheta = function(theta,param=1) theta
    } else { stop("prior must be either 'uniform', 'semi-uniform', 'lognormal', or 'semi-lognormal'")}
  } else if(mod.fit == "ext") {
    if(prior == "uniform")
    { 
      thetal = c(0.1400, 0.0900, 0.9500, 0.0001, 0.85, 0.0001, 0)
      thetau = c(0.5000, 0.1430, 1.3000, 0.5000, 1.15, 0.0100, 5)
      rtheta = function() u2theta(runif(7,thetal,thetau),thetal,thetau)
      ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
    } else if(prior == "semi-uniform") {
      thetal = c(0.1400, 0.0900, 0.9500, 0.0001, 0.85, 0.0001, 0)
      thetau = c(0.5000, 0.1430, 1.3000, 0.5000, 1.15, 0.0100, 5)
      rtheta = function() runif(7,thetal,thetau)
      ftheta = function(theta,param=1) theta
    } else if(prior == "lognormal") {
      theta.mean = c(-1.3296, -2.1764, 0.1055, -1.609, -0.0114, -7.0516, 2.5)
      theta.sd = c(.3248, .1183, .0800, .3536, .0771, .2803, 1)
      rtheta = function() rnorm(7, theta.mean, theta.sd)
      ftheta = function(theta,param=1) if(param == 7) theta else exp(theta)
    } else if(prior == "semi-lognormal") {
      theta.mean = c(-1.3296, -2.1764, 0.1055, -1.609, -0.0114, -7.0516, 2.5)
      theta.sd = c(.3248, .1183, .0800, .3536, .0771, .2803, 1)
      rtheta = function()
      {
        prior.draw = rnorm(7, theta.mean, theta.sd)
        return(c(exp(prior.draw[1:6]),prior.draw[7]))
      }
      ftheta = function(theta,param=1) theta
    } else { stop("prior must be either 'uniform', 'semi-uniform', 'lognormal', or 'semi-lognormal'")}
  } else { stop("mod must be either 'orig' or 'ext'")}

  # Set seed
  set.seed(seed)

  # Run one of the particle filters
  if(mod.fit == "orig")
  {
    if(filt == "BF")
    {
      # Run bootstrap filter, 
      dllik_bf = function(y, x) dllik(y, x, b, varsigma, sigma, eta)
      revo_bf = function(x) revo(x, ftheta(x[3:5],1:3), P)
      rprior_bf = function()
      { 
        myprior = rprior(rtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("bf.r")
      out = bf(y, dllik_bf, revo_bf, rprior_bf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else if(filt == "APF"){
      # Run auxiliary particle filter
      dllik_apf = function(y, x) dllik(y, x, b, varsigma, sigma, eta)
      pstate_apf = function(x) revo(x, ftheta(x[3:5],1:3), P, FALSE)
      revo_apf = function(x) revo(x, ftheta(x[3:5],1:3), P)
      rprior_apf = function()
      { 
        myprior = rprior(rtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("apf.r")
      out = apf(y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else if(filt == "KD"){
      # Run kernel density particle filter
      dllik_kd = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
#      pstate_kd = function(x, theta) revo(x, ftheta(theta,1:3), P, FALSE)
#      revo_kd = function(x, theta) revo(x, ftheta(theta,1:3), P)
      pstate_kd = function(x, theta) revo(x, c(ftheta(theta,1:2),nu), P, FALSE)
      revo_kd = function(x, theta) revo(x, c(ftheta(theta,1:2),nu), P)
      rprior_kd = function() rprior(rtheta)
      source("kd_pf.r")
#      sink(paste("kd_pf-",n,".txt",sep=""))
      out = kd_pf(y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
#      sink()
    } else if(filt == "RM"){
      source("sir_mcmc_functions.r")
      dllik_rm = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
      revo_rm = function(x, theta) revo(x, c(theta[1:2],nu), P)
      rtheta = function() exp(rnorm(2, theta.mean, theta.sd))
      rprior_rm = function() rprior(rtheta)
      rmove_rm = function(y, x, theta) rmove(y, x, c(theta[1:2], nu), list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P)) 
      source("rm_pf.r")
    
      # Run resample-move particle filter
#      dllik_rm = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
#      revo_rm = function(x, theta) revo(x, ftheta(theta,1:3), P)
#      revo_rm = function(x, theta) revo(x, c(ftheta(theta,1:2),nu), P)
#      rprior_rm = function() rprior(rtheta)
#      rmove_rm = function(y, x, theta) rmove(y, x, ftheta(theta,1:3), b, varsigma, sigma, eta, P) 
#      rmove_rm = function(y, x, theta) rmove(y, x, c(ftheta(theta,1:2), nu), b, varsigma, sigma, eta, P) 
#      source("rm_pf.r")
#      sink(paste("rm_pf-",n,".txt",sep=""))
      out = rm_pf(y, dllik_rm, revo_rm, rprior_rm, rmove_rm, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
#      sink()
    } else{ stop("filt must be one of 'BF', 'APF', 'KD', or 'RM'") } 
  } else if(mod.fit == "ext") {
    if(filt == "BF")
    {
      # Run bootstrap filter
      dllik_bf = function(y, x) dllik(y, x[1:2], ftheta(x[6],4), ftheta(x[7],5), ftheta(x[8],6), ftheta(x[9],7))
      revo_bf = function(x) revo(x, c(ftheta(x[3],4),ftheta(x[4],5),ftheta(x[5],6)), P)
      rprior_bf = function()
      {
        myprior = rprior(rtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("bf.r")
      out = bf(y, dllik_bf, revo_bf, rprior_bf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else if(filt == "APF"){
      # Run auxiliary particle filter
      dllik_apf = function(y, x) dllik(y, x[1:2], ftheta(x[6],4), ftheta(x[7],5), ftheta(x[8],6), ftheta(x[9],7))
      pstate_apf = function(x) revo(x, c(ftheta(x[3],4),ftheta(x[4],5),ftheta(x[5],6)), P, FALSE)
      revo_apf = function(x) revo(x, c(ftheta(x[3],4),ftheta(x[4],5),ftheta(x[5],6)), P)
      rprior_apf = function()
      { 
        myprior = rprior(rtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("apf.r")
      out = apf(y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else {
      # Run kernel density particle filter
      dllik_kd = function(y, x, theta) dllik(y, x, ftheta(theta[4],4), ftheta(theta[5],5), ftheta(theta[6],6), ftheta(theta[7],7))
      pstate_kd = function(x, theta) revo(x, c(ftheta(theta[1],1),ftheta(theta[2],2),ftheta(theta[3],3)), P, FALSE)
      revo_kd = function(x, theta) revo(x, c(ftheta(theta[1],1),ftheta(theta[2],2),ftheta(theta[3],3)), P)
      rprior_kd = function() rprior(rtheta)
      source("kd_pf.r")
      out = kd_pf(y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    }
  } else { stop("mod must be either 'orig' or 'ext'")}

  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  save(pf.out, file=paste(dpath,"PF-",mod.data,"-",mod.fit,"-",filt,"-",prior,"-",resamp,"-",n,"-",nonunif,"-",100*thresh,".rdata",sep=""))
}

pf.quant = function(n, filt, resamp, prior, mod.data, mod.fit, nonunif = "ess", thresh = 0.8)
{
  # Load data
  load(paste(dpath,"PF-",mod.data,"-",mod.fit,"-",filt,"-",prior,"-",resamp,"-",n,"-",nonunif,"-",100*thresh,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta
  
  # If resample-move particle filter, obtain last element of list of states
  if(filt == "RM")
  {
    require(plyr)
    out.t = laply(out$state, function(x) x[,,dim(x)[3]])
    out$state = aperm(out.t, c(2,3,1))
    ftheta = function(x,param=1) x
  }

  # Calculate 2.5%, 50%, and 97.5% quantiles of states over time
  probs = c(.5,.025,.975)
  tt = dim(out$state)[3]
  states = array(NA,dim=c(3,n,tt))
  states[1:2,,] = out$state[1:2,,]
  for(i in 1:tt)
  {
    states[3,,i] = 1 - apply(out$state[1:2,,i],2,sum)
  }
  state.quant = pf.quantile(states, out$weight, function(x,param=1) x, probs)

  # Calculate 2.5%, 50%, and 97.5% quantiles of parameters over time
  if("theta" %in% names(out))
  {
    theta.quant = pf.quantile(out$theta, out$weight, ftheta, probs)
  } else if(dim(out$state)[1] > 2){
    theta.quant = pf.quantile(out$state[3:(dim(out$state)[1]),,], out$weight, ftheta, probs)
  } else { theta.quant = NULL}

  # Save data
  pf.quant.out = list(state.quant=state.quant,theta.quant=theta.quant,probs=probs)
  save(pf.quant.out, file=paste(dpath,"PF-quant-",mod.data,"-",mod.fit,"-",filt,"-",prior,"-",resamp,"-",n,"-",nonunif,"-",100*thresh,".rdata",sep=""))
}
