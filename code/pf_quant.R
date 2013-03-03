source("pf_functions.R")

# Set data path
dpath = "../data/"

pf.quant = function(n, filt, resamp, prior, ...)
{
  # Load data
  load(paste(dpath,"PF-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta
  p = pf.out$p

  # Calculate 2.5%, 50%, and 97.5% quantiles of states over time
  probs = c(.5,.025,.975)
  state.quant = pf.quantile(out$state[1:2,,], out$weight, function(x,param=1) x, probs)

  # Calculate 2.5%, 50%, and 97.5% quantiles of parameters over time
  if("theta" %in% names(out))
  {
    theta.quant = pf.quantile(out$theta, out$weight, ftheta, probs)
  } else if(dim(out$state)[1] > 2){
    theta.quant = pf.quantile(out$state[3:(dim(out$state)[1]),,], out$weight, ftheta, probs)
  } else { theta.quant = NULL}

  # Save data
  pf.quant.out = list(state.quant=state.quant,theta.quant=theta.quant,p=p,probs=probs)
  save(pf.quant.out, file=paste(dpath,"PF-quant-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
}

# Calculate quantiles for particle filters run
require(plyr)
load(paste(dpath,"filt_instruc.rdata",sep=""))
m_ply(mydata,pf.quant)

# Clear objects
rm(list=ls(all=TRUE))