source("pf_functions.r")

# Set data path
dpath = "../data/"

# Which model
mod = ""

# Which data frame to use
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("BF","APF","KD"), resamp = c("multinomial","residual","stratified","systematic"), prior = c("normal","uniform"), nonuniformity="ess", threshold=0.8, stringsAsFactors=FALSE)

pf.quant = function(n, filt, resamp, prior, ...)
{
  # Load data
  load(paste(dpath,"PF",mod,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta

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
  pf.quant.out = list(state.quant=state.quant,theta.quant=theta.quant,probs=probs)
  save(pf.quant.out, file=paste(dpath,"PF-quant",mod,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
}

# Calculate quantiles for particle filters run
require(plyr)
m_ply(mydata,pf.quant)

# Clear objects
rm(list=ls(all=TRUE))