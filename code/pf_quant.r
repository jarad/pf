source("pf_functions.r")

# Set data path
dpath = "../data/"

pf.quant = function(n, mod, filt, resamp, prior, nonunif = "", thresh = "")
{
  # Load data
  load(paste(dpath,"PF",mod,"-",filt,"-",prior,"-",resamp,"-",n,nonunif,thresh,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta

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
  save(pf.quant.out, file=paste(dpath,"PF-quant",mod,"-",filt,"-",prior,"-",resamp,"-",n,nonunif,thresh,".rdata",sep=""))
}

# Calculate quantiles for particle filters run
require(plyr)
mydata = expand.grid(n = c(100, 1000, 10000, 20000, 40000), mod = "-2", filt = "KD", resamp = "stratified", prior = "normal", nonunif = "-ess", thresh = paste("-",100*.8,sep=""), stringsAsFactors=FALSE)
m_ply(mydata,pf.quant)

# Clear objects
rm(list=ls(all=TRUE))