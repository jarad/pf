# Set graphics and data path
gpath = "../graphs/"
dpath = "../data/"

# Load simulated data
load(paste(dpath,"sim-orig.rdata",sep=""))

pf_coverage.params <- function(filt, prior, resamp, transform, delta, seed)
{
  n = c(100, 1000, 10000, 20000, 40000)
  n.sim = 1:20
  if(prior == "lognormal")
  {
    load(paste(dpath,"PF-quant-",filt,"-",prior,"-",resamp,"-100-",transform,"-",delta,"-",seed,"-1.rdata",sep=""))
  } else {
    load(paste(dpath,"PF-quant-uniform-systematic-",filt,"-100-",transform,"-",seed,"-1.rdata",sep=""))
  }
  tt = dim(pf.quant.out$state.quant)[1]
  n.params = dim(pf.quant.out$theta.quant)[2]
  covered = array(NA, dim=c(length(n.sim),length(n),tt,n.params))
  for(i in 1:length(n)){
  for(j in 1:length(n.sim))
  {
    if(prior == "lognormal")
    {
      load(paste(dpath,"PF-quant-",filt,"-",prior,"-",resamp,"-",n[i],"-",transform,"-",delta,"-",seed,"-",n.sim[j],".rdata",sep=""))
    } else {
      load(paste(dpath,"PF-quant-uniform-systematic-",filt,"-",n[i],"-",transform,"-",seed,"-",n.sim[j],".rdata",sep=""))
    }
    for(k in 1:n.params) covered[j,i,,k] = pf.quant.out$theta.quant[,k,2] < mysim$true.params$theta[k] & pf.quant.out$theta.quant[,k,3] > mysim$true.params$theta[k]
  }}
  return(apply(covered, c(2,4), mean))
}

pf_coverage.states <- function(filt, prior, resamp, transform, delta, seed)
{
  n = c(100, 1000, 10000, 20000, 40000)
  n.sim = 1:20
  if(prior == "lognormal")
  {
    load(paste(dpath,"PF-quant-",filt,"-",prior,"-",resamp,"-100-",transform,"-",delta,"-",seed,"-1.rdata",sep=""))
  } else {
    load(paste(dpath,"PF-quant-uniform-systematic-",filt,"-100-",transform,"-",seed,"-1.rdata",sep=""))
  }
  tt = dim(pf.quant.out$state.quant)[1]
  n.states = 2
  covered = array(NA, dim=c(length(n.sim),length(n),tt,n.states))
  for(i in 1:length(n)){
  for(j in 1:length(n.sim))
  {
    if(prior == "lognormal")
    {
      load(paste(dpath,"PF-quant-",filt,"-",prior,"-",resamp,"-",n[i],"-",transform,"-",delta,"-",seed,"-",n.sim[j],".rdata",sep=""))
    } else {
      load(paste(dpath,"PF-quant-uniform-systematic-",filt,"-",n[i],"-",transform,"-",seed,"-",n.sim[j],".rdata",sep=""))
    }
    for(k in 1:n.states) covered[j,i,,k] = pf.quant.out$state.quant[,k,2] < mysim$sim[[n.sim[j]]]$x[k,] & pf.quant.out$state.quant[,k,3] > mysim$sim[[n.sim[j]]]$x[k,]
  }}
  return(apply(covered, c(2,4), mean))
}

# Calculate average coverage probabilities over all particle filter runs
data1 = expand.grid(filt=c("BF","APF","KD"), prior = "uniform", resamp = "systematic", transform = "logit", delta = 0.99, seed = 61, stringsAsFactors=FALSE)
data2 = expand.grid(filt = "KD", prior = "lognormal", resamp = c("multinomial","residual","stratified","systematic"), transform="orig", delta = .99, seed = 61, stringsAsFactors=FALSE)
data3 = expand.grid(filt = "KD", prior = "lognormal", resamp = "stratified", transform="orig", delta = c(.9,.95,.96,.97,.98), seed = 61, stringsAsFactors=FALSE)
mydata = rbind(data1, data2, data3)
require(plyr)
out.params = mlply(.data = mydata, .fun = pf_coverage.params)
out.states = mlply(.data = mydata, .fun = pf_coverage.states)

# Compare coverage probabilities versus delta
tab.delta.params = laply(out.params[c(8:12,6)], function(x) x)
tab.delta.states = laply(out.states[c(8:12,6)], function(x) x)
dimnames(tab.delta.params) = list(delta=c(.9,seq(.95,.99,.01)), J=c(100,1000,10000,20000,40000), param = c("beta","gamma","nu"))
dimnames(tab.delta.states) = list(delta=c(.9,seq(.95,.99,.01)), J=c(100,1000,10000,20000,40000), state = c("s","i"))

# Compare coverage probabilities versus resampling
tab.resamp.params = laply(out.params[4:7], function(x) x)
tab.resamp.states = laply(out.states[4:7], function(x) x)
dimnames(tab.resamp.params) = list(resamp=c("multinomial","residual","stratified","systematic"), J=c(100,1000,10000,20000,40000), param = c("beta","gamma","nu"))
dimnames(tab.resamp.states) = list(resamp=c("multinomial","residual","stratified","systematic"), J=c(100,1000,10000,20000,40000), state = c("s","i"))

