source("pf_functions.r")
source("sir_functions.r")

# Set data path
dpath = "../data/"

pf.quant = function(n.sim, n, filt, resamp, prior, transform, delta, seed)
{
  # Load data
  load(paste(dpath,"PF-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta
  
  # Calculate 2.5%, 50%, and 97.5% quantiles of states over time
  probs = c(.5, 0.25, 0.975)
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
  file = paste(dpath,"PF-quant-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep="")
  print(file)
  save(pf.quant.out, file=file)
}

# Create data frame and use plyr to run particle filters
data1 = expand.grid(n.sim = 1, n = 100, filt = c("BF","APF","KD"), resamp = "systematic", prior = "unif", transform = "logit", delta = 0.99, seed = 61, stringsAsFactors=FALSE)
data2 = expand.grid(n.sim = 1, n = c(100), filt = "KD", resamp = c("multinomial","residual","stratified","systematic"), prior="orig", transform="log", delta = .99, seed = 61, stringsAsFactors = FALSE)
data3 = expand.grid(n.sim = 1, n = c(100), filt = "KD", resamp = "stratified", prior="orig", transform="log", delta = c(0.9,0.95,0.96,0.97,0.98), seed = 61, stringsAsFactors = FALSE)
data4 = expand.grid(n.sim = 1, n = c(100), filt = "KD", resamp = "stratified", prior = "disp", transform = "log", delta = 0.99, seed = 61, stringsAsFactors=FALSE)
mydata = rbind(data1, data2, data3, data4)

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf.quant, .parallel = TRUE)
