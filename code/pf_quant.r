source("pf_functions.r")
source("sir_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

pf.quant = function(data = "", n.sim, n, filt, resamp, prior, transform, delta, seed)
{
  # Load data
  load(paste(dpath,"PF-",data,n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta
  
  # Which quantiles to calculate
  probs = c(0.5, 0.25, 0.75, 0.025, 0.975, 0.05, 0.95)
  
  # Calculate quantiles of states over time
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
  file = paste(dpath,"PF-",data,"quant-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep="")
  print(file)
  save(pf.quant.out, file=file)
}

# Create data frame and use plyr to run particle filters
data1 = expand.grid(data = "", n.sim = 21:40, n = c(100, 1000, 10000, 20000), filt = c("BF","APF","KD"), resamp = "systematic", prior = "unif", transform = "logit", delta = 0.99, seed = 61, stringsAsFactors=FALSE)
data2 = expand.grid(data = "", n.sim = 21:40, n = c(100, 1000, 10000, 20000), filt = "KD", resamp = c("multinomial","residual","stratified","systematic"), prior="orig", transform="log", delta = .99, seed = 61, stringsAsFactors = FALSE)
data3 = expand.grid(data = "", n.sim = 21:40, n = c(100, 1000, 10000, 20000), filt = "KD", resamp = "stratified", prior="orig", transform="log", delta = c(0.9,0.95,0.96,0.97,0.98), seed = 61, stringsAsFactors = FALSE)
data4 = expand.grid(data = "", n.sim = 21:40, n = c(100, 1000, 10000, 20000), filt = "KD", resamp = "stratified", prior = "disp", transform = "log", delta = 0.99, seed = 61, stringsAsFactors=FALSE)
data5 = expand.grid(data = "ext-", n.sim = 1:20, n = c(10000, 20000, 40000), filt = "KD", resamp = "stratified", prior = "orig", transform = "log", delta = 0.99, seed = 61, stringsAsFactors = FALSE)
data6 = expand.grid(data = "ext.orig-", n.sim = 1:20, n = c(10000, 20000, 40000), filt = "KD", resamp = "stratified", prior = "orig", transform = "log", delta = 0.99, seed = 61, stringsAsFactors = FALSE)
data7 = expand.grid(data = "ext.orig-", n.sim = 1, n = 60000, filt = "KD", resamp = "stratified", prior = "orig", transform = "log", delta = 0.99, seed = 61, stringsAsFactors = FALSE)
mydata = mydata = rbind(data1,data2,data3,data4,data5,data5,data7)

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf.quant, .parallel = TRUE)