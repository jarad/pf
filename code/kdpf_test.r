source("sir_functions.r")
source("sir_mcmc_functions.r")
source("pf_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# Set pf params
n = 20000
delta = 0.99
progress = "FALSE"
resamp = "stratified"
n.sim = 1

# Load simulated data
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysims[[n.sim]]$sim$y
P = mysims[[n.sim]]$true.params$P
b = mysims[[n.sim]]$true.params$b
varsigma = mysims[[n.sim]]$true.params$varsigma
sigma = mysims[[n.sim]]$true.params$sigma
eta = mysims[[n.sim]]$true.params$eta  
  
# Define functions to transform theta to original scale
ftheta = function(theta,param=1) exp(theta)
  
# Define functions to sample prior draws of fixed parameters
rtheta <- function()
{
  theta <- rep(NA, 3)
  log.params <- find.mu.sigma(c(1.5, .09, .95), c(3, .143, 1.3))
  theta[2:3] <- exp(rnorm(2, log.params[[1]][2:3], log.params[[2]][2:3]))
  theta[1] <- theta[2]*exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
  return(log(theta))
}
 
# Function to run kd_pf
kd.pf <- function(seed)
{
  # Set seed
  set.seed(seed)

  # Run KDPF
  mydllik = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
  myrevo = function(x, theta) revo(x, ftheta(theta), P)
  pstate = function(x, theta) revo(x, ftheta(theta), P, FALSE)
  myrprior = function() rprior(rtheta)
  source("kd_pf.r")
  time = system.time(out <- kd_pf(y, mydllik, pstate, myrevo, myrprior, n, delta, progress, method=resamp, nonuniformity = "ess", threshold = 0.8, log=F))

  # Save KDPF
  pf.out = list(out=out,ftheta=ftheta)
  file = paste(dpath,"PF-",n.sim,"-",n,"-KD-",resamp,"-orig-log-",delta,"-",seed,".rdata",sep="")
  cat(file,time[1:3],"\n")
  save(pf.out, file=file)
  
  # Calculate 95% credible intervals
  probs = c(0.5, 0.25, 0.75, 0.025, 0.975, 0.05, 0.95)
  tt = dim(out$state)[3]
  states = array(NA,dim=c(3,n,tt))
  states[1:2,,] = out$state[1:2,,]
  for(i in 1:tt) states[3,,i] = 1 - apply(out$state[1:2,,i],2,sum)
  state.quant = pf.quantile(states, out$weight, function(x,param=1) x, probs)
  theta.quant = pf.quantile(out$theta, out$weight, ftheta, probs)

  # Save data
  pf.quant.out = list(state.quant=state.quant,theta.quant=theta.quant,probs=probs)
  file = paste(dpath,"PF-quant-",n.sim,"-",n,"-KD-",resamp,"-orig-log-",delta,"-",seed,".rdata",sep="")
  print(file)
  save(pf.quant.out, file=file)
}

require(plyr)
require(doMC)
registerDoMC()
mydata = data.frame(seed=62:80,stringsAsFactors = FALSE)
m_ply(.data = mydata, .fun = kd.pf, .parallel = TRUE)