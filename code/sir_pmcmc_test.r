source("sir_pmcmc_functions.r")
require(pomp)

# Set data path
dpath = "../data/"

# Load simulated data and redefine data / known parameter values
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysims[[1]]$sim$y
P = mysims[[1]]$true.params$P
b = mysims[[1]]$true.params$b
varsigma = mysims[[1]]$true.params$varsigma
sigma = mysims[[1]]$true.params$sigma
eta = mysims[[1]]$true.params$eta

# Create pomp object
sir.pomp <- pomp(data = y, times = 1:dim(y)[2], t0 = 0, P=P, b=b, varsigma=varsigma, sigma=sigma, eta = eta, rprocess=rprocess, rmeasure=rmeasure, dmeasure=dmeasure, skeleton = skeleton, skeleton.type = "map", initializer=initializer, rprior=rprior.pomp, dprior=dprior)

# Run pmcmc
rw.sd = c(0.005, 0.001, 0.01); names(rw.sd) = c('beta','gamma','nu')
sir_pmcmc <- function(n.chain, niter, np, ymax)
{
  # Create pomp object
  sir.pomp <- pomp(data = y[,1:ymax], times = 1:dim(y)[2], t0 = 0, P=P, b=b, varsigma=varsigma, sigma=sigma, eta = eta, rprocess=rprocess, rmeasure=rmeasure, dmeasure=dmeasure, skeleton = skeleton, skeleton.type = "map", initializer=initializer, rprior=rprior.pomp, dprior=dprior)
  
  # Run pmcmc
  pars.init = mysims[[1]]$true.params$theta
  names(pars.init) = c('beta','gamma','nu')
  time = system.time(out <- pmcmc(sir.pomp, Nmcmc = niter, start = pars.init, pars = c('beta','gamma','nu'), rw.sd = rw.sd, Np = np, max.fail = Inf))
  file = paste(dpath,"sir_pmcmc_test-hsd-",paste(n.chain,niter,np,sep="-"),".rdata",sep="")
  save(out, file = file)
  print(time)
}
require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(n.chain=1:3, niter=20000, np=100, ymax = c(30,60,90,125))
m_ply(.data = mydata, .fun = sir_pmcmc, .parallel = TRUE)
