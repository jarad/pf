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

# Run pmcmc
rw.sd = c(0.005, 0.001, 0.01); names(rw.sd) = c('beta','gamma','nu')
sir_pmcmc <- function(n.chain, niter, np, y.max, cont = FALSE, n.cont = 1)
{
  if(cont)
  {
    load(paste(dpath,"sir_pmcmc_test-",paste(n.chain,niter,np,y.max,sep="-"),".rdata",sep=""))
    time = system.time(out <- continue(out, Nmcmc = n.cont, max.fail = Inf))
    file = paste(dpath,"sir_pmcmc_test-",paste(n.chain,niter+n.cont,np,y.max,sep="-"),".rdata",sep="")
    save(out, file = file)
    cat(file,time[1:3],"\n")
  } else {
    # Create pomp object
    y = y[,1:y.max]
    sir.pomp <- pomp(data = y, times = 1:dim(y)[2], t0 = 0, P=P, b=b, varsigma=varsigma, sigma=sigma, eta = eta, rprocess=rprocess, rmeasure=rmeasure, dmeasure=dmeasure, skeleton = skeleton, skeleton.type = "map", initializer=initializer, rprior=rprior.pomp, dprior=dprior)

    # Run pmcmc
    pars.init = mysims[[1]]$true.params$theta
    names(pars.init) = c('beta','gamma','nu')
    time = system.time(out <- pmcmc(sir.pomp, Nmcmc = niter, start = pars.init, pars = c('beta','gamma','nu'), rw.sd = rw.sd, Np = np, max.fail = Inf))
    file = paste(dpath,"sir_pmcmc_test-",paste(n.chain,niter,np,y.max,sep="-"),".rdata",sep="")
    save(out, file = file)
    cat(file,time[1:3],"\n")
  }
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(n.chain=1:3, niter=20000, np=100, y.max = c(30,60,90,125,65,70,75,80,85,95,100,105,110,115,120,5,10,15,20,25,35,40,45,50,55), cont = TRUE, n.cont = 10000)
m_ply(.data = mydata, .fun = sir_pmcmc, .parallel = TRUE)
