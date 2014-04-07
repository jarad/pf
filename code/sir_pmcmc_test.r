source("sir_functions.r")
myrprior = rprior
source("sir_mcmc_functions.r")

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

# Create functions for 'pomp' object
rprocess <- function(xstart, times, params, P, log.par = FALSE, ...)
{
  ns <- dim(xstart)[1]
  np <- dim(xstart)[2]
  tt <- length(times)
  if(log.par) params = exp(params)
  
  x <- array(NA, dim = c(ns, np, tt))
  x[,,1] = xstart
  for(i in 2:tt)
  {
    for(j in 1:np)
    {
      x[,j,i] = revo(x[,j,i-1], params, P)
    }
  }
  rownames(x) = c("s","i")
  return(x)
}

rmeasure <- function(x, t, params, b, varsigma, sigma, eta, ...)
{
  y <- robs(x,b,varsigma,sigma,eta)
  names(y) = 1:length(y)
  return(y)
}

dmeasure <- function(y, x, t, params, log = TRUE, b, varsigma, sigma, eta, ...)
{
  log.pdf = dllik(y, x, b, varsigma, sigma, eta)
  if(log) return(log.pdf) else return(exp(log.pdf))
}

skeleton <- function(x, t, params, P, log.par = FALSE, ...)
{
  if(log.par) params = exp(params)
  x.map = revo(x, params, P, random = FALSE)
  names(x.map) = c("s","i")
  return(x.map)
}

initializer <- function(params, t0, ...)
{
  i0 = rnorm(1,0.002,0.0005)
  while(i0 < 0 | i0 > 1) i0 = rnorm(1,0.002,0.0005)
  s0 = 1 - i0
  x.init = c(s0, i0)
  names(x.init) = c("s","i")
  return(x.init)
}

rprior <- function(params, log=FALSE, ...)
{
  theta <- rep(NA, 3)
  log.params <- find.mu.sigma(c(1.5, .09, .95), c(3, .143, 1.3))
  theta[2:3] <- exp(rnorm(2, log.params[[1]][2:3], log.params[[2]][2:3]))
  theta[1] <- theta[2]*exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
  names(theta) = c("beta","gamma","nu")
  if(!log) return(theta) else return(log(theta))
}

dprior <- function(params, log=FALSE, log.par = FALSE, ...)
{
  if(log.par) params = exp(params)
  log.prior <- dlbeta0gamma0(params[1], params[2]) + dlnu0(params[3])
  if(log) return(log.prior) else return(exp(log.prior))
}

par.trans = function(params, ...)
{
  log.params = log(params)
  names(log.params) = c('beta','gamma','nu')
  return(log.params)
}

par.inv.trans = function(params, ...)
{
  exp.params = exp(params)
  names(exp.params) = c('beta','gamma','nu')
  return(exp.params)
}

# Create pomp object
sir.pomp <- pomp(data = y, times = 1:dim(y)[2], t0 = 0, P=P, b=b, varsigma=varsigma, sigma=sigma, eta = eta, rprocess=rprocess, rmeasure=rmeasure, dmeasure=dmeasure, skeleton = skeleton, skeleton.type = "map", initializer=initializer, rprior=rprior, dprior=dprior)

# # Run pfilter
# np = 1000
# params = matrix(NA, nr=3, nc=np)
# for(i in 1:np) params[,i] = rprior(i)
# rownames(params) = c('beta','gamma','nu')
# sir.pf = pfilter(sir.pomp, params, np, filter.mean = TRUE, save.states=TRUE, save.params=TRUE, verbose=TRUE)
# 
# # Plot 95% credible intervals for beta
# lbeta = sapply(sir.pf$saved.params, function(x) quantile(x[1,],.025))
# ubeta = sapply(sir.pf$saved.params, function(x) quantile(x[1,],.975))
# ymin = min(lbeta, mysims[[1]]$true.params$theta[1])
# ymax = max(ubeta, mysims[[1]]$true.params$theta[1])
# plot(1:dim(y)[2], lbeta, ylim=c(ymin,ymax), type="l", col = 2)
# lines(1:dim(y)[2], ubeta, col = 2)
# abline(h=mysims[[1]]$true.params$theta[1], col = "gray70", lwd = 3)
# 
# # Liu and West filter
# rprocess.kd <- function(xstart, times, params, ...) rprocess(xstart, times, params, P, log.par = TRUE, ...)
# skeleton.kd <- function(x, t, params, ...) skeleton(x, t, params, P, log.par = TRUE, ...)
# rprior.kd <- function(params, ...) rprior(params, log = TRUE, ...)
# dprior.kd <- function(params, log=FALSE, ...) dprior(params, log=FALSE, log.par = TRUE, ...)
# sir.pomp.kd <- pomp(data = y, times = 1:dim(y)[2], t0 = 0, P=P, b=b, varsigma=varsigma, sigma=sigma, eta = eta, rprocess=rprocess.kd, rmeasure=rmeasure, dmeasure=dmeasure, skeleton = skeleton.kd, skeleton.type = "map", initializer=initializer, rprior=rprior.kd, dprior=dprior.kd, parameter.transform = par.trans, parameter.inv.transform = par.inv.trans)
# 
# np = 1000
# par.init = rprior.kd(1)
# #par.init = matrix(NA, nr=3, nc=np)
# #for(i in 1:np) par.init[,i] = log(rprior(i))
# sir.pf.kd = bsmc(sir.pomp.kd, params = par.init, Np = np, est = c('beta','gamma','nu'), smooth = 0.1, verbose=TRUE, transform = TRUE)
# 
# # Calculate 95% credible intervals for theta
# (ltheta = apply(sir.pf.kd$post, 1, function(x) quantile(x, .025)))
# (utheta = apply(sir.pf.kd$post, 1, function(x) quantile(x, .975)))
# mysims[[1]]$true.params$theta

# Run pmcmc
rw.sd = c(0.005, 0.001, 0.01); names(rw.sd) = c('beta','gamma','nu')
sir_pmcmc <- function(n.chain, niter, np)
{
  pars.init = rprior(1); names(pars.init) = c('beta','gamma','nu')
  out = pmcmc(sir.pomp, Nmcmc = niter, start = pars.init, pars = c('beta','gamma','nu'), rw.sd = rw.sd, Np = np, max.fail = Inf, verbose=TRUE)
  file = paste(dpath,"sir_pmcmc_test-",paste(n.chain,niter,np,sep="-"),".rdata",sep="")
  save(out, file = file)
}
require(plyr)
mydata = data.frame(n.chain=1:3, niter=1100, np=100)
m_ply(mydata, sir_pmcmc)
