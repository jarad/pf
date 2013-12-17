source("sir_functions.r")
source("sir_mcmc_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# Load simulated data and redefine data / known parameter values
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysim$sim[[1]]$y
P = mysim$true.params$P
b = mysim$true.params$b
varsigma = mysim$true.params$varsigma
sigma = mysim$true.params$sigma
eta = mysim$true.params$eta

# Set up initial values of MCMC
initial.chains = list()
n.chains = 3
nt = dim(y)[2]
for(j in 1:n.chains)
{
  set.seed(60 + j)
  log.params = find.mu.sigma(c(.14, .09, .95), c(.5, .143, 1.3))
  tmp = rprior(function() exp(rnorm(3, log.params[[1]], log.params[[2]])))
  theta = tmp$theta
  x = mysim$sim[[1]]$x
  initial.chains[[j]] = list(x = x, theta = theta)
}

# Test mcmc
psi = list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P)
tuning.1 = list(tuning.x = matrix(0.001, nr=2, nc = nt+1), tuning.theta = c(0.01, 0.001, 0.01))
tuning.2 = list(tuning.x = matrix(c(0.002, 0.001), nr=2, nc = nt+1), tuning.theta = c(0.01, 0.001, 0.01))
mcmc.details = list(n.thin=1000, n.sims=11000000, n.burn=1000000, tune = TRUE)
my_sir_mcmc <- function(n.chain, tune.type, x, beta, gamma, nu, progress, print.iter)
{
  if(tune.type == 1) tuning = tuning.1 else tuning = tuning.2
  steps = c('x','beta','gamma','nu')
  params.est <- which(as.logical(c(x,beta,gamma,nu)))
  steps = steps[params.est]
  if(!('x' %in% steps)) initial.chains[[n.chain]]$x = mysim$sim[[1]]$x
  if(!('beta' %in% steps)) initial.chains[[n.chain]]$theta[1] = mysim$true.params$theta[1]
  if(!('gamma' %in% steps)) initial.chains[[n.chain]]$theta[2] = mysim$true.params$theta[2]
  if(!('nu' %in% steps)) initial.chains[[n.chain]]$theta[3] = mysim$true.params$theta[3]
  cat(n.chain,tune.type,x,beta,gamma,nu,"\n",sep=" ")
  out = sir_mcmc(y, psi, initial.chains[[n.chain]], tuning, mcmc.details, steps, progress, print.iter)
  return(list(out=out,samp.details=list(n.chain=n.chain,tune.type=tune.type,x=x,beta=beta,gamma=gamma,nu=nu)))
}
require(plyr)
require(doMC)
registerDoMC()
mydata = data.frame(n.chain=1:3, tune.type=1, x=1, beta=1, gamma=1, nu=1, progress=FALSE, print.iter=TRUE)
mcmc.chains = mlply(.data = mydata, .fun = my_sir_mcmc, .parallel = TRUE)
for(i in 1:dim(mydata)[1])
{
  label = mcmc.chains[[i]]$samp.details
  file = paste(dpath,"sir_mcmc_test-",paste(label$n.chain,label$tune.type,label$x,label$beta,label$gamma,label$nu,sep="-"),".rdata",sep="")
  out = mcmc.chains[[i]]$out
  save(out, file = file)
}