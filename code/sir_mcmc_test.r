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
#  x = matrix(NA, 2, nt + 1)
#  x[,1] = tmp$x
#  beta.x <- exp(rnorm(nt, log.params[[1]][1], log.params[[2]][1]))
#  gamma.x <- exp(rnorm(nt, log.params[[1]][2], log.params[[2]][2]))
#  nu.x <- exp(rnorm(nt, log.params[[1]][3], log.params[[2]][3]))
#  for(i in 1:nt) x[,i+1] = revo(x[,i], c(beta.x[i], gamma.x[i], nu.x[i]), P)
  x = mysim$sim[[1]]$x
  initial.chains[[j]] = list(x = x, theta = theta)
}

# Test mcmc
psi = list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P)
tuning.1 = list(tuning.x = matrix(0.001, nr=2, nc = nt+1), tuning.theta = c(0.01, 0.001, 0.01))
tuning.2 = list(tuning.x = matrix(c(0.002, 0.001), nr=2, nc = nt+1), tuning.theta = c(0.01, 0.001, 0.01))
mcmc.details = list(n.thin=50, n.sims=600000, n.burn=100000, tune = TRUE)
my_sir_mcmc <- function(n.chain, tune.type, x, beta, gamma, nu, progress)
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
  out = sir_mcmc(y, psi, initial.chains[[n.chain]], tuning, mcmc.details, steps, progress)
  save(out, file = paste(dpath,"sir_mcmc_test-",paste(n.chain,tune.type,x,beta,gamma,nu,sep="-"),".rdata",sep=""))
}
require(plyr)
require(doMC)
registerDoMC()
mydata = matrix(nr=0,nc=0)
#data = data.frame(x=c(0,0,0,1,1),beta=c(1,0,0,0,1),gamma=c(0,1,0,0,1),nu=c(0,0,1,0,1),progress=FALSE)
data = data.frame(x=c(1,1),beta=c(0,1),gamma=c(0,1),nu=c(0,1),progress=FALSE)
for(k in 1:dim(data)[1])
{
  for(i in 1:2)
  {
    for(j in 1:n.chains)
    {
      mydata = rbind(mydata, data.frame(n.chain=j,tune.type=i,data[k,]))
    }
  }
}
rownames(mydata) = 1:dim(mydata)[1]
m_ply(.data = mydata, .fun = my_sir_mcmc, .parallel = TRUE)