source("sir_functions.r")
source("sir_mcmc_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# Load simulated data and redefine data / known parameter values
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysims[[1]]$sim$y
P = mysims[[1]]$true.params$P
b = mysims[[1]]$true.params$b
varsigma = mysims[[1]]$true.params$varsigma
sigma = mysims[[1]]$true.params$sigma
eta = mysims[[1]]$true.params$eta

# Test mcmc
psi = list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P)
mcmc.details = list(n.thin=1000, n.sims=10100000, n.burn=100000, tune = TRUE)
my_sir_mcmc <- function(n.chain, x, beta, gamma, nu, ymax, progress, print.iter)
{
  if(missing(ymax)) ymax = dim(y)[2]
  t = 1:ymax
  nt = length(t)
  
  # Set initial values
  set.seed(60 + n.chain)  
  rtheta <- function()
  {
    theta <- rep(NA, 3)
    log.params <- find.mu.sigma(c(1.5, .09, .95), c(3, .143, 1.3))
    theta[2:3] <- exp(rnorm(2, log.params[[1]][2:3], log.params[[2]][2:3]))
    theta[1] <- theta[2]*exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
    return(theta)
  }  
  tmp = rprior(rtheta)
  theta = tmp$theta
  x.init = mysims[[1]]$sim$x[,c(t,ymax+1)]
  initial.chains = list(x = x.init, theta = theta)
  
  tuning = list(tuning.x = matrix(0.001, nr=2, nc = nt+1), tuning.theta = c(0.01, 0.001, 0.01))
  steps = c('x','beta','gamma','nu')
  params.est <- which(as.logical(c(x,beta,gamma,nu)))
  steps = steps[params.est]
  if(!('x' %in% steps)) initial.chains$x = mysims[[1]]$sim$x
  if(!('beta' %in% steps)) initial.chains$theta[1] = mysims[[1]]$true.params$theta[1]
  if(!('gamma' %in% steps)) initial.chains$theta[2] = mysims[[1]]$true.params$theta[2]
  if(!('nu' %in% steps)) initial.chains$theta[3] = mysims[[1]]$true.params$theta[3]

  cat(n.chain,x,beta,gamma,nu,ymax,"\n",sep=" ")
  out = sir_mcmc(y[,t], psi, initial.chains, tuning, mcmc.details, steps, progress, print.iter)
  file = paste(dpath,"sir_mcmc_test-",paste(n.chain,x,beta,gamma,nu,ymax,sep="-"),".rdata",sep="")
  save(out, file = file)
}

require(plyr)
require(doMC)
registerDoMC()

mydata = data.frame(n.chain=1:3, x=1, beta=1, gamma=1, nu=1, ymax = 125, progress=FALSE, print.iter=TRUE)
m_ply(.data = mydata, .fun = my_sir_mcmc, .parallel = TRUE)