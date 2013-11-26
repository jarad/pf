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

# Test mcmc
psi = list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P)
mcmc.details = list(n.thin=1, n.sims=100000, n.burn=10000, tune = TRUE)
my_sir_mcmc <- function(n.chain) sir_mcmc(y, psi, mcmc.details=mcmc.details)
require(plyr)
require(doMC)
registerDoMC()
out = mlply(.data = data.frame(n.chain=1:3), .fun = my_sir_mcmc, .parallel = TRUE)
save(out, file = paste(dpath,"sir_mcmc_test.rdata",sep=""))