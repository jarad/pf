source("sir_functions.r")
source("sir_mcmc_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# Load simulated data and redefine data / known parameter values
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysim$sim$y
P = mysim$true.params$P
b = mysim$true.params$b
varsigma = mysim$true.params$varsigma
sigma = mysim$true.params$sigma
eta = mysim$true.params$eta
nu = mysim$true.params$theta[3]

# Test mcmc
psi = list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P)
mcmc.details = list(n.thin=1, n.sims=100000, n.burn=0)
n.chains = 3
my_sir_mcmc <- function(n) sir_mcmc(y, psi, mcmc.details, verbose=T)
mydata = data.frame(n=1:n.chains)
require(plyr)
require(doMC)
registerDoMC()
out = mlply(.data = mydata, .fun = my_sir_mcmc, .parallel = TRUE)
save(out, file = paste(dpath,"sir_mcmc_test.rdata",sep=""))