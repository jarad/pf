source("sir_functions.r")
source("sir_mcmc_functions.r")
source("pf_functions.r")

# Set data path
dpath = "../data/"

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
 
# Set seed
set.seed(61)

# Run KDPF
mydllik = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
myrevo = function(x, theta) revo(x, ftheta(theta), P)
pstate = function(x, theta) revo(x, ftheta(theta), P, FALSE)
myrprior = function() rprior(rtheta)
source("kd_pf.r")
time = system.time(out <- kd_pf(y, mydllik, pstate, myrevo, myrprior, n, delta, progress, method=resamp, nonuniformity = "ess", threshold = 0.8, log=F))

# Calculate 95% credible intervals
probs = c(0.025, 0.975)
tt = dim(out$state)[3]
states = array(NA,dim=c(3,n,tt))
states[1:2,,] = out$state[1:2,,]
for(i in 1:tt) states[3,,i] = 1 - apply(out$state[1:2,,i],2,sum)
state.quant = pf.quantile(states, out$weight, function(x,param=1) x, probs)
theta.quant = pf.quantile(out$theta, out$weight, ftheta, probs)

# Print output
file = paste(dpath,"PF-1-",n,"-KD-stratified-orig-log-0.99-61.rdata",sep="")
cat(file,time[1:3],"\n")
tab.quant <- aperm(theta.quant[c(31,61,91,126),,], c(3,2,1))
dimnames(tab.quant) = list(c("2.5%","97.5%"),c("beta","gamma","nu"),c(30,60,90,125))
print(tab.quant)
