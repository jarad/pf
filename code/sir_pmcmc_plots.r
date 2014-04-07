# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data and redefine data / known parameter values
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysims[[1]]$sim$y
P = mysims[[1]]$true.params$P
b = mysims[[1]]$true.params$b
varsigma = mysims[[1]]$true.params$varsigma
sigma = mysims[[1]]$true.params$sigma
eta = mysims[[1]]$true.params$eta

# Load PMCMC objects
n.chains = 3
load(paste(dpath,"sir_pmcmc_test-1-1100-100.rdata",sep=""))
niter = dim(conv.rec(out))[1] - 1
nthin = 1
nburn = 100
out.beta = matrix(NA, nr = (niter - nburn) %/% nthin, nc=n.chains)
out.gamma = matrix(NA, nr = (niter - nburn) %/% nthin, nc=n.chains)
out.nu = matrix(NA, nr = (niter - nburn) %/% nthin, nc=n.chains)
out.theta = list(); length(out.theta) = n.chains
iter = seq(nburn+nthin+1,niter+1,nthin)
for(i in 1:n.chains)
{
  load(paste(dpath,"sir_pmcmc_test-",i,"-1100-100.rdata",sep=""))
  out.theta[[i]] = conv.rec(out)[,4:6]
  theta.temp = conv.rec(out)[iter,4:6]
  out.beta[,i] = theta.temp[,1]
  out.gamma[,i] = theta.temp[,2]
  out.nu[,i] = theta.temp[,3]
}

# Compute 95% credible intervals
cred.int = matrix(NA, nr=2, nc = 3)
colnames(cred.int) = c("beta","gamma","nu")
for(i in 1:3) cred.int[,i] = quantile(sapply(out.theta, function(x) x[,i]),c(0.025,0.975))
print(cred.int)

# Construct traceplots
windows()
par(mfrow=c(n.chains,1),mar=c(5,7,4,2)+0.1)
plot(iter, out.beta[,1], type="l", ylim = c(min(out.beta),max(out.beta)), xlab="", ylab = expression(beta), cex.lab = 1.75)
lines(iter, out.beta[,2], col = 2)
lines(iter, out.beta[,3], col = 4)
plot(iter, out.gamma[,1], type="l", ylim = c(min(out.gamma),max(out.gamma)), xlab="", ylab = expression(gamma), cex.lab = 1.75)
lines(iter, out.gamma[,2], col = 2)
lines(iter, out.gamma[,3], col = 4)
plot(iter, out.nu[,1], type="l", ylim = c(min(out.nu),max(out.nu)), xlab="Iteration", ylab = expression(nu), cex.lab = 1.75)
lines(iter, out.nu[,2], col = 2)
lines(iter, out.nu[,3], col = 4)

