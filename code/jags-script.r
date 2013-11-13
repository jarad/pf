require(rjags)

# Load simulated data
load("../storage/sheinson_research/sim-orig.rdata")

# Set data values
y = mysim$sim$y[,1:60]
N = dim(y)[2]
L = length(mysim$true.params$b)
b = mysim$true.params$b
varsigma = mysim$true.params$varsigma
sigma = mysim$true.params$sigma
eta = mysim$true.params$eta
P = mysim$true.params$P

# Test JAGS model
d = list(y=y,N=N,L=L,b=b,varsigma=varsigma,sigma=sigma,eta=eta,P=P)
mod = jags.model("jags-model.txt", data=d, n.chains=1, n.adapt=1e3)
samps = coda.samples(mod, c("beta","gamma","nu"), n.iter=1e5)
save(samps, file="/storage/sheinson_research/sir-jags.rdata")

# Diagnostics
print(gelman.diag(samps))
print(summary(sapply(samps,function(x) x[90000:100000,])))

# Trace plots
pdf(file="../graphs/sir-jags-traceplots.pdf")
par(mfrow=c(3,1))
plot(1:100000,samps[[1]][,1],type="l")
lines(1:100000,samps[[2]][,1],col=2)
lines(1:100000,samps[[3]][,1],col=4)
plot(1:100000,samps[[1]][,2],type="l")
lines(1:100000,samps[[2]][,2],col=2)
lines(1:100000,samps[[3]][,2],col=4)
plot(1:100000,samps[[1]][,3],type="l")
lines(1:100000,samps[[2]][,3],col=2)
lines(1:100000,samps[[3]][,3],col=4)
dev.off()

