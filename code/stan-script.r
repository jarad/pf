library(rstan)
set_cppo("fast")  # for best running speed

# Load simulated data
load("../data/sim-orig.rdata")

# Set data values
N = 60
y = mysim$sim$y[,1:N]
obs1 = which(!is.na(y[1,])); N1 = length(obs1)
obs2 = which(!is.na(y[2,])); N2 = length(obs2)
obs3 = which(!is.na(y[3,])); N3 = length(obs3)
obs4 = which(!is.na(y[4,])); N4 = length(obs4)
y1 = y[1,obs1]
y2 = y[2,obs2]
y3 = y[3,obs3]
y4 = y[4,obs4]
b = mysim$true.params$b
varsigma = mysim$true.params$varsigma
sigma = mysim$true.params$sigma
eta = mysim$true.params$eta
P = mysim$true.params$P
nu = mysim$true.params$theta[3]

# Initial values for each chain
n.chains = 3
d = list(y1=y1,y2=y2,y3=y3,y4=y4,obs1=obs1,obs2=obs2,obs3=obs3,obs4=obs4,N=N,N1=N1,N2=N2,N3=N3,N4=N4,b=b,varsigma=varsigma,sigma=sigma,eta=eta,P=P,nu=nu)
beta.init = rlnorm(N, -1.3296, 0.1183)
gamma.init = rlnorm(N, -2.1764, 0.1055)
i0 = rnorm(n.chains, 0.002, 0.0005)
while(!all(i0 > 0 & i0 < 1)) i0 = rnorm(n.chains, 0.002, 0.0005)
x = array(NA, c(2,N,n.chains))
for(j in 1:n.chains)
{
  x[1,1,j] = 1 - i0[j] - beta.init[1]*i0[j]*(1-i0[j])^nu + rnorm(1, 0, sqrt(beta.init[1] / P^2))
  x[2,1,j] = i0[j] + beta.init[1]*i0[j]*(1-i0[j])^nu - gamma.init[1]*i0[j] + rnorm(1, -1*(x[1,1,j] - ((1-i0[j]) - beta.init[1]*i0[j]*(1-i0[j])^nu)), sqrt(gamma.init[1] / P^2))
  for(i in 1:(N-1)){
    x[1,i+1,j] = x[1,i,j] - beta.init[i]*x[2,i,j]*x[1,i,j]^nu + rnorm(1, 0, sqrt(beta.init[i] / P^2))
    x[2,i+1,j] = x[2,i,j] + beta.init[i]*x[2,i,j]*x[1,i,j]^nu - gamma.init[i]*x[2,i,j] + rnorm(1, -1*(x[1,i+1,j] - (x[1,i,j] - beta.init[i]*x[2,i,j]*x[1,i,j]^nu)), sqrt(gamma.init[i] / P^2))
  }
}
init = list()
for(i in 1:n.chains) init[[i]] = list(i0 = i0[i], beta = beta.init[i], gamma = gamma.init[i], x = x[,,i])

# Fit stan model
fit <- stan(file = 'sir-model.stan', init = init, data = d, iter = 1000, chains = 3)