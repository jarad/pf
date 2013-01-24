# Set known parameter values
P = 5000
d = 5
b = c(.25,.27,.23,.29)
varsigma = c(1.07,1.05,1.01,.98)
sigma = c(.0012,.0008,.001,.0011)

# Set unknown parameter values
beta = 0.2399
gamma = 0.1066
nu = 1.2042

# Functions for simulating data 
rstate = function(x)
{
  tau = 1/d
  x1var = (x[3] + x[4])*tau^2 / P^2
  x2var = x[3]*tau^2 / P^2
  for(i in 1:d)
  { 
    is = -1
    ss = -1
    while(!(is >= 0 & is <= 1 & ss >= 0 & ss <= 1))
    {
      is = x[1] + tau*x[1]*(x[3]*x[2]^x[5] - x[4]) + rnorm(1,0,sqrt(x1var))
      ss = x[2] - tau*x[3]*x[1]*x[2]^x[5] + rnorm(1,0,sqrt(x2var))
    }
    x[1] = is
    x[2] = ss
  }
  return(x)
}
robs = function(x)
{
  j = length(b)
  if(!(j == length(varsigma) & j == length(sigma))) stop("b, varsigma, and sigma must all have same length")
  J = sample(1:j,1)
  y = rep(NA,j)
  y[J] = rnorm(1,b[J]*x[1]^varsigma[J],sigma[J])
  return(y)
}
rinit = function() return(c(10/P,1-10/P,beta,gamma,nu))

# Generate data
source("ss.sim.R")
nt = 154
sim = ss.sim(nt,rstate,robs,rinit)

# Functions to reparameterize theta to and from the real line
theta2u = function(theta,a,b)
{
  etheta = exp(theta)
  return((b*etheta + a) / (1 + etheta))
}
u2theta = function(u,a,b)
{
  U = (u - a) / (b - a)
  return(log(U / (1 - U)))
}

# Set bounds on unknown parameters
betal = .14; betau = .5
gammal = .09; gammau = .143
nul = .95; nuu = 1.3

# Functions for particle filters
revo = function(x)
{
  tau = 1/d
  x1var = (x[3] + x[4])*tau^2 / P^2
  x2var = x[3]*tau^2 / P^2
  for(i in 1:d)
  { 
    is = -1
    ss = -1
    while(!(is >= 0 & is <= 1 & ss >= 0 & ss <= 1))
    {
      is = x[1] + tau*x[1]*(x[3]*x[2]^x[5] - x[4]) + rnorm(1,0,sqrt(x1var))
      ss = x[2] - tau*x[3]*x[1]*x[2]^x[5] + rnorm(1,0,sqrt(x2var))
    }
    x[1] = is
    x[2] = ss
  }
  return(x)
}
revo_kd = function(x,theta=NULL)
{
  theta[1] = theta2u(theta[1],betal,betau)
  theta[2] = theta2u(theta[2],gammal,gammau)
  theta[3] = theta2u(theta[3],nul,nuu)

  tau = 1/d
  x1var = (theta[1] + theta[2])*tau^2 / P^2
  x2var = theta[1]*tau^2 / P^2
  for(i in 1:d)
  { 
    is = -1
    ss = -1
    while(!(is >= 0 & is <= 1 & ss >= 0 & ss <= 1))
    {
      is = x[1] + tau*x[1]*(theta[1]*x[2]^theta[3] - theta[2]) + rnorm(1,0,sqrt(x1var))
      ss = x[2] - tau*theta[1]*x[1]*x[2]^theta[3] + rnorm(1,0,sqrt(x2var))
    }
    x[1] = is
    x[2] = ss
  }
  return(x)
}
dllik = function(y,x,theta=NULL)
{
  J = which(!is.na(y))
  if (!(length(J) == 1)) stop("y must be a vector with all but 1 element NA")
  h = b[J]*x[1]^varsigma[J]
  return(dnorm(y[J],h,sigma[J],log=T))
}
rprior = function(stateonly=TRUE)
{
  require(msm)
  J = which(!is.na(sim$y[,1]))
  if (!(length(J) == 1)) stop("y must be a vector with all but 1 element NA")
  
  i0 = rtnorm(1,sim$y[J,1]/b[J],sigma[J],0,1)
  s0 = rtnorm(1,1-sim$y[J,1]/b[J],sigma[J],0,1)
  beta0 = runif(1,.14,.500)
  gamma0 = runif(1,.09,.143)
  nu0 = runif(1,.95,1.3)
  if(stateonly){
    return(c(i0,s0,beta0,gamma0,nu0))
  }else{
    return(list(x=c(i0,s0),theta=u2theta(c(beta0,gamma0,nu0),c(betal,gammal,nul),c(betau,gammau,nuu))))
  }
}
pstate = function(x)
{
  tau = 1/d
  for(i in 1:d)
  { 
    x[1] = x[1] + tau*x[1]*(x[3]*x[2]^x[5] - x[4])
    x[2] = x[2] - tau*x[3]*x[1]*x[2]^x[5]
  }
  return(x)
}
pstate_kd = function(x,theta=NULL)
{
  theta[1] = theta2u(theta[1],betal,betau)
  theta[2] = theta2u(theta[2],gammal,gammau)
  theta[3] = theta2u(theta[3],nul,nuu)

  tau = 1/d
  for(i in 1:d)
  { 
    x[1] = x[1] + tau*x[1]*(theta[1]*x[2]^theta[3] - theta[2])
    x[2] = x[2] - tau*theta[1]*x[1]*x[2]^theta[3]
  }
  return(x)
}

# Run bootstrap filter
source("bf.R")
n=100
out = bf(sim$y, dllik, revo, rprior, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run auxilliary particle filter
source("apf.R")
out2 = apf(sim$y, dllik, pstate, revo, rprior, n,
	   method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run kernel density particle filter
source("kd_pf.R")
rprior_kd = function() rprior(FALSE)
out3 = kd_pf(sim$y, dllik, pstate_kd, revo_kd, rprior_kd, n,
 	   method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Comparison
# Plot % infected
tt = nt + 1
plot(sim$x[1,],type="l",ylim=c(0,.25),xlab="Time (days)",ylab="% Population",
	main="95% Credible Intervals of %Pop Infected")
bf.i = apply(out$state[1,,]*out$weight,2,sum)
apf.i = apply(out2$state[1,,]*out2$weight,2,sum)
kd.i = apply(out3$state[1,,]*out3$weight,2,sum)
bf.li = bf.ui = apf.li = apf.ui = kd.li = kd.ui = numeric(tt)
for(i in 1:tt){
	bf.li[i] = wtd.quantile(out$state[1,,i], out$weight[,i], normwt=T, probs=.025)
	bf.ui[i] = wtd.quantile(out$state[1,,i], out$weight[,i], normwt=T, probs=.975)
	apf.li[i] = wtd.quantile(out2$state[1,,i], out2$weight[,i], normwt=T, probs=.025)
	apf.ui[i] = wtd.quantile(out2$state[1,,i], out2$weight[,i], normwt=T, probs=.975)
	kd.li[i] = wtd.quantile(out3$state[1,,i], out3$weight[,i], normwt=T, probs=.025)
	kd.ui[i] = wtd.quantile(out3$state[1,,i], out3$weight[,i], normwt=T, probs=.975)
}
lines(bf.i,col=2)
lines(bf.li,col=2,lty=2)
lines(bf.ui,col=2,lty=2)
lines(apf.i,col=4)
lines(apf.li,col=4,lty=2)
lines(apf.ui,col=4,lty=2)
lines(kd.i,col=3)
lines(kd.li,col=3,lty=2)
lines(kd.ui,col=3,lty=2)
legend("topright",c("Truth","BF","AFP","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2))

# Histograms of unknown parameters
# bootstrap filter
require(plotrix)
cutoff = 25
windows(width=5,height=8)
par(mfrow=c(3,1))
weighted.hist(out$state[3,,cutoff],out$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate")
mtext(expression(paste(beta," = ",0.2399,sep="")),side=3,cex=.8)
weighted.hist(out$state[4,,cutoff],out$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time")
mtext(expression(paste(gamma," = ",0.1066,sep="")),side=3)
weighted.hist(out$state[5,,cutoff],out$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity")
mtext(expression(paste(nu," = ",1.2042,sep="")),side=3)

# auxillary particle filter
windows(width=5,height=8)
par(mfrow=c(3,1))
weighted.hist(out2$state[3,,cutoff],out2$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate")
mtext(expression(paste(beta," = ",0.2399,sep="")),side=3,cex=.8)
weighted.hist(out2$state[4,,cutoff],out2$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time")
mtext(expression(paste(gamma," = ",0.1066,sep="")),side=3)
weighted.hist(out2$state[5,,cutoff],out2$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity")
mtext(expression(paste(nu," = ",1.2042,sep="")),side=3)

# kernel density particle filter
windows(width=5,height=8)
par(mfrow=c(3,1))
weighted.hist(theta2u(out3$theta[1,,cutoff],betal,betau),out3$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate")
mtext(expression(paste(beta," = ",0.2399,sep="")),side=3,cex=.8)
weighted.hist(theta2u(out3$theta[2,,cutoff],gammal,gammau),out3$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time")
mtext(expression(paste(gamma," = ",0.1066,sep="")),side=3)
weighted.hist(theta2u(out3$theta[3,,cutoff],nul,nuu),out3$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity")
mtext(expression(paste(nu," = ",1.2042,sep="")),side=3)
