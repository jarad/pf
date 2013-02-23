# Initialize .Random.seed
set.seed(sample(1:1000,1))

# Set graphics and data paths, param label
gpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF-D2/"
dpath = "C:/Users/Danny/My Documents/UCSB - Research/pf/data/D2/"
param = "-3P-normpriors"

# Set known parameter values
P = 5000
d = 5
b = c(.25, .27, .23, .29)
varsigma = c(1.07, 1.05, 1.01, .98)
sigma = c(.0012, .0008, .0010, .0011)
dpower = 2

# Set unknown parameter values
theta = c(0.2399, 0.1066, 1.2042)

# Functions to simulate epidemic
source("sir_functions.R")
revo_sim = function(x){ revo(x, P, d, theta)}
robs_sim = function(x){ robs(x, b, varsigma, sigma, dpower)}
rinit_sim = function(){ rinit(10/P)}
source("ss.sim.R")

# Create function to sample from prior distribution of theta
theta.mean = c(-1.3296, -2.1764, 0.1055)
theta.var = c(0.1055, 0.0140, 0.0064)
rtheta = function(){ rnorm(3,theta.mean,theta.var)}

# Create function to map theta from real line to original scale
ftheta = function(theta) exp(theta)

# Functions for bootstrap filter
dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
revo_bf = function(x){ revo(x, P, d, ftheta(x[3:5]))}
rprior_bf = function()
{ 
  myprior = rprior(sim$y[,1], rtheta, b, varsigma, sigma, dpower)
  return(c(myprior$x,myprior$theta))
}
source("bf.R")

# Functions for auxiliary particle filter
dllik_apf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
pstate_apf = function(x) { revo(x, P, d, ftheta(x[3:5]), FALSE)}
revo_apf = function(x){ revo(x, P, d, ftheta(x[3:5]))}
rprior_apf = function()
{ 
  myprior = rprior(sim$y[,1], rtheta, b, varsigma, sigma, dpower)
  return(c(myprior$x,myprior$theta))
}
source("apf.R")

# Functions for kernel density particle filter
dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma, dpower)}
pstate_kd = function(x, theta) { revo(x, P, d, ftheta(theta), FALSE)}
revo_kd = function(x, theta){ revo(x, P, d, ftheta(theta))}
rprior_kd = function(){ rprior(sim$y[,1], rtheta, b, varsigma, sigma, dpower)}
source("kd_pf.R")

# Simulate N epidemics of length nt days and run particle filters with n particles
N = 100
n = 100
nt = 125
np = length(theta)
ns = length(rinit_sim())
current.seed = .Random.seed
no = length(robs_sim(rinit_sim()))
.Random.seed = current.seed
bf.med = apf.med = kd.med = array(NA,dim=c(N,nt+1,ns+np))
sims_y = array(NA,dim=c(N,no,nt))
require(Hmisc)
for(j in 1:N)
{
  # Simulate epidemic
  sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
  sims_y[j,,] = sim$y
 
  # Run bootstrap filter
  out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n, method="stratified", nonuniformity="ess", threshold=0.8, log=F)

  # Run auxilliary particle filter
  out2 = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, method="stratified", nonuniformity="ess", threshold=0.8, log=F)

  # Run kernel density particle filter
  out3 = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, method="stratified", nonuniformity="ess", threshold=0.8, log=F)

  # Calculate median of states and unknown parameters
  for(k in 1:(ns+np))
  {  
    for(i in 1:(nt+1))
    { 
      if(k <= 2)
      {
        bf.med[j,i,k] = wtd.quantile(out$state[k,,i], out$weight[,i], normwt=T, probs=.5)
        apf.med[j,i,k] = wtd.quantile(out2$state[k,,i], out2$weight[,i], normwt=T, probs=.5)
        kd.med[j,i,k] = wtd.quantile(out3$state[k,,i], out3$weight[,i], normwt=T, probs=.5)
      } else {
        bf.med[j,i,k] = wtd.quantile(ftheta(out$state[k,,i]), out$weight[,i], normwt=T, probs=.5)
        apf.med[j,i,k] = wtd.quantile(ftheta(out2$state[k,,i]), out2$weight[,i], normwt=T, probs=.5)
	kd.med[j,i,k] = wtd.quantile(ftheta(out3$theta[k-ns,,i]), out3$weight[,i], normwt=T, probs=.5)
      }
    }
  }
  cat("\n",j,"\n")
}

# Calculate median of the medians and lower/upper quantiles
bf.m = apf.m = kd.m = bf.l = apf.l = kd.l = bf.u = apf.u = kd.u = matrix(NA,nr=nt+1,nc=ns+np)
for(k in 1:(ns+np))
{
  bf.m[,k] = apply(bf.med[,,k],2,median)
  bf.l[,k] = apply(bf.med[,,k],2,function(x) quantile(x,probs=.025))
  bf.u[,k] = apply(bf.med[,,k],2,function(x) quantile(x,probs=.975))
  apf.m[,k] = apply(apf.med[,,k],2,median)
  apf.l[,k] = apply(apf.med[,,k],2,function(x) quantile(x,probs=.025))
  apf.u[,k] = apply(apf.med[,,k],2,function(x) quantile(x,probs=.975))
  kd.m[,k] = apply(kd.med[,,k],2,median)
  kd.l[,k] = apply(kd.med[,,k],2,function(x) quantile(x,probs=.025))
  kd.u[,k] = apply(kd.med[,,k],2,function(x) quantile(x,probs=.975))
}

# Save data
save.image(paste(dpath,"sir.pf.med",param,"-",n,"-",N,".rdata",sep=""))
#load(paste(dpath,"sir.pf.med",param,"-",n,"-",N,".rdata",sep=""))

# Plot medians and lower/upper quantiles of parameters over time
expr = expression(beta,gamma,nu)
xlabs = c("Time (days)","","")
ylabs = c(paste("J = ",n,sep=""),"","")
msize = 2.4
labsize = 2.2
legendsize = 1.5
pdf(paste(gpath,"PF-params-med",param,"-",n,"-",N,".pdf",sep=""),width=15,height=5)
par(mfrow=c(1,3),mar=c(5,5,3,1)+.1)
for(j in 1:np)
{
  ymax = max(bf.u[,j+ns],apf.u[,j+ns],kd.u[,j+ns],sim$x[j+ns,])
  ymin = min(bf.l[,j+ns],apf.l[,j+ns],kd.l[,j+ns],sim$x[j+ns,])
  if(n == 100)
  {
    plot(0:nt,sim$x[j+ns,],type="l",ylim=c(ymin,ymax),col=2,,xlab="",ylab=ylabs[j],main=expr[j],cex.lab=labsize,cex.main=msize)
    if(j == 1) legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2),cex=legendsize)
  } else if(n == 10000) {
    plot(0:nt,sim$x[j+ns,],type="l",ylim=c(ymin,ymax),col=2,,xlab=xlabs[j],ylab=ylabs[j],cex.lab=labsize)
  } else {
    plot(0:nt,sim$x[j+ns,],type="l",ylim=c(ymin,ymax),col=2,ylab=ylabs[j],xlab="",cex.lab=labsize)
  }
  lines(0:nt,bf.m[,j+ns],col=2)
  lines(0:nt,bf.l[,j+ns],col=2,lty=2)
  lines(0:nt,bf.u[,j+ns],col=2,lty=2)
  lines(0:nt,apf.m[,j+ns],col=4)
  lines(0:nt,apf.l[,j+ns],col=4,lty=2)
  lines(0:nt,apf.u[,j+ns],col=4,lty=2)
  lines(0:nt,kd.m[,j+ns],col=3)
  lines(0:nt,kd.l[,j+ns],col=3,lty=2)
  lines(0:nt,kd.u[,j+ns],col=3,lty=2)
}
dev.off()

# Plot medians and upper/lower quantiles of %infected/susceptible over time
mlabs = expression(i,s)
xlabs = c("Time (days)","")
ylabs = c(paste("J = ",n,sep=""),"")
pdf(paste(gpath,"PF-states-med",param,"-",n,"-",N,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2),mar=c(5,5,3,1)+0.1)
for(k in 1:2)
{
  if(k == 1) {
    ymax = max(bf.u[,k],apf.u[,k],kd.u[,k],sim$x[k,])
    ymin = 0
  } else {
    ymax = 1
    ymin = min(bf.l[,k],apf.l[,k],kd.l[,k],sim$x[k,])
  }
  plot(0:nt,sim$x[k,],type="l",ylim=c(ymin,ymax),xlab=xlabs[k],ylab=ylabs[k],main=mlabs[k],cex.lab=1.5,cex.main=1.75)
  legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2))
  lines(0:nt,bf.m[,k],col=2)
  lines(0:nt,bf.l[,k],col=2,lty=2)
  lines(0:nt,bf.u[,k],col=2,lty=2)
  lines(0:nt,apf.m[,k],col=4)
  lines(0:nt,apf.l[,k],col=4,lty=2)
  lines(0:nt,apf.u[,k],col=4,lty=2)
  lines(0:nt,kd.m[,k],col=3)
  lines(0:nt,kd.l[,k],col=3,lty=2)
  lines(0:nt,kd.u[,k],col=3,lty=2)
}
dev.off()