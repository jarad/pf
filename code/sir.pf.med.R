# Initialize .Random.seed
set.seed(sample(1:1000,1))

# Set graphics and data paths
gpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/"
dpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/"

# How many unknown parameters? Set p = 3 or p = 6
p = 3
if(p == 3) param = "" else param = "-6P"

# Set known parameter values
P = 5000
d = 5
b = c(.25, .27, .23, .29)
varsigma = c(1.07, 1.05, 1.01, .98)
sigma = c(.0012, .0008, .0010, .0011)

# Set unknown parameter values and prior bounds
theta = c(0.2399, 0.1066, 1.2042, b[1], varsigma[1], sigma[1])
thetal = c(0.1400, 0.0900, 0.9500, 0.00, 0.85, 0.0001)
thetau = c(0.5000, 0.1430, 1.3000, 0.50, 1.15, 0.0050)
theta = theta[1:p]
thetal = thetal[1:p]
thetau = thetau[1:p]

# Functions to simulate epidemic
source("sir_functions.R")
revo_sim = function(x){ revo(x, P, d)}
robs_sim = function(x){ robs(x, b[1:s], varsigma[1:s], sigma[1:s])}
rinit_sim = function(){ rinit(10/P, theta)}

if(p == 3)
{
  # 3 parameter PF functions
  # Functions for bootstrap filter
  dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma)}
  revo_bf = function(x){ revo(x, P, d)}
  rprior_bf = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma)}

  # Functions for auxiliary particle filter
  dllik_apf = function(y, x){ dllik(y, x, b, varsigma, sigma)}
  pstate_apf = function(x) { revo(x, P, d, random = FALSE)}
  revo_apf = function(x){ revo(x, P, d)}
  rprior_apf = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma)}

  # Functions for kernel density particle filter
  dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma)}
  pstate_kd = function(x, theta) { revo(x, P, d, theta, thetal, thetau, FALSE, FALSE)}
  revo_kd = function(x, theta){ revo(x, P, d, theta, thetal, thetau, FALSE)}
  rprior_kd = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma, FALSE)}
} else {
  # 6 parameter PF functions
  # Run bootstrap filter
  dllik_bf = function(y, x){ dllik(y, x, NULL, NULL, NULL, addparam=TRUE)}
  revo_bf = function(x){ revo(x, P, d)}
  rprior_bf = function(){ rprior(sim, thetal, thetau, addparam=TRUE)}

  # Run auxiliary particle filter
  dllik_apf = function(y, x){ dllik(y, x, NULL, NULL, NULL, addparam=TRUE)}
  pstate_apf = function(x) { revo(x, P, d, random = FALSE)}
  revo_apf = function(x){ revo(x, P, d)}
  rprior_apf = function(){ rprior(sim, thetal, thetau, addparam=TRUE)}

  # Run kernel density particle filter
  dllik_kd = function(y, x, theta){ dllik(y, x, NULL, NULL, NULL, theta, thetal, thetau, FALSE, TRUE)}
  pstate_kd = function(x, theta) { revo(x, P, d, theta, thetal, thetau, FALSE, FALSE)}
  revo_kd = function(x, theta){ revo(x, P, d, theta, thetal, thetau, FALSE)}
  rprior_kd = function(){ rprior(sim, thetal, thetau, stateonly=FALSE, addparam=TRUE)}
}

# Run particle filters N times with n particles
N = 100
n = 100
nt = 125
ns = length(rinit_sim())
current.seed = .Random.seed
no = length(robs_sim(rinit_sim()))
.Random.seed = current.seed
bf.med = apf.med = kd.med = array(NA,dim=c(N,nt+1,p+2))
sims = list(x=array(NA,dim=c(N,ns,nt+1)),y=array(NA,dim=c(N,no,nt)))
source("ss.sim.R")
source("bf.R")
source("apf.R")
source("kd_pf.R")
require(Hmisc)
for(j in 1:N)
{
  # Simulate epidemic
  sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
  sims$x[j,,] = sim$x
  sims$y[j,,] = sim$y
 
  # Run bootstrap filter
  out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

  # Run auxilliary particle filter
  out2 = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

  # Run kernel density particle filter
  out3 = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

  # Calculate median of states and unknown parameters
  for(k in 1:(p+2))
  {  
    for(i in 1:(nt+1)){
	  bf.med[j,i,k] = wtd.quantile(out$state[k,,i], out$weight[,i], normwt=T, probs=.5)
	  apf.med[j,i,k] = wtd.quantile(out2$state[k,,i], out2$weight[,i], normwt=T, probs=.5)
	  if(k <= 2)
        {
          kd.med[j,i,k] = wtd.quantile(out3$state[k,,i], out3$weight[,i], normwt=T, probs=.5)
        } else {
	    kd.med[j,i,k] = wtd.quantile(theta2u(out3$theta[k-2,,i],thetal[k-2],thetau[k-2]), out3$weight[,i], normwt=T, probs=.5)
	  }
    }
  }
  cat("\n",j,"\n")
}

# Calculate median of the medians and lower/upper quantiles
bf.m = apf.m = kd.m = bf.l = apf.l = kd.l = bf.u = apf.u = kd.u = matrix(NA,nr=nt+1,nc=p+2)
truth.m = truth.l = truth.u = matrix(NA,nr=nt+1,nc=p+2)
for(k in 1:(p+2))
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
  truth.m[,k] = apply(sims$x[,k,],2,median)
  truth.l[,k] = apply(sims$x[,k,],2,function(x) quantile(x,probs=.025))
  truth.u[,k] = apply(sims$x[,k,],2,function(x) quantile(x,probs=.975))
}

# Save data
save.image(paste(dpath,"sir.pf.med",param,"-",n,"-",N,".rdata",sep=""))
#load(paste(dpath,"sir.pf.med",param,"-",n,"-",N,".rdata",sep=""))

# Plot medians and lower/upper quantiles of parameters over time
expr = expression(beta,gamma,nu,b,varsigma,sigma)
xlabs = c("Time (days)",rep("",p-1))
labsize = 2.5
pdf(paste(gpath,"PF-params-med",param,"-",n,"-",N,".pdf",sep=""),width=15,height=5*((p == 6)+1))
par(mfrow=c(1+(p == 6),3),mar=c(5,6,4,0)+0.1)
for(k in 1:p)
{
  ymax = max(bf.u[,k+2],apf.u[,k+2],kd.u[,k+2],truth.u[,k+2])
  ymin = min(bf.l[,k+2],apf.l[,k+2],kd.l[,k+2],truth.l[,k+2])
  plot(0:nt,truth.m[,k+2],type="l",ylim=c(ymin,ymax),
  	xlab=xlabs[k],ylab=expr[k],cex.lab=labsize)
  lines(0:nt,truth.l[,k+2],lty=2)
  lines(0:nt,truth.u[,k+2],lty=2)
  lines(0:nt,bf.m[,k+2],col=2)
  lines(0:nt,bf.l[,k+2],col=2,lty=2)
  lines(0:nt,bf.u[,k+2],col=2,lty=2)
  lines(0:nt,apf.m[,k+2],col=4)
  lines(0:nt,apf.l[,k+2],col=4,lty=2)
  lines(0:nt,apf.u[,k+2],col=4,lty=2)
  lines(0:nt,kd.m[,k+2],col=3)
  lines(0:nt,kd.l[,k+2],col=3,lty=2)
  lines(0:nt,kd.u[,k+2],col=3,lty=2)
}
dev.off()

# Plot medians and upper/lower quantiles of %infected/susceptible over time
ylabs = expression(i,s)
xlabs = c("Time (days)","")
pdf(paste(gpath,"PF-states-med",param,"-",n,"-",N,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2),mar=c(5,5,4,1)+0.1)
for(k in 1:2)
{
  if(k == 1) {
    ymax = max(bf.u[,k],apf.u[,k],kd.u[,k],truth.u[,k])
    ymin = 0
  } else {
    ymax = 1
    ymin = min(bf.l[,k],apf.l[,k],kd.l[,k],truth.l[,k])
  }
  plot(0:nt,truth.m[,k],type="l",ylim=c(ymin,ymax),
  	xlab=xlabs[k],ylab=ylabs[k],cex.lab=2)
  lines(0:nt,truth.l[,k],lty=2)
  lines(0:nt,truth.u[,k],lty=2)
  lines(0:nt,bf.m[,k],col=2)
  lines(0:nt,bf.l[,k],col=2,lty=2)
  lines(0:nt,bf.u[,k],col=2,lty=2)
  lines(0:nt,apf.m[,k],col=4)
  lines(0:nt,apf.l[,k],col=4,lty=2)
  lines(0:nt,apf.u[,k],col=4,lty=2)
  lines(0:nt,kd.m[,k],col=3)
  lines(0:nt,kd.l[,k],col=3,lty=2)
  lines(0:nt,kd.u[,k],col=3,lty=2)
  if(k == 1 & n == 100) legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2))
}
dev.off()