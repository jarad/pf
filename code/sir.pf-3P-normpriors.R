##### Set true values and simulate data #####
# Skip to (*) if loading previous data instead
#############################################

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

# Simulate epidemic
source("sir_functions.R")
revo_sim = function(x){ revo(x, P, d, theta)}
robs_sim = function(x){ robs(x, b, varsigma, sigma, dpower)}
rinit_sim = function(){ rinit(10/P)}
nt = 125
source("ss.sim.R")
sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
save.image(paste(dpath,"sim-xy.rdata",sep=""))

###### (*) ######
dpath = "C:/Users/Danny/My Documents/UCSB - Research/pf/data/D2/"
load(paste(dpath,"sim-xy.rdata",sep=""))
################# 

# Plot the data
y = apply(sim$y,2,function(x) x[which(!is.na(x))])
J = apply(sim$y,2,function(x) which(!is.na(x)))
no = dim(sim$y)[1]
pdf(paste(gpath,"sim-y",param,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2),mar=c(5,5,4,1)+.1)
plot(which(J==1),y[J==1],ylim=c(min(y),max(y)),xlim=c(0,nt),
	xlab="Time (Days)",ylab=expression(y),cex.lab=1.5)
if(no > 1) for(i in 2:no) points(which(J==i),y[J==i],col=i)
plot(which(J==1),exp(y[J==1]),ylim=c(0,max(exp(y))),xlim=c(0,nt),
	xlab="Time (Days)",ylab=expression(e^y),cex.lab=1.5)
if(no > 1) for(i in 2:no) points(which(J==i),exp(y[J==i]),col=i)
dev.off()

# Create function to sample from prior distribution of theta
theta.mean = c(-1.3296, -2.1764, 0.1055)
theta.sd = sqrt(c(0.1055, 0.0140, 0.0064))
rtheta = function(){ rnorm(3,theta.mean,theta.sd)}

# How many particles?
n = 1000

# Run bootstrap filter
dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
revo_bf = function(x){ revo(x, P, d, exp(x[3:5]))}
rprior_bf = function()
{ 
  myprior = rprior(sim$y[,1], rtheta, b, varsigma, sigma, dpower)
  return(c(myprior$x,myprior$theta))
}
source("bf.R")
out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n, method="stratified", nonuniformity="ess", threshold=0.8, log=F)

# Run auxiliary particle filter
dllik_apf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
pstate_apf = function(x) { revo(x, P, d, exp(x[3:5]), FALSE)}
revo_apf = function(x){ revo(x, P, d, exp(x[3:5]))}
rprior_apf = function()
{ 
  myprior = rprior(sim$y[,1], rtheta, b, varsigma, sigma, dpower)
  return(c(myprior$x,myprior$theta))
}
source("apf.R")
out2 = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, method="stratified", nonuniformity="ess", threshold=0.8, log=F)

# Run kernel density particle filter
dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma, dpower)}
pstate_kd = function(x, theta) { revo(x, P, d, exp(theta), FALSE)}
revo_kd = function(x, theta){ revo(x, P, d, exp(theta))}
rprior_kd = function(){ rprior(sim$y[,1], rtheta, b, varsigma, sigma, dpower)}
source("kd_pf.R")
out3 = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, method="stratified", nonuniformity="ess", threshold=0.8, log=F)

# Calculate 2.5%, 50%, and 97.5% quantiles of states over time
require(Hmisc)
tt = nt + 1
bf.mi = apf.mi = kd.mi = bf.ms = apf.ms = kd.ms = numeric(tt)
bf.li = bf.ui = apf.li = apf.ui = kd.li = kd.ui = numeric(tt)
bf.ls = bf.us = apf.ls = apf.us = kd.ls = kd.us = numeric(tt)
for(i in 1:tt)
{
  bf.li[i] = wtd.quantile(out$state[1,,i], out$weight[,i], normwt=T, probs=.025)
  bf.ui[i] = wtd.quantile(out$state[1,,i], out$weight[,i], normwt=T, probs=.975)
  bf.mi[i] = wtd.quantile(out$state[1,,i], out$weight[,i], normwt=T, probs=.500)
  bf.ls[i] = wtd.quantile(out$state[2,,i], out$weight[,i], normwt=T, probs=.025)
  bf.us[i] = wtd.quantile(out$state[2,,i], out$weight[,i], normwt=T, probs=.975)
  bf.ms[i] = wtd.quantile(out$state[2,,i], out$weight[,i], normwt=T, probs=.500)
  apf.li[i] = wtd.quantile(out2$state[1,,i], out$weight[,i], normwt=T, probs=.025)
  apf.ui[i] = wtd.quantile(out2$state[1,,i], out$weight[,i], normwt=T, probs=.975)
  apf.mi[i] = wtd.quantile(out2$state[1,,i], out$weight[,i], normwt=T, probs=.500)
  apf.ls[i] = wtd.quantile(out2$state[2,,i], out$weight[,i], normwt=T, probs=.025)
  apf.us[i] = wtd.quantile(out2$state[2,,i], out$weight[,i], normwt=T, probs=.975)
  apf.ms[i] = wtd.quantile(out2$state[2,,i], out$weight[,i], normwt=T, probs=.500)
  kd.li[i] = wtd.quantile(out3$state[1,,i], out$weight[,i], normwt=T, probs=.025)
  kd.ui[i] = wtd.quantile(out3$state[1,,i], out$weight[,i], normwt=T, probs=.975)
  kd.mi[i] = wtd.quantile(out3$state[1,,i], out$weight[,i], normwt=T, probs=.500)
  kd.ls[i] = wtd.quantile(out3$state[2,,i], out$weight[,i], normwt=T, probs=.025)
  kd.us[i] = wtd.quantile(out3$state[2,,i], out$weight[,i], normwt=T, probs=.975)
  kd.ms[i] = wtd.quantile(out3$state[2,,i], out$weight[,i], normwt=T, probs=.500)
}

# Calculate 2.5%, 50%, and 97.5% quantiles of unknown parameters over time
bf.m = apf.m = kd.m = matrix(NA,nr=tt,nc=3)
bf.l = bf.u = apf.l = apf.u = kd.l = kd.u = matrix(NA,nr=tt,nc=3)
st = 3:5
for(i in 1:tt)
{
  for(j in 1:3)
  {
    bf.m[i,j] = wtd.quantile(exp(out$state[st[j],,i]),out$weight[,i],normwt=T,probs=.5)
    apf.m[i,j] = wtd.quantile(exp(out2$state[st[j],,i]),out$weight[,i],normwt=T,probs=.5)
    kd.m[i,j] = wtd.quantile(exp(out3$theta[j,,i]),out$weight[,i],normwt=T,probs=.5)
    bf.l[i,j] = wtd.quantile(exp(out$state[st[j],,i]),out$weight[,i],normwt=T,probs=.025)
    apf.l[i,j] = wtd.quantile(exp(out2$state[st[j],,i]),out$weight[,i],normwt=T,probs=.025)
    kd.l[i,j] = wtd.quantile(exp(out3$theta[j,,i]),out$weight[,i],normwt=T,probs=.025)
    bf.u[i,j] = wtd.quantile(exp(out$state[st[j],,i]),out$weight[,i],normwt=T,probs=.975)
    apf.u[i,j] = wtd.quantile(exp(out2$state[st[j],,i]),out$weight[,i],normwt=T,probs=.975)
    kd.u[i,j] = wtd.quantile(exp(out3$theta[j,,i]),out$weight[,i],normwt=T,probs=.975)
  }
}

# Save data
save.image(paste(dpath,"sir.pf",param,"-",n,".rdata",sep=""))
#load(paste(dpath,"sir.pf",param,"-",n,".rdata",sep=""))

# Plot % infected and % susceptible quantiles over time
pdf(paste(gpath,"PF-states",param,"-",n,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2),mar=c(5,5,3,1)+0.1)
ymax = max(bf.ui[-1],apf.ui[-1],kd.ui[-1],sim$x[1,])
if(n == 100)
{
  plot(0:nt,sim$x[1,],type="l",ylim=c(0,ymax),xlab="",ylab=paste("J = ",n,sep=""),main=expression(i),cex.lab=1.5,cex.main=1.75)
  legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2))
} else if(n == 1000){
  plot(0:nt,sim$x[1,],type="l",ylim=c(0,ymax),xlab="",ylab=paste("J = ",n,sep=""),cex.lab=1.5)
} else {
  plot(0:nt,sim$x[1,],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab=paste("J = ",n,sep=""),cex.lab=1.5)
}
lines(0:nt,bf.mi,col=2)
lines(0:nt,bf.li,col=2,lty=2)
lines(0:nt,bf.ui,col=2,lty=2)
lines(0:nt,apf.mi,col=4)
lines(0:nt,apf.li,col=4,lty=2)
lines(0:nt,apf.ui,col=4,lty=2)
lines(0:nt,kd.mi,col=3)
lines(0:nt,kd.li,col=3,lty=2)
lines(0:nt,kd.ui,col=3,lty=2)
ymin = min(bf.ls[-1],apf.ls[-1],kd.ls[-1],sim$x[2,])
if(n == 100)
{
  plot(0:nt,sim$x[2,],type="l",ylim=c(ymin,1),xlab="",ylab="",main=expression(s),cex.main=1.75)
} else {
  plot(0:nt,sim$x[2,],type="l",ylim=c(ymin,1),xlab="",ylab="")
}
lines(0:nt,bf.ms,col=2)
lines(0:nt,bf.ls,col=2,lty=2)
lines(0:nt,bf.us,col=2,lty=2)
lines(0:nt,apf.ms,col=4)
lines(0:nt,apf.ls,col=4,lty=2)
lines(0:nt,apf.us,col=4,lty=2)
lines(0:nt,kd.ms,col=3)
lines(0:nt,kd.ls,col=3,lty=2)
lines(0:nt,kd.us,col=3,lty=2)
dev.off()

# Plot 95% credible bounds and medians of parameters over time
expr = expression(beta,gamma,nu)
xlabs = c("Time (days)","","")
ylabs = c(paste("J = ",n,sep=""),"","")
msize = 2.4
labsize = 2.2
legendsize = 1.5
pdf(paste(gpath,"PF-params",param,"-",n,".pdf",sep=""),width=15,height=5)
par(mfrow=c(1,3),mar=c(5,5,3,1)+.1)
for(j in 1:3)
{
  ymin = min(bf.l[,j],apf.l[,j],kd.l[,j],theta[j])
  ymax = max(bf.u[,j],apf.u[,j],kd.u[,j],theta[j])
  if(n == 100)
  {
    plot(0:nt,bf.m[,j],type="l",ylim=c(ymin,ymax),col=2,,xlab="",ylab=ylabs[j],main=expr[j],cex.lab=labsize,cex.main=msize)
    if(j == 1) legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2),cex=legendsize)
  } else if(n == 10000) {
    plot(0:nt,bf.m[,j],type="l",ylim=c(ymin,ymax),col=2,xlab=xlabs[j],ylab=ylabs[j],cex.lab=labsize)
  } else {
    plot(0:nt,bf.m[,j],type="l",ylim=c(ymin,ymax),col=2,ylab=ylabs[j],xlab="",cex.lab=labsize)
  }
  lines(0:nt,bf.l[,j],col=2,lty=2)
  lines(0:nt,bf.u[,j],col=2,lty=2)
  lines(0:nt,apf.m[,j],col=4)
  lines(0:nt,apf.l[,j],col=4,lty=2)
  lines(0:nt,apf.u[,j],col=4,lty=2)
  lines(0:nt,kd.m[,j],col=3)
  lines(0:nt,kd.l[,j],col=3,lty=2)
  lines(0:nt,kd.u[,j],col=3,lty=2)
  abline(h=theta[j])
}
dev.off()