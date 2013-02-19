# Run particle filters with one of beta, gamma, or nu unknown.
# Uses same simulated data generated from sir.pf.R - make sure data and graphics path set to where desired data sets are and graphs go.

# Set graphics and data paths
gpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/"
dpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/"

# Load data
load(paste(dpath,"sir.pf-100.rdata",sep=""))

# Which parameter unknown?
p = 2
# Label for data and graphs
param = "-1Pgamma"
# How many particles?
n = 10000

# Create parameter-specific variables
s = rep(0,3); s[p] = 1; t = rep(1,3); t[p] = 0

# Run bootstrap filter
dllik_bf = function(y, x){ dllik(y, c(x[1], x[2], x[3]*s + theta*t), b, varsigma, sigma)}
revo_bf = function(x){ revo(c(x[1], x[2], x[3]*s + theta*t), P, d)[c(1,2,2+p)]}
rprior_bf = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma)[c(1,2,2+p)]}
source("bf.R")
out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run auxiliary particle filter
dllik_apf = function(y, x){ dllik(y, c(x[1], x[2], x[3]*s + theta*t), b, varsigma, sigma)}
pstate_apf = function(x) { revo(c(x[1], x[2], x[3]*s + theta*t), P, d, random = FALSE)[c(1,2,2+p)]}
revo_apf = function(x){ revo(c(x[1], x[2], x[3]*s + theta*t), P, d)[c(1,2,2+p)]}
rprior_apf = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma)[c(1,2,2+p)]}
source("apf.R")
out2 = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run kernel density particle filter
dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma)}
pstate_kd = function(x, mytheta) { revo(x, P, d, mytheta*s + u2theta(theta,thetal,thetau)*t, thetal, thetau, FALSE, FALSE)}
revo_kd = function(x, mytheta){ revo(x, P, d, mytheta*s + u2theta(theta,thetal,thetau)*t, thetal, thetau, FALSE)}
rprior_kd = function()
{ 
  mylist = rprior(sim, thetal, thetau, b, varsigma, sigma, FALSE)
  return(list(x=mylist$x,theta=mylist$theta[p]))
}
source("kd_pf.R")
out3 = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Save data
save.image(paste(dpath,"sir.pf",param,"-",n,".rdata",sep=""))
#load(paste(dpath,"sir.pf",param,"-",n,".rdata",sep=""))

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
bf.m = apf.m = kd.m = numeric(tt)
bf.l = bf.u = apf.l = apf.u = kd.l = kd.u = numeric(tt)
for(i in 1:tt)
{
  bf.m[i] = wtd.quantile(out$state[3,,i],out$weight[,i],normwt=T,probs=.5)
  apf.m[i] = wtd.quantile(out2$state[3,,i],out$weight[,i],normwt=T,probs=.5)
  kd.m[i] = wtd.quantile(theta2u(out3$theta[1,,i],thetal[p],thetau[p]),out$weight[,i],normwt=T,probs=.5)
  bf.l[i] = wtd.quantile(out$state[3,,i],out$weight[,i],normwt=T,probs=.025)
  apf.l[i] = wtd.quantile(out2$state[3,,i],out$weight[,i],normwt=T,probs=.025)
  kd.l[i] = wtd.quantile(theta2u(out3$theta[1,,i],thetal[p],thetau[p]),out$weight[,i],normwt=T,probs=.025)
  bf.u[i] = wtd.quantile(out$state[3,,i],out$weight[,i],normwt=T,probs=.975)
  apf.u[i] = wtd.quantile(out2$state[3,,i],out$weight[,i],normwt=T,probs=.975)
  kd.u[i] = wtd.quantile(theta2u(out3$theta[1,,i],thetal[p],thetau[p]),out$weight[,i],normwt=T,probs=.975)
}

# Save quantiles data
save.image(paste(dpath,"sir.pf.quant",param,"-",n,".rdata",sep=""))
#load(paste(dpath,"sir.pf.quant",param,"-",n,".rdata",sep=""))

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
expr = expression(beta,gamma,nu,b,varsigma,sigma)
msize = 1.75
labsize = 1.5
legendsize = 1.0
pdf(paste(gpath,"PF-params",param,"-",n,".pdf",sep=""))
par(mar=c(5,5,3,1)+.1)
ymin = min(bf.l,apf.l,kd.l,theta[p])
ymax = max(bf.u,apf.u,kd.u,theta[p])
if(n == 100)
{
  plot(0:nt,bf.m,type="l",ylim=c(ymin,ymax),col=2,xlab="Time (days)",ylab=paste("J = ",n,sep=""),main=expr[p],cex.lab=labsize,cex.main=msize)
  legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c (1,1,1,1,2),cex=legendsize)
} else {
  plot(0:nt,bf.m,type="l",ylim=c(ymin,ymax),col=2,xlab="",ylab=paste("J = ",n,sep=""),main=expr[p],cex.lab=labsize,cex.main=msize)
}
lines(0:nt,bf.l,col=2,lty=2)
lines(0:nt,bf.u,col=2,lty=2)
lines(0:nt,apf.m,col=4)
lines(0:nt,apf.l,col=4,lty=2)
lines(0:nt,apf.u,col=4,lty=2)
lines(0:nt,kd.m,col=3)
lines(0:nt,kd.l,col=3,lty=2)
lines(0:nt,kd.u,col=3,lty=2)
abline(h=theta[p])
dev.off()