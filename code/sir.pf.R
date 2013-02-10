# Initialize .Random.seed
set.seed(sample(1:1000,1))

# Set graphics and data paths
gpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/"
dpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/"

# How many unknown parameters? Set p = 3 or p = 6
p = 6
if(p == 3)
{
  s = 4
  param = ""
} else { 
  s = 1
  param = "-6P"
}

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

# Simulate epidemic
source("sir_functions.R")
revo_sim = function(x){ revo(x, P, d)}
robs_sim = function(x){ robs(x, b[1:s], varsigma[1:s], sigma[1:s])}
rinit_sim = function(){ rinit(10/P, theta)}
nt = 125
source("ss.sim.R")
sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)

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

# How many particles?
n = 100

# 3 parameter PFs
# Run bootstrap filter
dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma)}
revo_bf = function(x){ revo(x, P, d)}
rprior_bf = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma)}
source("bf.R")
out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run auxiliary particle filter
dllik_apf = function(y, x){ dllik(y, x, b, varsigma, sigma)}
pstate_apf = function(x) { revo(x, P, d, random = FALSE)}
revo_apf = function(x){ revo(x, P, d)}
rprior_apf = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma)}
source("apf.R")
out2 = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run kernel density particle filter
dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma)}
pstate_kd = function(x, theta) { revo(x, P, d, theta, thetal, thetau, FALSE, FALSE)}
revo_kd = function(x, theta){ revo(x, P, d, theta, thetal, thetau, FALSE)}
rprior_kd = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma, FALSE)}
source("kd_pf.R")
out3 = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# 6 parameter PFs
# Run bootstrap filter
dllik_bf = function(y, x){ dllik(y, x, NULL, NULL, NULL, addparam=TRUE)}
revo_bf = function(x){ revo(x, P, d)}
rprior_bf = function(){ rprior(sim, thetal, thetau, addparam=TRUE)}
source("bf.R")
out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run auxiliary particle filter
dllik_apf = function(y, x){ dllik(y, x, NULL, NULL, NULL, addparam=TRUE)}
pstate_apf = function(x) { revo(x, P, d, random = FALSE)}
revo_apf = function(x){ revo(x, P, d)}
rprior_apf = function(){ rprior(sim, thetal, thetau, addparam=TRUE)}
source("apf.R")
out2 = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run kernel density particle filter
dllik_kd = function(y, x, theta){ dllik(y, x, NULL, NULL, NULL, theta, thetal, thetau, FALSE, TRUE)}
pstate_kd = function(x, theta) { revo(x, P, d, theta, thetal, thetau, FALSE, FALSE)}
revo_kd = function(x, theta){ revo(x, P, d, theta, thetal, thetau, FALSE)}
rprior_kd = function(){ rprior(sim, thetal, thetau, stateonly=FALSE, addparam=TRUE)}
source("kd_pf.R")
out3 = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Save data
save.image(paste(dpath,"sir.pf",param,"-",n,".rdata",sep=""))
#load(paste(dpath,"sir.pf-",n,".rdata",sep=""))

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

# Plot % infected and % susceptible quantiles over time
pdf(paste(gpath,"PF-states",param,"-",n,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2))
ymax = max(bf.ui[-1],apf.ui[-1],kd.ui[-1],sim$x[1,])
plot(0:nt,sim$x[1,],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab="% Population",
	main="95% Credible Intervals of %Pop Infected")
lines(0:nt,bf.mi,col=2)
lines(0:nt,bf.li,col=2,lty=2)
lines(0:nt,bf.ui,col=2,lty=2)
lines(0:nt,apf.mi,col=4)
lines(0:nt,apf.li,col=4,lty=2)
lines(0:nt,apf.ui,col=4,lty=2)
lines(0:nt,kd.mi,col=3)
lines(0:nt,kd.li,col=3,lty=2)
lines(0:nt,kd.ui,col=3,lty=2)
if(n == 100) legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2))
ymin = min(bf.ls[-1],apf.ls[-1],kd.ls[-1],sim$x[2,])
plot(0:nt,sim$x[2,],type="l",ylim=c(ymin,1),xlab="Time (days)",ylab="% Population",
	main="95% Credible Intervals of %Pop Susceptible")
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

# Calculate 2.5%, 50%, and 97.5% quantiles of unknown parameters over time
bf.m = apf.m = kd.m = matrix(NA,nr=tt,nc=p)
bf.l = bf.u = apf.l = apf.u = kd.l = kd.u = matrix(NA,nr=tt,nc=p)
st = 3:8
for(i in 1:tt)
{
  for(j in 1:p)
  {
    bf.m[i,j] = wtd.quantile(out$state[st[j],,i],out$weight[,i],normwt=T,probs=.5)
    apf.m[i,j] = wtd.quantile(out2$state[st[j],,i],out$weight[,i],normwt=T,probs=.5)
    kd.m[i,j] = wtd.quantile(theta2u(out3$theta[j,,i],thetal[j],thetau[j]),out$weight[,i],normwt=T,probs=.5)
    bf.l[i,j] = wtd.quantile(out$state[st[j],,i],out$weight[,i],normwt=T,probs=.025)
    apf.l[i,j] = wtd.quantile(out2$state[st[j],,i],out$weight[,i],normwt=T,probs=.025)
    kd.l[i,j] = wtd.quantile(theta2u(out3$theta[j,,i],thetal[j],thetau[j]),out$weight[,i],normwt=T,probs=.025)
    bf.u[i,j] = wtd.quantile(out$state[st[j],,i],out$weight[,i],normwt=T,probs=.975)
    apf.u[i,j] = wtd.quantile(out2$state[st[j],,i],out$weight[,i],normwt=T,probs=.975)
    kd.u[i,j] = wtd.quantile(theta2u(out3$theta[j,,i],thetal[j],thetau[j]),out$weight[,i],normwt=T,probs=.975)
  }
}

# Plot 95% credible bounds and medians of parameters over time
expr = expression(beta,gamma,nu,b,varsigma,sigma)
pdf(paste(gpath,"PF-params",param,"-",n,".pdf",sep=""),width=15,height=5*((p == 6)+1))
par(mfrow=c((p == 6)+1,3),mar=c(5,5,4,1)+.1)
for(j in 1:p)
{
  ymin = min(bf.l[,j],apf.l[,j],kd.l[,j],theta[j])
  ymax = max(bf.u[,j],apf.u[,j],kd.u[,j],theta[j])
  plot(0:nt,bf.m[,j],type="l",xlab="Time (Days)",ylab=expr[j],ylim=c(ymin,ymax),col=2,
  	cex.lab=1.75)
  lines(0:nt,bf.l[,j],col=2,lty=2)
  lines(0:nt,bf.u[,j],col=2,lty=2)
  lines(0:nt,apf.m[,j],col=4)
  lines(0:nt,apf.l[,j],col=4,lty=2)
  lines(0:nt,apf.u[,j],col=4,lty=2)
  lines(0:nt,kd.m[,j],col=3)
  lines(0:nt,kd.l[,j],col=3,lty=2)
  lines(0:nt,kd.u[,j],col=3,lty=2)
  abline(h=theta[j])
  if(n == 100 & j == 1) legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2))
  par(mar=c(5,5,4,1)+.1)
}
dev.off()

# Save quantiles data
save.image(paste(dpath,"sir.pf.quant",param,"-",n,".rdata",sep=""))