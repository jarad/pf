# Set known parameter values
P = 5000
d = 5
b = c(.25, .27, .23, .29)
varsigma = c(1.07, 1.05, 1.01, .98)
sigma = 40*c(.0012, .0008, .0010, .0011)

# Set unknown parameter values
theta = c(0.2399, 0.1066, 1.2042)

# Set prior bounds on unknown parameters
thetal = c(0.140, 0.090, 0.950)
thetau = c(0.500, 0.143, 1.300)

# Simulate epidemic
source("sir_functions.R")
revo_sim = function(x){ revo(x, P, d)}
robs_sim = function(x){ robs(x, b, varsigma, sigma)}
rinit_sim = function(){ rinit(10/P, theta)}
nt = 125
source("ss.sim.R")
sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)

# Run bootstrap filter
dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma)}
revo_bf = function(x){ revo(x, P, d)}
rprior_bf = function(){ rprior(sim, thetal, thetau, b, varsigma, sigma)}
source("bf.R")
n = 10000
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

# Save data
save.image(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/sir.pf.test-",n,".rdata",
	sep=""))

# Comparison
# Plot % infected
require(Hmisc)
tt = nt + 1
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
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/PF-quant-",n,".pdf",sep=""))
plot(sim$x[1,],type="l",ylim=c(0,.28),xlab="Time (days)",ylab="% Population",
	main="95% Credible Intervals of %Pop Infected")
lines(bf.i,col=2)
lines(bf.li,col=2,lty=2)
lines(bf.ui,col=2,lty=2)
lines(apf.i,col=4)
lines(apf.li,col=4,lty=2)
lines(apf.ui,col=4,lty=2)
lines(kd.i,col=3)
lines(kd.li,col=3,lty=2)
lines(kd.ui,col=3,lty=2)
legend("topright",c("Truth","BF","APF","KD","95% bounds"),col=c(1,2,4,3,1),lty=c(1,1,1,1,2))
dev.off()

# Histograms of unknown parameters
# bootstrap filter
require(plotrix)
cutoff = nt + 1
msize = labsize = 1.5
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Hist-BF-",n,"-day",cutoff-1,".pdf",sep=""),width=5,height=8)
par(mfrow=c(3,1)) 
weighted.hist(out$state[3,,cutoff],out$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=theta[1])),side=3)
weighted.hist(out$state[4,,cutoff],out$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=theta[2])),side=3)
weighted.hist(out$state[5,,cutoff],out$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=theta[3])),side=3)
dev.off()

# auxillary particle filter
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Hist-APF-",n,"-day",cutoff-1,".pdf",sep=""),width=5,height=8)
par(mfrow=c(3,1))
weighted.hist(out2$state[3,,cutoff],out2$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=theta[1])),side=3)
weighted.hist(out2$state[4,,cutoff],out2$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=theta[2])),side=3)
weighted.hist(out2$state[5,,cutoff],out2$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=theta[3])),side=3)
dev.off()

# kernel density particle filter
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Hist-KD-",n,"-day",cutoff-1,".pdf",sep=""),width=5,height=8)
par(mfrow=c(3,1))
weighted.hist(theta2u(out3$theta[1,,cutoff],thetal[1],thetau[1]),out3$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=theta[1])),side=3)
weighted.hist(theta2u(out3$theta[2,,cutoff],thetal[2],thetau[2]),out3$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=theta[2])),side=3)
weighted.hist(theta2u(out3$theta[3,,cutoff],thetal[3],thetau[3]),out3$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=theta[3])),side=3)
dev.off()
