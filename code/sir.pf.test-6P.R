# Set known parameter values
P = 5000
d = 5

# Set unknown parameter values
theta = c(0.2399, 0.1066, 1.2042, 0.25, 1.07, 0.050)

# Set prior bounds on unknown parameters
thetal = c(0.140, 0.090, 0.950, 0.05, 0.85, 0.0050)
thetau = c(0.500, 0.143, 1.300, 0.45, 1.15, 0.1000)

# Simulate epidemic
source("sir_functions.R")
revo_sim = function(x){ revo(x, P, d)}
robs_sim = function(x){ robs(x, theta[4], theta[5], theta[6])}
rinit_sim = function(){ rinit(10/P, theta)}
nt = 125
source("ss.sim.R")
sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)

# Run bootstrap filter
dllik_bf = function(y, x){ dllik(y, x, NULL, NULL, NULL, addparam=TRUE)}
revo_bf = function(x){ revo(x, P, d)}
rprior_bf = function(){ rprior(sim, thetal, thetau, addparam=TRUE)}
source("bf.R")
n = 10000
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
save.image(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/sir.pf.test-6P-",n,".rdata",
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
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/PF-quant-6P-",n,".pdf",sep=""))
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
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Hist-BF-6P-",n,"-day",cutoff-1,".pdf",sep=""),width=10,height=8)
par(mfrow=c(2,3))
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
weighted.hist(out$state[6,,cutoff],out$weight[,cutoff],
	xlab="b",main=expression(paste("Histogram of ",b,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste("b"," = ",aa,sep=""),list(aa=theta[4])),side=3)
weighted.hist(out$state[7,,cutoff],out$weight[,cutoff],
	xlab=expression(varsigma),main=expression(paste("Histogram of ",varsigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=theta[5])),side=3)
weighted.hist(out$state[8,,cutoff],out$weight[,cutoff],
	xlab=expression(sigma),main=expression(paste("Histogram of ",sigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=theta[6])),side=3)
dev.off()

# auxillary particle filter
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Hist-APF-6P-",n,"-day",cutoff-1,".pdf",sep=""),width=10,height=8)
par(mfrow=c(2,3))
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
weighted.hist(out2$state[6,,cutoff],out2$weight[,cutoff],
	xlab="b",main=expression(paste("Histogram of ",b,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste("b"," = ",aa,sep=""),list(aa=theta[4])),side=3)
weighted.hist(out2$state[7,,cutoff],out2$weight[,cutoff],
	xlab=expression(varsigma),main=expression(paste("Histogram of ",varsigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=theta[5])),side=3)
weighted.hist(out2$state[8,,cutoff],out2$weight[,cutoff],
	xlab=expression(sigma),main=expression(paste("Histogram of ",sigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=theta[6])),side=3)
dev.off()

# kernel density particle filter
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Hist-KD-6P-",n,"-day",cutoff-1,".pdf",sep=""),width=10,height=8)
par(mfrow=c(2,3))
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
weighted.hist(theta2u(out3$theta[4,,cutoff],thetal[4],thetau[4]),out3$weight[,cutoff],
	xlab="b",main=expression(paste("Histogram of ",b,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste("b"," = ",aa,sep=""),list(aa=theta[4])),side=3)
weighted.hist(theta2u(out3$theta[5,,cutoff],thetal[5],thetau[5]),out3$weight[,cutoff],
	xlab=expression(varsigma),main=expression(paste("Histogram of ",varsigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=theta[5])),side=3)
weighted.hist(theta2u(out3$theta[6,,cutoff],thetal[6],thetau[6]),out3$weight[,cutoff],
	xlab=expression(sigma),main=expression(paste("Histogram of ",sigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=theta[6])),side=3)
dev.off()