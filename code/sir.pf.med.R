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

# Functions to simulate epidemic
source("sir_functions.R")
revo_sim = function(x){ revo(x, P, d)}
robs_sim = function(x){ robs(x, b, varsigma, sigma)}
rinit_sim = function(){ rinit(10/P, theta)}

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

# Run particle filters N times
N = 100
n = 100
nt = 125
bf.med = apf.med = kd.med = matrix(NA,nr=N,nc=nt+1)
source("ss.sim.R")
source("bf.R")
source("apf.R")
source("kd_pf.R")
require(Hmisc)
for(j in 1:N)
{
  # Simulate epidemic
  sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
 
  # Run bootstrap filter
  out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

  # Run auxilliary particle filter
  out2 = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

  # Run kernel density particle filter
  out3 = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n,
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

  for(i in 1:(nt+1)){
	bf.med[j,i] = wtd.quantile(out$state[1,,i], out$weight[,i], normwt=T, probs=.5)
	apf.med[j,i] = wtd.quantile(out2$state[1,,i], out2$weight[,i], normwt=T, probs=.5)
	kd.med[j,i] = wtd.quantile(out3$state[1,,i], out3$weight[,i], normwt=T, probs=.5)
  }
  cat("\n",j,"\n")
}

bf.i = apply(bf.med,2,mean)
apf.i = apply(apf.med,2,mean)
kd.i = apply(kd.med,2,mean)
bf.li = bf.ui = apf.li = apf.ui = kd.li = kd.ui = numeric(nt+1)
for(i in 1:(nt+1)){
	bf.li[i] = quantile(bf.med[,i],probs=.025)
	bf.ui[i] = quantile(bf.med[,i],probs=.975)
	apf.li[i] = quantile(apf.med[,i],probs=.025)
	apf.ui[i] = quantile(apf.med[,i],probs=.975)
	kd.li[i] = quantile(kd.med[,i],probs=.025)
	kd.ui[i] = quantile(kd.med[,i],probs=.975)
}
ymax = max(bf.ui,apf.ui,kd.ui)
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/PF-med-",n,"-",n,".pdf",sep=""))
plot(sim$x[1,],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab="% Population",
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

save.image(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/sir.pf.med-",n,"-",n,".rdata",sep=""))