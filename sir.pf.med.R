load("../Data/sir.pf.test.rdata") # Use same simulated data

N = 100
n = 100
bf.med = apf.med = kd.med = matrix(NA,nr=N,nc=tt)

for(j in 1:N)
{
  # Run bootstrap filter
  source("bf.R")
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

  for(i in 1:tt){
	bf.med[j,i] = wtd.quantile(out$state[1,,i], out$weight[,i], normwt=T, probs=.5)
	apf.med[j,i] = wtd.quantile(out2$state[1,,i], out2$weight[,i], normwt=T, probs=.5)
	kd.med[j,i] = wtd.quantile(out3$state[1,,i], out3$weight[,i], normwt=T, probs=.5)
  }
  print(j)
}

bf.i = apply(bf.med,2,mean)
apf.i = apply(apf.med,2,mean)
kd.i = apply(kd.med,2,mean)
bf.li = bf.ui = apf.li = apf.ui = kd.li = kd.ui = numeric(tt)
for(i in 1:tt){
	bf.li[i] = quantile(bf.med[,i],probs=.025)
	bf.ui[i] = quantile(bf.med[,i],probs=.975)
	apf.li[i] = quantile(apf.med[,i],probs=.025)
	apf.ui[i] = quantile(apf.med[,i],probs=.975)
	kd.li[i] = quantile(kd.med[,i],probs=.025)
	kd.ui[i] = quantile(kd.med[,i],probs=.975)
}
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

save.image("../Data/sir.pf.med.rdata")