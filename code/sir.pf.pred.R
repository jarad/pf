# Load particle filter data
param = "-6P"
num = "-10000"
load(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/sir.pf.test",param,num,".rdata",sep=""))
ns = dim(sim$x)[1]

# Function needed to generate predicted SIR curves
ss.pred = function(nt,revo,rinit)
{
  # Find dimension of state
  current.seed = .Random.seed
  ns = length(rinit())
  .Random.seed = current.seed

  # Initialize state matrix
  x = matrix(NA,ns,nt+1)
  x[,1] = rinit()

  # Generate states and observations
  for(i in 1:nt)
  {
    x[,i+1] = revo(x[,i])
  }

  return(x)
}

# Predict epidemic by propagating particles after being filtered for different numbers of days
require(Hmisc)
days <- seq(10,35,5)
bf.li <- bf.i <- bf.ui <- list()
apf.li <- apf.i <- apf.ui <- list()
kd.li <- kd.i <- kd.ui <- list()
bf.peak.prop <- bf.peak.calc <- list()
apf.peak.prop <- apf.peak.calc <- list()
kd.peak.prop <- kd.peak.calc <- list()
length(bf.li) = length(bf.i) = length(bf.ui) = length(days)
length(apf.li) = length(apf.i) = length(apf.ui) = length(days)
length(kd.li) = length(kd.i) = length(kd.ui) = length(days)
length(bf.peak.prop) = length(bf.peak.calc) = length(days)
length(apf.peak.prop) = length(apf.peak.calc) = length(days)
length(kd.peak.prop) = length(kd.peak.calc) = length(days)
pb = txtProgressBar(0,n*length(days),style=3)
for(j in 1:length(days))
{
  # Grab  particles from each filter at time days[j]
  mystate = out$state[,,days[j]+1]
  mystate2 = out2$state[,,days[j]+1]
  mybeta = theta2u(out3$theta[1,,days[j]+1],thetal[1],thetau[1])
  mygamma = theta2u(out3$theta[2,,days[j]+1],thetal[2],thetau[2])
  mynu = theta2u(out3$theta[3,,days[j]+1],thetal[3],thetau[3])
  if(ns == 8) {
    myb = theta2u(out3$theta[4,,days[j]+1],thetal[4],thetau[4])
    myvarsigma = theta2u(out3$theta[5,,days[j]+1],thetal[5],thetau[5])
    mysigma = theta2u(out3$theta[6,,days[j]+1],thetal[6],thetau[6])
    mystate3 = rbind(out3$state[,,days[j]+1],mybeta,mygamma,mynu,myb,myvarsigma,mysigma)
  } else if(ns == 5) {
    mystate3 = rbind(out3$state[,,days[j]+1],mybeta,mygamma,mynu)  
  } else { stop("Wrong number of parameters") }

  # Propagate particles forward to end of epidemic
  pred = pred2 = pred3 = array(NA,dim=c(8,n,nt-days[j]))
  for(i in 1:n)
  {
    setTxtProgressBar(pb,i + n*(j-1))
    pred[,i,] = ss.pred(nt-(days[j]+1),pstate_apf,function() return(mystate[,i]))
    pred2[,i,] = ss.pred(nt-(days[j]+1),pstate_apf,function() return(mystate2[,i]))
    pred3[,i,] = ss.pred(nt-(days[j]+1),pstate_apf,function() return(mystate3[,i]))
  }

  # Calculate .025 and .975 quantiles of predicted particle samples
  bf.li[[j]] = bf.i[[j]] = bf.ui[[j]] = apf.li[[j]] = apf.i[[j]] = apf.ui[[j]] = kd.li[[j]] = kd.i[[j]] = kd.ui[[j]] = rep(NA,nt-days[j])
  for(i in 1:(nt-days[j]))
  {
    bf.li[[j]][i] = wtd.quantile(pred[1,,i], out$weight[,days[j]+1], normwt=T, probs=.025)
    bf.i[[j]][i] = wtd.quantile(pred[1,,i], out$weight[,days[j]+1], normwt=T, probs=.5)
    bf.ui[[j]][i] = wtd.quantile(pred[1,,i], out$weight[,days[j]+1], normwt=T, probs=.975)
    apf.li[[j]][i] = wtd.quantile(pred2[1,,i], out2$weight[,days[j]+1], normwt=T, probs=.025)
    apf.i[[j]][i] = wtd.quantile(pred2[1,,i], out2$weight[,days[j]+1], normwt=T, probs=.5)
    apf.ui[[j]][i] = wtd.quantile(pred2[1,,i], out2$weight[,days[j]+1], normwt=T, probs=.975)
    kd.li[[j]][i] = wtd.quantile(pred3[1,,i], out3$weight[,days[j]+1], normwt=T, probs=.025)
    kd.i[[j]][i] = wtd.quantile(pred3[1,,i], out3$weight[,days[j]+1], normwt=T, probs=.5)
    kd.ui[[j]][i] = wtd.quantile(pred3[1,,i], out3$weight[,days[j]+1], normwt=T, probs=.975)
  }

  # Get peak times and intensities of particles
  # By propagation
  peak_time = matrix(NA,nr=n,nc=3)
  peak_intensity = matrix(NA,nr=n,nc=3)
  get.peak.time = function(x){which(x == max(x))}
  peak_intensity[,1] = apply(pred[1,,],1,max)
  peak_intensity[,2] = apply(pred2[1,,],1,max)
  peak_intensity[,3] = apply(pred3[1,,],1,max)
  peak_time[,1] = apply(pred[1,,],1,get.peak.time) - 1 + days[j] 
  peak_time[,2] = apply(pred2[1,,],1,get.peak.time) - 1 + days[j] 
  peak_time[,3] = apply(pred3[1,,],1,get.peak.time) - 1 + days[j] 
  # By calculation
  peak_time_calc = matrix(NA,nr=n,nc=3)
  peak_intensity_calc = matrix(NA,nr=n,nc=3)
  calc.peak.time = function(x) log(1/x[ns+1])/(x[4]*(x[3]/x[4] - 1)) + 0.4/x[4]
  calc.peak.intensity = function(x) 1 - x[4]/x[3] + log(x[4]/x[3])/(x[3]/x[4])
  peak_intensity_calc[,1] = apply(rbind(mystate,out$state[1,,1]),2,calc.peak.intensity)
  peak_intensity_calc[,2] = apply(rbind(mystate2,out2$state[1,,1]),2,calc.peak.intensity)
  peak_intensity_calc[,3] = apply(rbind(mystate3,out3$state[1,,1]),2,calc.peak.intensity)
  peak_time_calc[,1] = apply(rbind(mystate,out$state[1,,1]),2,calc.peak.time)
  peak_time_calc[,2] = apply(rbind(mystate2,out2$state[1,,1]),2,calc.peak.time)
  peak_time_calc[,3] = apply(rbind(mystate3,out3$state[1,,1]),2,calc.peak.time)

  # Calculate .025, .5, and .975 quantiles of peak times and intensities
  bf.peak.prop[[j]] = apf.peak.prop[[j]] = kd.peak.prop[[j]] = matrix(NA,nr=2,nc=3)
  bf.peak.prop[[j]][1,] = wtd.quantile(peak_time[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  apf.peak.prop[[j]][1,] = wtd.quantile(peak_time[,2],out2$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  kd.peak.prop[[j]][1,] = wtd.quantile(peak_time[,3],out3$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  bf.peak.prop[[j]][2,] = wtd.quantile(peak_intensity[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  apf.peak.prop[[j]][2,] = wtd.quantile(peak_intensity[,2],out2$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  kd.peak.prop[[j]][2,] = wtd.quantile(peak_intensity[,3],out3$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  
  bf.peak.calc[[j]] = apf.peak.calc[[j]] = kd.peak.calc[[j]] = matrix(NA,nr=2,nc=3)
  bf.peak.calc[[j]][1,] = wtd.quantile(peak_time_calc[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  apf.peak.calc[[j]][1,] = wtd.quantile(peak_time_calc[,2],out2$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  kd.peak.calc[[j]][1,] = wtd.quantile(peak_time_calc[,3],out3$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  bf.peak.calc[[j]][2,] = wtd.quantile(peak_intensity_calc[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  apf.peak.calc[[j]][2,] = wtd.quantile(peak_intensity_calc[,2],out2$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
  kd.peak.calc[[j]][2,] = wtd.quantile(peak_intensity_calc[,3],out3$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
}

save.image(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/sir.pf.pred",param,"-",n,".rdata",sep=""))

# Plot predicted % population infected
# Bootstrap filter
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/EpidPred-BF",param,"-",n,".pdf",sep=""),
	width=10,height=8)
par(mfrow=c(2,3))
ymax = max(sapply(bf.ui,max),max(sim$x[1,]))
for(j in 1:length(days))
{
  plot(sim$x[1,],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab="% Population",
	main="Predicted %Pop Infected")
  mtext(paste("k = ",days[j]," Days",sep=""),side=3)
  lines((days[j]+1):nt,bf.i[[j]],col=2)
  lines((days[j]+1):nt,bf.li[[j]],lty=2,col=2)
  lines((days[j]+1):nt,bf.ui[[j]],lty=2,col=2)
  if(j == 1){ legend("topright",legend=c("Truth","Median","95% CI"),
	lty=c(1,1,2),col=c(1,2,2))}
}
dev.off()

# Auxiliary particle filter
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/EpidPred-APF",param,"-",n,".pdf",sep=""),
	width=10,height=8)
par(mfrow=c(2,3))
ymax = max(sapply(apf.ui,max),max(sim$x[1,]))
for(j in 1:length(days))
{
  plot(sim$x[1,],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab="% Population",
	main="Predicted %Pop Infected")
  mtext(paste("k = ",days[j]," Days",sep=""),side=3)
  lines((days[j]+1):nt,apf.i[[j]],col=4)
  lines((days[j]+1):nt,apf.li[[j]],lty=2,col=4)
  lines((days[j]+1):nt,apf.ui[[j]],lty=2,col=4)
  if(j == 1){ legend("topright",legend=c("Truth","Median","95% CI"),
	lty=c(1,1,2),col=c(1,4,4))}
}
dev.off()

# Kernel density particle filter
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/EpidPred-KD",param,"-",n,".pdf",sep=""),
	width=10,height=8)
par(mfrow=c(2,3))
ymax = max(sapply(kd.ui,max),max(sim$x[1,]))
for(j in 1:length(days))
{
  plot(sim$x[1,],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab="% Population",
	main="Predicted %Pop Infected")
  mtext(paste("k = ",days[j]," Days",sep=""),side=3)
  lines((days[j]+1):nt,kd.i[[j]],col=3)
  lines((days[j]+1):nt,kd.li[[j]],lty=2,col=3)
  lines((days[j]+1):nt,kd.ui[[j]],lty=2,col=3)
  if(j == 1){ legend("topright",legend=c("Truth","Median","95% CI"),
	lty=c(1,1,2),col=c(1,3,3))}
}
dev.off()

# Plot predicted peak times and intensities
tpeak = which(sim$x[1,]==max(sim$x[1,]))-1
ipeak = max(sim$x[1,])

# Bootstrap filter
times = sapply(bf.peak.prop,function(x) x[1,2])
times.l = sapply(bf.peak.prop,function(x) x[1,1])
times.u = sapply(bf.peak.prop,function(x) x[1,3])
timesc = sapply(bf.peak.calc,function(x) x[1,2])
timesc.l = sapply(bf.peak.calc,function(x) x[1,1])
timesc.u = sapply(bf.peak.calc,function(x) x[1,3])
ints = sapply(bf.peak.prop,function(x) x[2,2])
ints.l = sapply(bf.peak.prop,function(x) x[2,1])
ints.u = sapply(bf.peak.prop,function(x) x[2,3])
intsc = sapply(bf.peak.calc,function(x) x[2,2])
intsc.l = sapply(bf.peak.calc,function(x) x[2,1])
intsc.u = sapply(bf.peak.calc,function(x) x[2,3])
ymin = min(times.l,timesc.l,tpeak)
ymax = max(times.u,timesc.u,tpeak)
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Peaks-BF",param,"-",n,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2))
plot(days,times,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),
	xlab="Days",ylab="Peak (Days)",main="Predicted Time of Epidemic Peak")
points(days+0.5,timesc,pch=2)
segments(days,times.l,days,times.u)
segments(days+.5,timesc.l,days+.5,timesc.u)
abline(h=tpeak,lty=2)
legend("top",legend=c("Prop","Calc","95% Cred Int","Truth"),lty=c(NA,NA,1,2),
	pch=c(1,2,NA,NA))
ymin = min(ints.l,intsc.l,ipeak)
ymax = max(ints.u,intsc.u,ipeak)
plot(days,ints,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),
	xlab="Days",ylab="Proportion Infected",
	main="Predicted Intensity of Epidemic Peak")
points(days+0.5,intsc,pch=2)
segments(days,ints.l,days,ints.u)
segments(days+.5,intsc.l,days+.5,intsc.u)
abline(h=ipeak,lty=2)
dev.off()

# Auxiliary particle filter
times = sapply(apf.peak.prop,function(x) x[1,2])
times.l = sapply(apf.peak.prop,function(x) x[1,1])
times.u = sapply(apf.peak.prop,function(x) x[1,3])
timesc = sapply(apf.peak.calc,function(x) x[1,2])
timesc.l = sapply(apf.peak.calc,function(x) x[1,1])
timesc.u = sapply(apf.peak.calc,function(x) x[1,3])
ints = sapply(apf.peak.prop,function(x) x[2,2])
ints.l = sapply(apf.peak.prop,function(x) x[2,1])
ints.u = sapply(apf.peak.prop,function(x) x[2,3])
intsc = sapply(apf.peak.calc,function(x) x[2,2])
intsc.l = sapply(apf.peak.calc,function(x) x[2,1])
intsc.u = sapply(apf.peak.calc,function(x) x[2,3])
ymin = min(times.l,timesc.l,tpeak)
ymax = max(times.u,timesc.u,tpeak)
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Peaks-APF",param,"-",n,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2))
plot(days,times,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),
	xlab="Days",ylab="Peak (Days)",main="Predicted Time of Epidemic Peak")
points(days+0.5,timesc,pch=2)
segments(days,times.l,days,times.u)
segments(days+.5,timesc.l,days+.5,timesc.u)
abline(h=tpeak,lty=2)
legend("top",legend=c("Prop","Calc","95% Cred Int","Truth"),lty=c(NA,NA,1,2),
	pch=c(1,2,NA,NA))
ymin = min(ints.l,intsc.l,ipeak)
ymax = max(ints.u,intsc.u,ipeak)
plot(days,ints,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),
	xlab="Days",ylab="Proportion Infected",
	main="Predicted Intensity of Epidemic Peak")
points(days+0.5,intsc,pch=2)
segments(days,ints.l,days,ints.u)
segments(days+.5,intsc.l,days+.5,intsc.u)
abline(h=ipeak,lty=2)
dev.off()

# Kernel density particle filter
times = sapply(kd.peak.prop,function(x) x[1,2])
times.l = sapply(kd.peak.prop,function(x) x[1,1])
times.u = sapply(kd.peak.prop,function(x) x[1,3])
timesc = sapply(kd.peak.calc,function(x) x[1,2])
timesc.l = sapply(kd.peak.calc,function(x) x[1,1])
timesc.u = sapply(kd.peak.calc,function(x) x[1,3])
ints = sapply(kd.peak.prop,function(x) x[2,2])
ints.l = sapply(kd.peak.prop,function(x) x[2,1])
ints.u = sapply(kd.peak.prop,function(x) x[2,3])
intsc = sapply(kd.peak.calc,function(x) x[2,2])
intsc.l = sapply(kd.peak.calc,function(x) x[2,1])
intsc.u = sapply(kd.peak.calc,function(x) x[2,3])
ymin = min(times.l,timesc.l,tpeak)
ymax = max(times.u,timesc.u,tpeak)
pdf(paste("C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/Peaks-KD",param,"-",n,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2))
plot(days,times,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),
	xlab="Days",ylab="Peak (Days)",main="Predicted Time of Epidemic Peak")
points(days+0.5,timesc,pch=2)
segments(days,times.l,days,times.u)
segments(days+.5,timesc.l,days+.5,timesc.u)
abline(h=tpeak,lty=2)
legend("top",legend=c("Prop","Calc","95% Cred Int","Truth"),lty=c(NA,NA,1,2),
	pch=c(1,2,NA,NA))
ymin = min(ints.l,intsc.l,ipeak)
ymax = max(ints.u,intsc.u,ipeak)
plot(days,ints,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),
	xlab="Days",ylab="Proportion Infected",
	main="Predicted Intensity of Epidemic Peak")
points(days+0.5,intsc,pch=2)
segments(days,ints.l,days,ints.u)
segments(days+.5,intsc.l,days+.5,intsc.u)
abline(h=ipeak,lty=2)
dev.off()