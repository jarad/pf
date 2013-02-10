# Set graphics and data path
gpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/"
dpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/"

# How many particles? Set n = 100, 1000, or 10000
n = 10000
# How many parameters? Set p = 3 or p = 6
p = 3
if(p == 3) param = "" else param = "-6P"

# Load data
load(paste(dpath,"sir.pf",param,"-",n,".rdata",sep=""))

# Resample particles at pre-specified times so they have equal weights
require(smcUtils)
M = 10000
cutoff = seq(16,121,15)
w = out3$weight[,cutoff]
tmps = matrix(NA,nr=M,nc=length(cutoff))
for(i in 1:length(cutoff)) tmps[,i] = resample(w[,i], M, method="stratified", nonuniformity="none")$indices

# Plot histograms of unknown parameters over time
msize = labsize = 1.5
expr = expression(beta,gamma,nu,b,varsigma,sigma)
parnames = c("beta","gamma","nu","b","varsigma","sigma")
st = 3:8
for(j in 1:p)
{
  x = list(bf=out$state[st[j],,cutoff],apf=out2$state[st[j],,cutoff],kd=theta2u(out3$theta[j,,cutoff],thetal[j],thetau[j]))
  pf = c("BF","APF","KD")
  for(k in 1:length(pf))
  {
    xrw = matrix(NA,nr=M,nc=length(cutoff))
    for(i in 1:length(cutoff)) xrw[,i] = x[[k]][tmps[,i],i]
    xmin = min(apply(xrw,2,function(x) min(hist(x,plot=FALSE)$breaks)))
    xmax = max(apply(xrw,2,function(x) max(hist(x,plot=FALSE)$breaks)))
    ymin = min(apply(xrw,2,function(x) min(hist(x,plot=FALSE)$density)))
    ymax = max(apply(xrw,2,function(x) max(hist(x,plot=FALSE)$density)))
    pdf(paste(gpath,"Hist-",pf[k],param,"-",n,"-",parnames[j],".pdf",sep=""),width=10,height=5)
    par(mfrow=c(2,4)) # dimensions should depend on cutoff
    for(i in 1:length(cutoff))
    {
      if(i==1){
        hist(xrw[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),freq=FALSE,
          main=paste("k = ",cutoff[i]-1,sep=""),xlab=expr[j],
	    cex.main=msize,cex.lab=labsize)
        abline(v=theta[j],col=2,lwd=2)
        if(j == 1) mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=theta[1])),side=3,cex=.85)
        if(j == 2) mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=theta[2])),side=3,cex=.85)
        if(j == 3) mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=theta[3])),side=3,cex=.85)
        if(j == 4) mtext(substitute(paste(b," = ",aa,sep=""),list(aa=theta[4])),side=3,cex=.85)
        if(j == 5) mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=theta[5])),side=3,cex=.85)
        if(j == 6) mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=theta[6])),side=3,cex=.85)
      } else {
        hist(xrw[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),freq=FALSE,
          main=paste("k = ",cutoff[i]-1,sep=""),xlab="",ylab="",
	    cex.main=msize,cex.lab=labsize)
        abline(v=theta[j],col=2,lwd=2)
      }
    }
    dev.off()
  }
}

# Scatterplot of beta v gamma over time
xb = list(bf=out$state[3,,cutoff],apf=out2$state[3,,cutoff],kd=theta2u(out3$theta[1,,cutoff],thetal[1],thetau[1]))
xg = list(bf=out$state[4,,cutoff],apf=out2$state[4,,cutoff],theta2u(out3$theta[2,,cutoff],thetal[2],thetau[2]))
for(k in 1:length(pf))
{
  xrwb = matrix(NA,nr=M,nc=length(cutoff))
  xrwg = matrix(NA,nr=M,nc=length(cutoff))
  xmin = Inf; xmax = -Inf; ymin = Inf; ymax = -Inf
  for(i in 1:length(cutoff))
  {
    xrwb[,i] = xb[[k]][tmps[,i],i]
    xrwg[,i] = xg[[k]][tmps[,i],i]
    xmin = min(xmin,xrwb[,i])
    xmax = max(xmax,xrwb[,i])
    ymin = min(ymin,xrwg[,i])
    ymax = max(ymax,xrwg[,i])
  }
  pdf(paste(gpath,"Hist-",pf[k],param,"-",n,"-betagamma.pdf",sep=""),width=10,height=5)
  par(mfrow=c(2,4))
  for(i in 1:length(cutoff))
  {
    if(i == 1)
    {
      plot(xrwb[,i],xrwg[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=expression(beta),ylab=expression(gamma),main=paste("k = ",cutoff[i]-1,sep=""),cex.main=msize,cex.lab=labsize)
      mtext(substitute(paste(beta," = ",aa,", ",gamma," = ",ab,sep=""),list(aa=theta[1],ab=theta[2])),side=3,cex=.85)
	points(theta[1],theta[2],col=4,pch=20)
    } else {
      plot(xrwb[,i],xrwg[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",main=paste("k = ",cutoff[i]-1,sep=""),cex.main=msize)
	points(theta[1],theta[2],col=4,pch=20)
    }
  }
  dev.off()
}