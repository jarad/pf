# Set graphics and data path
gpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF-D2/"
dpath = "C:/Users/Danny/My Documents/UCSB - Research/pf/data/D2/"

# Load PF data
load(paste(dpath,"sir.pf-2P-unifpriors-1000.rdata",sep=""))

# Which parameters unknown, labels
p = 1:2
expr = expression(beta,gamma,nu)
parnames = c("beta","gamma","nu")
all.theta = c(theta,b,varsigma,sigma)[p]

# Create function to map theta from real line to original scale
ftheta = function(theta,thetal,thetau)
{
  theta2u(theta,thetal,thetau)
}

# Resample particles at pre-specified times so they have equal weights
require(smcUtils)
M = 10000
cutoff = seq(16,121,15)
w = array(NA,dim=c(3,n,length(cutoff)))
w[1,,] = out$weight[,cutoff]
w[2,,] = out2$weight[,cutoff]
w[3,,] = out3$weight[,cutoff]
tmps = array(NA,dim=c(3,M,length(cutoff)))
for(i in 1:length(cutoff))
{
  tmps[1,,i] = resample(w[1,,i], M, method="stratified", nonuniformity="none")$indices
  tmps[2,,i] = resample(w[2,,i], M, method="stratified", nonuniformity="none")$indices
  tmps[3,,i] = resample(w[3,,i], M, method="stratified", nonuniformity="none")$indices
}

# Plot histograms of unknown parameters over time
msize = labsize = 1.5
st = p+2
for(j in p)
{
  x = list(bf=ftheta(out$state[st[j],,cutoff],thetal[j],thetau[j]),apf=ftheta(out2$state[st[j],,cutoff],thetal[j],thetau[j]),kd=ftheta(out3$theta[j,,cutoff],thetal[j],thetau[j]))
  pf = c("BF","APF","KD")
  for(k in 1:length(pf))
  {
    xrw = matrix(NA,nr=M,nc=length(cutoff))
    for(i in 1:length(cutoff)) xrw[,i] = x[[k]][tmps[k,,i],i]
    xmin = min(apply(xrw,2,function(x) min(hist(x,plot=FALSE)$breaks)))
    xmax = max(apply(xrw,2,function(x) max(hist(x,plot=FALSE)$breaks)))
    ymin = min(apply(xrw,2,function(x) min(hist(x,plot=FALSE)$density)))
    ymax = max(apply(xrw,2,function(x) max(hist(x,plot=FALSE)$density)))
    pdf(paste(gpath,"Hist-",pf[k],param,"-",n,"-",parnames[j],".pdf",sep=""),width=10,height=5)
    par(mfrow=c(2,4),mar=c(4,5.2,3,.5)+.1) # dimensions should depend on cutoff
    for(i in 1:length(cutoff))
    {
      if(i==1){
        hist(xrw[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),freq=FALSE,main=paste("t = ",cutoff[i]-1,sep=""),xlab=expr[j],cex.main=msize,cex.lab=labsize)
        abline(v=all.theta[j],col=2,lwd=2)
        mtext(paste("Truth = ",all.theta[j],sep=""),side=3,cex=.65)
      } else {
        hist(xrw[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),freq=FALSE,main=paste("t = ",cutoff[i]-1,sep=""),xlab="",ylab="",cex.main=msize)
        abline(v=all.theta[j],col=2,lwd=2)
      }
    }
    dev.off()
  }
}

# Scatterplot of 2 unknown parameters over time
p1 = 1; p2 = 2
xb = list(bf=ftheta(out$state[p1+2,,cutoff],thetal[p1],thetau[p1]),apf=ftheta(out2$state[p1+2,,cutoff],thetal[p1],thetau[p1]),kd=ftheta(out3$theta[p1,,cutoff],thetal[p1],thetau[p1]))
xg = list(bf=ftheta(out$state[p2+2,,cutoff],thetal[p2],thetau[p2]),apf=ftheta(out2$state[p2+2,,cutoff],thetal[p2],thetau[p2]),ftheta(out3$theta[p2,,cutoff],thetal[p2],thetau[p2]))
for(k in 1:length(pf))
{
  xrwb = matrix(NA,nr=M,nc=length(cutoff))
  xrwg = matrix(NA,nr=M,nc=length(cutoff))
  xmin = Inf; xmax = -Inf; ymin = Inf; ymax = -Inf
  for(i in 1:length(cutoff))
  {
    xrwb[,i] = xb[[k]][tmps[k,,i],i]
    xrwg[,i] = xg[[k]][tmps[k,,i],i]
    xmin = min(xmin,xrwb[,i])
    xmax = max(xmax,xrwb[,i])
    ymin = min(ymin,xrwg[,i])
    ymax = max(ymax,xrwg[,i])
  }
  pdf(paste(gpath,"Hist-",pf[k],param,"-",n,"-",parnames[p1],parnames[p2],".pdf",sep=""),width=10,height=5)
  par(mfrow=c(2,4),mar=c(4,6,2,.5)+.1)
  for(i in 1:length(cutoff))
  {
    if(i == 1)
    {
      plot(xrwb[,i],xrwg[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=expr[p1],ylab=expr[p2],main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize,cex.lab=labsize)
      mtext(paste("Truth = ",all.theta[p1],sep=""),side=1,cex=.65)
	mtext(paste("Truth = ",all.theta[p2],sep=""),side=2,cex=.65)
      points(all.theta[p1],all.theta[p2],col=2,pch=20)
    } else {
      plot(xrwb[,i],xrwg[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize)
      points(all.theta[p1],all.theta[p2],col=2,pch=20)
    }
  }
  dev.off()
}