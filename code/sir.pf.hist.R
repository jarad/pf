# Set data path
dpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Data/"

# Load data
load(paste(dpath,"sir.pf.test-1000.rdata",sep=""))

# How many parameters?
p = 3

# Set graphics path
gpath = "C:/Users/Danny/Dropbox/SIR_Particle_Filtering/Graphs/PF/"

# Plot histograms of unknown parameters over time
require(smcUtils)
cutoff = seq(15,120,15)
msize = labsize = 1.5
expr = expression(beta,gamma,nu,b,varsigma,sigma)
parnames = c("beta","gamma","nu","b","varsigma","sigma")
for(j in 1:p)
{
  x = theta2u(out3$theta[j,,cutoff],thetal[j],thetau[j])
  w = out3$weight[,cutoff]
  for(i in 1:length(cutoff))
  {
    tmp = resample(w[,i], method="stratified", nonuniformity="none")
    x[,i] = x[tmp$indices,i]
  }
  xmin = min(apply(x,2,function(x) min(hist(x,plot=FALSE)$breaks)))
  xmax = max(apply(x,2,function(x) max(hist(x,plot=FALSE)$breaks)))
  ymin = min(apply(x,2,function(x) min(hist(x,plot=FALSE)$density)))
  ymax = max(apply(x,2,function(x) max(hist(x,plot=FALSE)$density)))
  pdf(paste(gpath,"Hist-Time-KD-",p,"P-",n,"-",parnames[j],".pdf",sep=""),width=10,height=5)
  par(mfrow=c(2,4))
  for(i in 1:length(cutoff))
  {
    if(i==1){
      hist(x[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),freq=FALSE,
        main=paste("k = ",cutoff[i],sep=""),xlab=expr[j],
	  cex.main=msize,cex.lab=labsize)
      abline(v=theta[j])
      if(j == 1) mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=theta[1])),side=3,cex=.85)
      if(j == 2) mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=theta[2])),side=3,cex=.85)
      if(j == 3) mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=theta[3])),side=3,cex=.85)
      if(j == 4) mtext(substitute(paste(b," = ",aa,sep=""),list(aa=theta[4])),side=3,cex=.85)
      if(j == 5) mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=theta[5])),side=3,cex=.85)
      if(j == 6) mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=theta[6])),side=3,cex=.85)
    } else {
      hist(x[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),freq=FALSE,
        main=paste("k = ",cutoff[i],sep=""),xlab="",ylab="",
	  cex.main=msize,cex.lab=labsize)
      abline(v=theta[j])
    }
  }
  dev.off()
}

# Plot 2D histograms of beta and gamma over time
require(smcUtils)
require(KernSmooth)
cutoff = seq(15,120,15)
msize = labsize = 1.5
xb = theta2u(out3$theta[1,,cutoff],thetal[1],thetau[1])
xg = theta2u(out3$theta[2,,cutoff],thetal[2],thetau[2])
w = out3$weight[,cutoff]
xmin = Inf; xmax = -Inf; ymin = Inf; ymax = -Inf
h = 0.05
for(i in 1:length(cutoff))
{
  tmp = resample(w[,i], method="stratified", nonuniformity="none")
  xb[,i] = xb[tmp$indices,i]
  xg[,i] = xg[tmp$indices,i]
  a = bkde2D(cbind(xb[,i],xg[,i]),h)
  xmin = min(xmin,a$x1)
  xmax = max(xmax,a$x1)
  ymin = min(ymin,a$x2)
  ymax = max(ymax,a$x2)
}

pdf(paste(gpath,"Hist-Time-KD-",p,"P-",n,"-betagamma.pdf",sep=""),width=10,height=5)
par(mfrow=c(2,4))
for(i in 1:length(cutoff))
{
  a = bkde2D(cbind(xb[,i],xg[,i]),h,range.x=list(x1=c(xmin,xmax),x2=c(ymin,ymax)))
  if(i==1){
    image(a$x1,a$x2,a$fhat,xlim=c(xmin,xmax),ylim=c(ymin,ymax),
        main=paste("k = ",cutoff[i],sep=""),xlab=expression(beta),
	  ylab=expression(gamma),cex.main=msize,cex.lab=labsize)
    contour(a$x1,a$x2,a$fhat,add=TRUE)
    points(theta[1],theta[2],pch=20,col=4)
    mtext(substitute(paste(beta," = ",aa,", ",gamma," = ",ab,sep=""),list(aa=theta[1],ab=theta[2])),
	side=3,cex=.85)
  } else {
    image(a$x1,a$x2,a$fhat,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",
        main=paste("k = ",cutoff[i],sep=""),cex.main=msize,cex.lab=labsize)
    contour(a$x1,a$x2,a$fhat,add=TRUE)
    points(theta[1],theta[2],pch=20,col=4)
  }
}
dev.off()
