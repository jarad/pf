# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

## Plot contours of p-values for time v threshold

# Set graphical parameters
cex.lab = 5
cex.main = 5
cex.axis = 3
cex.sub = 2.5
cex.leg = 1.7

# 3 by 3 figure of contours of pvals by threshold by time (rows = nonunif, cols = param)
params = expression(beta,gamma,nu)
nonunifs = c("cov","entropy","ess")
threshs = seq(.05,.95,.05)
breaks = c(0,.01,.05,.1,.5,.9)
cols = c("red",heat.colors(length(breaks)-1)[3:(length(breaks)-1)],"white")
legend = levels(cut(seq(0,1,.005),breaks=breaks))
n = 10000
load(paste(dpath,"PF-KD-normal-stratified-",n,"-ks.rdata",sep=""))

# Construct plot
pdf(paste(gpath,"PF-KD-normal-stratified-",n,"-ks.pdf",sep=""),width=20,height=20)
par(mfrow=c(3,3),mar=c(9,13,7,1)+.1,mgp=c(7,1.5,0))
for(i in 1:length(nonunifs))
{
  for(j in 1:length(params))
  {
      # Plot images
      x1 = 1:dim(mypvals)[3]
      x2 = threshs
      z = t(mypvals[i,,,j])
      if(i == 1 & j == 1)
      { 
        image(x1,x2,z,col=cols,breaks=breaks,xlab="Time (days)",ylab=nonunifs[i],main=params[j],cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis)
        mtext("Threshold",side=2,line=3.75,cex=cex.sub)
        legend("topright",legend=legend,fill=cols,cex=cex.leg)
      } else if(i == 1) {
        image(x1,x2,z,col=cols,breaks=breaks,xlab="",ylab="",main=params[j],cex.main=cex.main,cex.axis=cex.axis)
      } else if(j == 1) {
        image(x1,x2,z,col=cols,breaks=breaks,xlab="",ylab=nonunifs[i],cex.lab=cex.lab,cex.axis=cex.axis)
      } else {
        image(x1,x2,z,col=cols,breaks=breaks,xlab="",ylab="",cex.axis=cex.axis)
      }

      # Calculate maximum p-value at each time point
      max_ind = apply(mypvals[i,,,j],2,function(x) which(x == max(x)))
      max_ind = as.list(max_ind)

      # Plot points where maximum p-value reached
      for(k in 1:dim(mypvals)[3])
      {
	  x = rep(k,length(max_ind[[k]]))
        y = threshs[max_ind[[k]]]
        points(x,y,pch=3)
      }
  }
}
dev.off()

## Plot p-values versus threshold for selected time points

# Set graphical parameters
cex.lab = 4
cex.main = 4
cex.axis = 2
cex.leg = 1.5

# 2 by 4 figure of panels of pvals by threshold, each panel a diff time point
params = expression(beta,gamma,nu)
nonunifs = c("cov","entropy","ess")
threshs = seq(.05,.95,.05)
cutoff = seq(from=15,by=15,length=8)
cols = c(2,4,3)
n = 10000
load(paste(dpath,"PF-KD-normal-stratified-",n,"-ks.rdata",sep=""))

# Construct plot
for(k in 1:length(params))
{
  pdf(paste(gpath,"PF-KD-normal-stratified-",n,"-kstime-",params[k],".pdf",sep=""),width=20,height=10)
  par(mfrow=c(2,4),mar=c(9,13,7,1)+.1,mgp=c(7,1.5,0))
  for(i in 1:length(cutoff))
  {
    for(j in 1:length(nonunifs))
    {
      # Plot images
      x1 = threshs
      x2 = mypvals[j,,cutoff[i],k]
      if(j == 1)
      {
        if(i == 1)
        {   
          plot(x1,x2,type="l",ylim=c(0,1),col=cols[j],xlab="Threshold",ylab="P-value",main=paste("t = ",cutoff[i],sep=""),cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis)
          legend("topright",legend=nonunifs,lty=rep(1,length(nonunifs)),col=cols,cex=cex.leg)
        } else {
          plot(x1,x2,type="l",ylim=c(0,1),col=cols[j],xlab="",ylab="",main=paste("t = ",cutoff[i],sep=""),cex.main=cex.main,cex.axis=cex.axis)
        }
      } else {
        lines(x1,x2,col=cols[j])
      }
    }
  }
  dev.off()
}