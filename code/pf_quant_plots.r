# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

## Plot quantiles for different nonuniformities/thresholds/params over time

# Set graphical parameters
params = expression(beta,gamma,nu)
ymins = c(0.20,0.075,0.8)
ymaxs = c(0.30,0.165,1.45)
cols = c(2,4,3)
cex.lab = 5
cex.main = 5
cex.axis = 3
cex.leg = 2.5

# 3 by 3 figure of plot panels for each param (panel = thresh, lines = nonunif)
nonunifs = c("ess","cov","entropy")
threshs = seq(.1,.9,.1)
n = 20000

# Construct plots
for(k in 1:length(params))
{
  pdf(paste(gpath,"PF-KD-normal-stratified-",n,"-",params[k],".pdf",sep=""),width=20,height=20)
  par(mfrow=c(3,3),mar=c(7,10,7,1)+.1,mgp=c(5,1.5,0))
  for(i in 1:length(threshs))
  {
    for(j in 1:length(nonunifs))
    {
      load(paste(dpath,"PF-quant-KD-normal-stratified-",n,"-",nonunifs[j],"-",100*threshs[i],".rdata",sep=""))
      out = pf.quant.out$theta.quant
      tt = dim(out)[1]; nt = tt - 1
      if(j == 1) # call plot function
      {
        if(i == 1) # label x, y axes, title, and legend
        {
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab=params[k],main=paste("Threshold = ",threshs[i],sep=""),cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
          legend("topright",legend=c("Truth",nonunifs),col=c("gray47",cols),lty=c(1,1,1,1),cex=cex.leg)
        } else { # title, axis ticks only
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",main=paste("Threshold = ",threshs[i],sep=""),cex.axis=cex.axis,cex.main=cex.main)
          lines(1:nt,out[-1,k,3],col=cols[j])
        }
      } else { # lines only
        lines(1:nt,out[-1,k,2],col=cols[j])
        lines(1:nt,out[-1,k,3],col=cols[j])
      }
      load(paste(dpath,"sim-xy.rdata",sep=""))
      abline(h=theta[k],col="gray47")
    }  
  }
  dev.off()
}
