source("pf_functions.r")

# Set graphics and data path
gpath = "../graphs/"
dpath = "/storage/sheinson_research/"

## Figure 1 - Compare particle filters over different # particles for systematic resampling, uniform priors

# Set graphical parameters
params = expression(beta,gamma,nu)
ymins = c(0.15,0.085,0.9)
ymaxs = c(0.35,0.145,1.35)
cols = c(2,4,3)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4

# length(ns) by length(params) figure of plot panels (rows = num particles, cols = params) with length(filts.load) credible bounds
filts.load = c("KD","RM")
ns = c(100)

# Construct plots
pdf(paste(gpath,"PF-systematic-uniform-",ns[length(ns)],".pdf",sep=""),width=10*length(params),height=10*length(ns))
par(mfrow=c(length(ns),length(params)),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:length(ns))
{
  for(k in 1:length(params))
  {
    for(j in 1:length(filts.load))
    {
      load(paste(dpath,"PF-quant-orig-orig-",filts.load[j],"-uniform-systematic-",ns[i],"-ess-80.rdata",sep=""))
      out = pf.quant.out$theta.quant
      tt = dim(out)[1]; nt = tt - 1
      if(j == 1) # call plot function
      {
        if(k == 1 & i == 1) # label y axis and title
        {
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = 
          ",ns[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(k == 1 & i == length(ns)) { # label x and y axes
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab=paste("J = 
          ",ns[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(k == 1) { # label y axis only
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = 
          ",ns[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(i == 1) { # label title only
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(i == length(ns)) { # label x axis only
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time 
          (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else { # label nothing
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        }
      } else { # lines only
        lines(1:nt,out[-1,k,2],col=cols[j])
        lines(1:nt,out[-1,k,3],col=cols[j])
      }
    }
    load(paste(dpath,"sim-orig.rdata",sep=""))
    abline(h=mysim$true.params$theta[k],col="gray47",lwd=6)
    if(k == 1 & i == 1) # add legend
    {
      legend("topright",legend=c("Truth",filts.load),col=c("gray47",cols),lty=c(1,rep(1,length(filts.load))),lwd=c(6,rep(1,length(filts.load))),cex=cex.leg)
    }
  }
}
dev.off()
