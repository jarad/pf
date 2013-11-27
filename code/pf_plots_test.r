# Set graphics and data path
gpath = "../graphs/"
dpath = "../data/"

## Figure 1 - Compare particle filters over different # particles for systematic resampling, uniform priors

# Set graphical parameters
params = expression(beta,gamma,nu)
ymins = c(0.14,0.09,0.95)
ymaxs = c(0.50,0.143,1.3)
cols = c(2,4,3)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4
pic.fac = 10
burn = 0

# length(n) by length(params) figure of plot panels (rows = num particles, cols = params) with length(filts.load) credible bounds
filt = c("BF","APF","KD")
n = c(100, 1000, 10000, 20000, 40000)
n.sim = 1

# Construct plots
pdf(paste(gpath,"PF-uniform-systematic-filt-",n.sim,".pdf",sep=""),width=pic.fac*length(params),height=pic.fac*length(n))
par(mfrow=c(length(n),length(params)),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:length(n))
{
  for(k in 1:length(params))
  {
    for(j in 1:length(filt))
    {
      load(paste(dpath,"PF-quant-uniform-systematic-",filt[j],"-",n[i],"-logit-61-",n.sim,".rdata",sep=""))
      out = pf.quant.out$theta.quant
      tt = dim(out)[1]; nt = tt - 1
      if(burn > 0)
      {
        quant = out[-burn,k,]
        x = (0:nt)[-burn]
      } else {
        quant = out[,k,]
        x = 0:nt
      }
      if(j == 1) # call plot function
      {
        if(k == 1 & i == 1) # label y axis and title
        {
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(k == 1 & i == length(n)) { # label x and y axes
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(k == 1) { # label y axis only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(i == 1) { # label title only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(i == length(n)) { # label x axis only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else { # label nothing
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
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
      legend("topright",legend=c("Truth",filt),col=c("gray47",cols),lty=c(1,rep(1,length(filt))),lwd=c(6,rep(1,length(filt))),cex=cex.leg)
    }
  }
}
dev.off()

## Figure 2 - Compare particle filters over different # particles for different resampling methods (KD pf, lognormal priors, delta = .99)

# Set graphical parameters
params = expression(beta,gamma,nu)
ymins = c(0.14,0.09,0.95)
ymaxs = c(0.50,0.15,1.4)
cols = c(2,4,3,6)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4
pic.fac = 10
burn = 0

# length(n) by length(params) figure of plot panels (rows = num particles, cols = params) with length(filts.load) credible bounds
resamp = c("multinomial","residual","stratified","systematic")
n = c(100, 1000, 10000, 20000, 40000)
n.sim = 1

# Construct plots
pdf(paste(gpath,"PF-KD-lognormal-resamp-",n.sim,".pdf",sep=""),width=pic.fac*length(params),height=pic.fac*length(n))
par(mfrow=c(length(n),length(params)),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:length(n))
{
  for(k in 1:length(params))
  {
    for(j in 1:length(resamp))
    {
      load(paste(dpath,"PF-quant-KD-lognormal-",resamp[j],"-",n[i],"-orig-0.99-61-",n.sim,".rdata",sep=""))
      out = pf.quant.out$theta.quant
      tt = dim(out)[1]; nt = tt - 1
      if(burn > 0)
      {
        quant = out[-burn,k,]
        x = (0:nt)[-burn]
      } else {
        quant = out[,k,]
        x = 0:nt
      }
      if(j == 1) # call plot function
      {
        if(k == 1 & i == 1) # label y axis and title
        {
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(k == 1 & i == length(n)) { # label x and y axes
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(k == 1) { # label y axis only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(i == 1) { # label title only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(i == length(n)) { # label x axis only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else { # label nothing
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
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
      legend("topright",legend=c("Truth",resamp),col=c("gray47",cols),lty=c(1,rep(1,length(resamp))),lwd=c(6,rep(1,length(resamp))),cex=cex.leg)
    }
  }
}
dev.off()

## Figure 3 - Compare particle filters over different # particles for different delta (KD pf, stratified resampling, lognormal priors)

# Set graphical parameters
params = expression(beta,gamma,nu)
ymins = c(0.14,0.09,0.95)
ymaxs = c(0.50,0.143,1.45)
cols = rainbow(6)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4
pic.fac = 10
burn = 0

# length(n) by length(params) figure of plot panels (rows = num particles, cols = params) with length(filts.load) credible bounds
delta = c(0.9, 0.95, 0.96, 0.97, 0.98, 0.99)
n = c(100, 1000, 10000, 20000, 40000)
n.sim = 1

# Construct plots
pdf(paste(gpath,"PF-KD-lognormal-stratified-delta-",n.sim,".pdf",sep=""),width=pic.fac*length(params),height=pic.fac*length(n))
par(mfrow=c(length(n),length(params)),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:length(n))
{
  for(k in 1:length(params))
  {
    for(j in 1:length(delta))
    {
      load(paste(dpath,"PF-quant-KD-lognormal-stratified-",n[i],"-orig-",delta[j],"-61-",n.sim,".rdata",sep=""))
      out = pf.quant.out$theta.quant
      tt = dim(out)[1]; nt = tt - 1
      if(burn > 0)
      {
        quant = out[-burn,k,]
        x = (0:nt)[-burn]
      } else {
        quant = out[,k,]
        x = 0:nt
      }
      if(j == 1) # call plot function
      {
        if(k == 1 & i == 1) # label y axis and title
        {
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(k == 1 & i == length(n)) { # label x and y axes
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(k == 1) { # label y axis only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(i == 1) { # label title only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else if(i == length(n)) { # label x axis only
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
        } else { # label nothing
          plot(x,quant[,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",cex.axis=cex.axis)
          lines(x,quant[,3],col=cols[j])
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
      legend("topright",legend=c("Truth",delta),col=c("gray47",cols),lty=c(1,rep(1,length(delta))),lwd=c(6,rep(1,length(delta))),cex=cex.leg)
    }
  }
}
dev.off()