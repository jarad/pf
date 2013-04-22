source("pf_functions.r")

# Set graphics and data path
gpath = "../graphs/"
dpath = "../data/"

## Figure 1 - Compare particle filters over different # particles for systematic resampling, uniform priors

# Set graphical parameters
params = expression(beta,gamma,nu)
ymins = c(0.15,0.075,0.8)
ymaxs = c(0.35,0.165,1.45)
cols = c(2,4,3)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4

# 4 by 3 figure of plot panels (rows = num particles, cols = params)
filts = c("BF","APF","KD")
ns = c(100,1000,10000,20000)

# Construct plots
pdf(paste(gpath,"PF-systematic-uniform.pdf",sep=""),width=30,height=40)
par(mfrow=c(4,3),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:length(ns))
{
  for(k in 1:length(params))
  {
    for(j in 1:length(filts))
    {
      load(paste(dpath,"PF-quant-",filts[j],"-uniform-systematic-",ns[i],".rdata",sep=""))
      out = pf.quant.out$theta.quant
      tt = dim(out)[1]; nt = tt - 1
      if(j == 1) # call plot function
      {
        if(k == 1 & i == 1) # label y axis and title
        {
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",ns[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(k == 1 & i == length(ns)) { # label x and y axes
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab=paste("J = ",ns[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(k == 1) { # label y axis only
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",ns[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(i == 1) { # label title only
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
          lines(1:nt,out[-1,k,3],col=cols[j])
        } else if(i == length(ns)) { # label x axis only
          plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
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
    load(paste(dpath,"sim-xy.rdata",sep=""))
    abline(h=theta[k],col="gray47")
    if(k == 1 & i == 1) # add legend
    {
      legend("topright",legend=c("Truth","BF","APF","KD"),col=c("gray47",cols),lty=c(1,1,1,1),cex=cex.leg)
    }
  }
}
dev.off()

## Figure 2 - Create scatterplots of beta v gamma over time

# Set graphical parameters
xlim=c(.17,.33); ylim=c(.07,.17)
borderx=c(.14,.5); bordery=c(.09,.143)
msize=5; labsize=5; axsize=3; ptsize=3
ptsty=20; ptcol="gray70"; rline = -2.3; rsize = 1.8

# Load simulated data to get true values of beta and gamma
load(paste(dpath,"sim-xy.rdata",sep=""))

# Panel 1 - uniform prior draws, logit transformation
# Panel 2 - uniform prior draws, no transformation
# Panel 3 - normal prior draws, log transformation
priors = c("uniform","semi-uniform","normal")
for(k in 1:length(priors))
{
  # Load particle filtered data
  load(paste(dpath,"PF-KD-",priors[k],"-systematic-10000.rdata",sep=""))

  # Resample particles at cutoff points to have equal weights
  cutoff = seq(16, 61, len=4)
  betas = pf.out$ftheta(pf.out$out$theta[1,,],1)
  gammas = pf.out$ftheta(pf.out$out$theta[2,,],2)
  myout = array(NA,dim=c(2,dim(betas)[1],dim(betas)[2]))
  myout[1,,] = betas
  myout[2,,] = gammas
  myscat = pf.scat(myout,pf.out$out$weight,cutoff)

  # Scatterplots over time
  file = paste(gpath,"Hist-KD-",priors[k],"-systematic-10000-betagamma.pdf",sep="")
  pdf(file,width=20,height=5)
  par(mfrow=c(1,4),mar=c(7,10,5,1)+.1,mgp=c(6,1.55,0))
  for(i in 1:length(cutoff))
  {
    if(priors[k] == "uniform")
    {
      if(i == 1)
      {
        plot(myscat$xrw1[,i],myscat$xrw2[,i],col="white",xlim=xlim,ylim=ylim,xlab=expression(beta),ylab=expression(gamma),main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize,cex.lab=labsize,cex.axis=axsize)
      } else {
        plot(myscat$xrw1[,i],myscat$xrw2[,i],col="white",xlim=xlim,ylim=ylim,xlab="",ylab="",main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize,axes=FALSE)
        box()
      }
      abline(v=borderx,lty=2)
      abline(h=bordery,lty=2)
      points(myscat$xrw1[,i],myscat$xrw2[,i],col=ptcol,pch=ptsty,cex=ptsize)
      mtext(line=rline,cex=rsize,paste("r = ",round(myscat$r[i],2),sep=""))
      points(theta[1],theta[2],col=2,pch=3,cex=1.5*ptsize)
    } else {
      plot(myscat$xrw1[,i],myscat$xrw2[,i],col="white",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
      box()
      abline(v=borderx,lty=2)
      abline(h=bordery,lty=2)
      points(myscat$xrw1[,i],myscat$xrw2[,i],col=ptcol,pch=ptsty,cex=ptsize)
      mtext(line=rline,cex=rsize,paste("r = ",round(myscat$r[i],2),sep=""))
      points(theta[1],theta[2],col=2,pch=3,cex=1.5*ptsize)
    }
  }
  dev.off()
}

## Figure 2a - plot correlations over time for uniform, semi-uniform, normal priors

# Load data and compute correlations
load(paste(dpath,"sim-xy.rdata",sep=""))
tt = nt + 1
priors = c("uniform","semi-uniform","normal")
cutoff = 2:tt
r = matrix(NA,nr=length(cutoff),length(priors))
for(k in 1:length(priors))
{
  # Load particle filtered data
  load(paste(dpath,"PF-KD-",priors[k],"-systematic-10000.rdata",sep=""))

  # Resample particles at cutoff points to have equal weights
  betas = pf.out$ftheta(pf.out$out$theta[1,,],1)
  gammas = pf.out$ftheta(pf.out$out$theta[2,,],2)
  myout = array(NA,dim=c(2,dim(betas)[1],dim(betas)[2]))
  myout[1,,] = betas
  myout[2,,] = gammas
  myscat = pf.scat(myout,pf.out$out$weight,cutoff)
 
  # Save correlations
  r[,k] = myscat$r
}

# Set graphical parameters
msize=1; labsize=1; axsize=1
col = c(1,2,4)

# Construct plot
file = paste(gpath,"Hist-KD-systematic-10000-betagamma-corrs.pdf",sep="")
pdf(file)
#par(mar=c(7,10,5,1)+.1,mgp=c(6,1.55,0))
for(k in 1:dim(r)[2])
{
  if(k == 1)
  {
    plot(cutoff-1,r[,k],type="l",main="Correlation Coefficient vs Time",xlab="Time (days)",ylab=expression(r),cex.main=msize,cex.lab=labsize,cex.axis=axsize)
  } else {
    lines(cutoff-1,r[,k],col=col[k])
  }
  legend("topright",c("Unif draws, logit transf.","Unif draws, no transf.","Normal draws, log transf."),lty=c(1,1,1),col=c(1,2,4))
}
dev.off()

## Figure 3 - Compare resampling schemes over different # particles using normal priors, kernel density PF

# Set graphical parameters
ymins = c(0.15,0.075,0.8)
ymaxs = c(0.35,0.165,1.45)
cols = c(2,3,4,5)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4

# 4 by 3 figure of panels (rows = num particles, cols = params)
require(graphics)
resamps = c("multinomial","residual","stratified","systematic")
ns = c(100,1000,10000,20000)
params = expression(beta,gamma,nu)

# Construct plot
pdf(paste(gpath,"PF-KD-normal.pdf",sep=""),width=30,height=40)
par(mfrow=c(4,3),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:length(ns))
{
  for(k in 1:length(params))
  {
    # Calculate average quantiles of 4 resampling methods at most particles
    load(paste(dpath,"sim-xy.rdata",sep=""))
    nt = dim(sim$y)[2]
    out.avg = matrix(0,nr=nt,nc=2)
    for(j in 1:length(resamps))
    {
      load(paste(dpath,"PF-quant-KD-normal-",resamps[j],"-",ns[length(ns)],".rdata",sep=""))
      out.avg[,1] = out.avg[,1] + pf.quant.out$theta.quant[-1,k,2]
      out.avg[,2] = out.avg[,2] + pf.quant.out$theta.quant[-1,k,3]
    }
    out.avg = out.avg / length(resamps)
    x = 1:nt
    for(j in 0:length(resamps))
    {
      if(j == 0) # call plot function
      {
        if(k == 1 & i == 1) # label y axis and title
        {
          plot(1:nt,out.avg[,1],type="l",ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab=paste("J = ",ns[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          lines(1:nt,out.avg[,2],col="gray")
        } else if(k == 1 & i == length(ns)) { # label x and y axes
          plot(1:nt,out.avg[,1],type="l",ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="Time (days)",ylab=paste("J = ",ns[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out.avg[,2],col="gray")
        } else if(k == 1) { # label y axis only
          plot(1:nt,out.avg[,1],type="l",ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab=paste("J = ",ns[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out.avg[,2],col="gray")
        } else if(i == 1) { # label title only
          plot(1:nt,out.avg[,1],type="l",ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
          lines(1:nt,out.avg[,2],col="gray")
        } else if(i == length(ns)) { # label x axis only
          plot(1:nt,out.avg[,1],type="l",ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
          lines(1:nt,out.avg[,2],col="gray")
        } else { # label nothing
          plot(1:nt,out.avg[,1],type="l",ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab="",cex.axis=cex.axis)
          lines(1:nt,out.avg[,2],col="gray")
        }
        y = out.avg[,1]
        polygon(c(x[length(x)],x[1],x[1],x,x[length(x)]),c(ymins[k],ymins[k],y[1],y,y[length(y)]),col="gray",border=NA)
        y = out.avg[,2]
	  polygon(c(x[length(x)],x[1],x[1],x,x[length(x)]),c(ymaxs[k],ymaxs[k],y[1],y,y[length(y)]),col="gray",border=NA)
      } else { # lines only
        load(paste(dpath,"PF-quant-KD-normal-",resamps[j],"-",ns[i],".rdata",sep=""))
        out = pf.quant.out$theta.quant
        tt = dim(out)[1]; nt = tt - 1
        lines(1:nt,out[-1,k,2],col=cols[j])
        lines(1:nt,out[-1,k,3],col=cols[j])
      }
    }
    # Whiten borders
    edgex = 0
    edgey = 1
    polygon(c(x[length(x)],x[1],x[1],x[length(x)]),c(ymins[k],ymins[k],ymins[k]-edgey,ymins[k]-edgey),col="white",border=NA) # bottom border
    polygon(c(x[length(x)],x[1],x[1],x[length(x)]),c(ymaxs[k],ymaxs[k],ymaxs[k]+edgey,ymaxs[k]+edgey),col="white",border=NA) # top border
    polygon(c(x[length(x)],x[length(x)]+edgex,x[length(x)]+edgex,x[length(x)]),c(ymins[k]-edgey,ymins[k]-edgey,ymaxs[k]+edgey,ymaxs[k]+edgey),col="white",border=NA) # right border
    polygon(c(x[1],x[1]-edgex,x[1]-edgex,x[1]),c(ymins[k]-edgey,ymins[k]-edgey,ymaxs[k]+edgey,ymaxs[k]+edgey),col="white",border=NA) # left border
    box()
    if(k == 1 & i == 1) # add legend
    {
      legend("topright",legend=resamps,col=cols,lty=rep(1,4),cex=cex.leg,bg="white")
    }
  }
}
dev.off()

## Figure 4 - Extended model: KD PF for n particles with normal priors and stratified resampling

# Load data and get points where resampling done
n = 40000
load(paste(dpath,"sim-xy-ext.rdata",sep=""))
tt = dim(sim$x)[2]; nt = tt - 1
resampled = rep(0,nt)
load(paste(dpath,"PF-ext-KD-normal-stratified-",n,".rdata",sep=""))
parents = pf.out$out$parent
for(i in 1:125) resampled[i] = !all(parents[,i] == 1:n)
spts = which(as.logical(resampled))

# Set graphical parameters
params = expression(beta,gamma,nu,b,varsigma,sigma)
ymins = c(0.15,0.08,0.85,0,.8,0)
ymaxs = c(0.35,0.15,1.4,.5,1.2,.002)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4

# Construct plot
pdf(paste(gpath,"PF-ext-KD-stratified-normal-",n,".pdf",sep=""),width=30,height=30)
par(mfrow=c(3,3),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(k in 1:length(params))
{
  load(paste(dpath,"PF-quant-ext-KD-normal-stratified-",n,".rdata",sep=""))
  out = pf.quant.out$theta.quant
  tt = dim(out)[1]; nt = tt - 1
  if(k == 1) # label y axis, title, legend
  {
     plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=4,xlab="",ylab=paste("J = ",n,sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,3],col=4)
     points(spts,rep(ymins[k],length(spts)),pch="|",cex=2)
     myexpr = expression(paste("Known ",b,", ",varsigma,", ",sigma,sep=""))
     legend("topright",legend=c("Truth",myexpr,"Ext. Model"),col=c("gray47",2,4),lty=c(1,1,1),cex=cex.leg)
  } else { # label title only
     plot(1:nt,out[-1,k,2],type="l",ylim=c(ymins[k],ymaxs[k]),col=4,xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,3],col=4)
     points(spts,rep(ymins[k],length(spts)),pch="|",cex=2)
  }
  load(paste(dpath,"sim-xy-ext.rdata",sep=""))
  theta = c(theta,b,varsigma,sigma)
  abline(h=theta[k],col="gray47")
  if(k %in% 1:3)
  {
    load(paste(dpath,"PF-quant-2-KD-normal-stratified-",n,".rdata",sep=""))
    out = pf.quant.out$theta.quant
    tt = dim(out)[1]; nt = tt - 1
    lines(1:nt,out[-1,k,2],col=2)
    lines(1:nt,out[-1,k,3],col=2)
  }
}
params = expression(i,s,r)
ymin = 0; ymax = 1
for(k in 2:1)
{
  load(paste(dpath,"PF-quant-ext-KD-normal-stratified-",n,".rdata",sep=""))
  out = pf.quant.out$state.quant
  tt = dim(out)[1]; nt = tt - 1
  if(k == 2) # label x axis, title
  {
     plot(1:nt,out[-1,k,2],type="l",ylim=c(ymin,1),col=4,xlab="Time (days)",ylab="",main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,3],col=4)
     lines(1:nt,sim$x[k,-1],col="gray47")
     points(spts,rep(0,length(spts)),pch="|",cex=2)
  } else { # label title only
     plot(1:nt,out[-1,k,2],type="l",ylim=c(0,ymax),col=4,xlab="Time (days)",ylab="",main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,3],col=4)
     lines(1:nt,sim$x[k,-1],col="gray47")
     points(spts,rep(0,length(spts)),pch="|",cex=2)
  }
  load(paste(dpath,"PF-quant-2-KD-normal-stratified-",n,".rdata",sep=""))
  out = pf.quant.out$state.quant
  tt = dim(out)[1]; nt = tt - 1
  lines(1:nt,out[-1,k,2],col=2)
  lines(1:nt,out[-1,k,3],col=2)
}
k = 3
load(paste(dpath,"PF-quant-ext-KD-normal-stratified-",n,".rdata",sep=""))
out = pf.quant.out$state.quant
tt = dim(out)[1]; nt = tt - 1
truex = 1 - apply(sim$x[,-1],2,sum)
plot(1:nt,out[-1,k,2],type="l",ylim=c(0,ymax),col=4,xlab="Time (days)",ylab="",main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
lines(1:nt,out[-1,k,3],col=4)
lines(1:nt,truex,col="gray47")
points(spts,rep(0,length(spts)),pch="|",cex=2)
load(paste(dpath,"PF-quant-2-KD-normal-stratified-",n,".rdata",sep=""))
out = pf.quant.out$state.quant
tt = dim(out)[1]; nt = tt - 1
lines(1:nt,out[-1,k,2],col=2)
lines(1:nt,out[-1,k,3],col=2)
dev.off()