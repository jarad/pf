# pf_functions - functions to help calculate quantiles of particle filtered distributions and construct plots
#
# pf.quantile - function that calculates quantiles of filtered distributions; returns a 3-D array with dimensions number of observations, number of parameters and/or states, and number of quantiles desired
# Arguments:
# out - an 3-D dimensional array of particles with dimensions number of parameters and/or states, number of particles, and number of observations; this can be the state or theta component of the output list from functions bf, apf, or kd_pf
# wts - a matrix of wts with rows corresponding to particles and columns corresponding to observations; can be the weight component of the output list from functions bf, apf, or kd_pf; number of rows must match the 2nd dimension of out and number of columns must match the third dimension of out
# ftheta - a function that maps a vector from one scale to another, e.g. exp(theta) to map theta from log(scale) to original scale; can vary with parameter
# probs - a vector of quantiles to be calculated
# normwt - TRUE or FALSE corresponding to the normwt argument in wtd.quantile
pf.quantile = function(out, wts, ftheta, probs=.5, normwt=TRUE)
{
  require(Hmisc)
  nq = length(probs)
  np = dim(out)[1]
  tt = dim(out)[3]
  if(dim(wts)[1] != dim(out)[2] | dim(wts)[2] != tt) stop("dimensions of out and wts don't match")
  quantiles = array(NA,dim=c(tt,np,nq))
  for(i in 1:tt)
  {
    for(j in 1:np)
    {
      quantiles[i,j,] = wtd.quantile(ftheta(out[j,,i], j), wts[,i], normwt=normwt, probs=probs)
    }
  }
  return(quantiles)
}

# pf.plot - function to construct plots of quantiles of particle filtered distributions
# Arguments:
# quantiles - list of arrays returned by pf.quantile
# x - x variable, vector of length equal to first dimension of elements of quantiles
# mr - number of rows in plot panel window
# mc - number of columns in plot panel window
# mar - corresponds to mar argument in par() for resizing plot windows
# ymins - a vector whose jth element corresponds to the minimum y-axis value for the jth plot panel; can have NAs
# ymaxs - a vector whose jth element corresponds to the maximum y-axis value for the jth plot panel; can have NAs
# mlabs - a vector of length equal to the second dimension of elements of quantiles giving title labels for each panel
# xlabs - a vector of length equal to the second dimension of elements of quantiles giving x-axis labels for each panel
# ylabs - a vector of length equal to the second dimension of elements of quantiles giving y-axis labels for each panel
# col - a vector of length equal to the length of quantiles; gives the color to be used for each particle filter
# lty - a vector length equal to the third dimension of elements of quantiles; gives the line type to be used for each quantile
# legend - TRUE/FALSE, should a legend be added to the upper left-most plot?
# labs - a vector giving labels for the legend
# col.leg - a vector giving colors for legend labels
# lty.leg - a vector giving line types for legend labels
# location - gives location for the legend; corresponds to first argument in legend()
# legendsize - scalar, gives size of legend labels
# file - name of outputted pdf file
# width, height - correspond to width and height arguments in pdf()
# ... - additional arguments to pass into plot()
pf.plot = function(quantiles, x, mr=1, mc=1, mar=c(5,5,3,1)+.1, ymins, ymaxs, mlabs, xlabs, ylabs, col, lty, legend=TRUE, labs="Truth", col.leg = 1, lty.leg = 1, location="topright", legendsize=1, file, width=10, height=5, ...)
{
  if(!is.list(quantiles)) stop("quantiles must be a list of 3-D arrays")
  nf = length(quantiles)
  if(nf < 1) stop("quantiles must have at least 1 element")
  np = dim(quantiles[[1]])[2]
  tt = dim(quantiles[[1]])[1]
  nt = tt - 1
  if(length(x) != tt) stop("x should be of length equal to first dimension of elements of quantiles")
  if(!(length(col) == nf)) stop("length of col should match length of quantiles")
  for(i in 1:np)
  {
    if(is.na(ymins[i])) ymins[i] = min(sapply(quantiles, function(x) min(x[,i,],na.rm=TRUE)))
    if(is.na(ymaxs[i])) ymaxs[i] = max(sapply(quantiles, function(x) max(x[,i,],na.rm=TRUE)))
  }
  pdf(file, width, height)
  par(mfrow=c(mr,mc),mar=mar)
  for(j in 1:np)
  {
    nq = dim(quantiles[[1]])[3]
    if(length(lty) < nq) stop("lty must have length at least equal to the largest of the 3rd dimensions of the elements of quantiles")
    plot(x,quantiles[[1]][,j,1],type="l",ylim=c(ymins[j],ymaxs[j]),col=col[1],lty=lty[1],xlab=xlabs[j],ylab=ylabs[j],main=mlabs[j],...)
    if(legend & j == 1) legend(location,labs,col=col.leg,lty=lty.leg,cex=legendsize)
    if(nq > 1)
    {
      for(i in 2:nq)
      {
       lines(x,quantiles[[1]][,j,i],col=col[1],lty=lty[i])
      }
    }
    if(nf > 1)
    {
      for(k in 2:nf)
      {
        nq = dim(quantiles[[k]])[3]
        if(length(lty) < nq) stop("lty must have length at least equal to the largest of the 3rd dimensions of the elements of quantiles")
        for(i in 1:nq)
        {
          lines(x,quantiles[[k]][,j,i],col=col[k],lty=lty[i])
        }
      }
    }
  }
  dev.off()
}

# pf.hist - function to plot histograms of unknown parameters at specified time points for particle filtered samples
# Arguments:
# out - np by n by tt array of particles
# wts - n by tt matrix of particle weights
# cutoff - k length vector of times at which to plot histograms
# M - number of samples to use to resample particles to get equal weights for constructing histograms
# ftheta - function that maps an n by tt matrix of particles to original scale and takes an additional integer argument to optionally depend on parameter
# plabs - vector of length np giving labels for parameters
# truth - np length vector of true parameter values
# tsize - size of subtitle displaying true parameter values
# mr - number of rows in plot panel window
# mc - number of columns in plot panel window
# mar - corresponds to mar argument in par() for resizing plot windows
# msize - size of titles of upper-left histogram
# labsize - size of axis labels on upper-left histogram
# file - np length character string vector of names of outputted pdf files
# width, height - arguments to pdf()
# ... - additional arguments passed to resample()
pf.hist = function(out,wts,cutoff,ftheta,plabs,truth,file,M=10000,tsize=.65,mr=1,mc=1,mar=c(4,5.2,3,.5)+.1,msize=1.5,labsize=1.5,width=10,height=5,...)
{
  require(smcUtils)

  # Check dimensions
  n = dim(wts)[1]
  tt = dim(wts)[2]
  np = dim(out)[1]
  if(n != dim(out)[2] | tt != dim(out)[3]) stop("dimensions of out and wts don't match")
  nd = length(cutoff)
  if(nd > tt) stop("Too many cutoff points")
  if(length(plabs) != np) stop("labels should be length equal to number of parameters")

  # Resample particles to get equal weights
  tmps = matrix(NA,nr=M,nc=length(cutoff))
  for(i in 1:length(cutoff))
  {
    tmps[,i] = resample(wts[,cutoff[i]], M, nonuniformity="none",...)$indices
  }

  # Plot histograms of unknown parameters over time
  for(j in 1:np)
  {
    x = ftheta(out[j,,cutoff],j)
    xrw = matrix(NA,nr=M,nc=length(cutoff))
    for(i in 1:length(cutoff)) xrw[,i] = x[tmps[,i],i]
    xmin = min(apply(xrw,2,function(x) min(hist(x,plot=FALSE)$breaks)))
    xmax = max(apply(xrw,2,function(x) max(hist(x,plot=FALSE)$breaks)))
    pdf(file[j],width=width,height=height)
    par(mfrow=c(mr,mc),mar=mar) # dimensions should depend on cutoff
    for(i in 1:length(cutoff))
    {
      if(i==1){
        hist(xrw[,i],xlim=c(xmin,xmax),freq=FALSE,main=paste("t = ",cutoff[i]-1,sep=""),xlab=plabs[j],cex.main=msize,cex.lab=labsize)
        abline(v=truth[j],col=2,lwd=2)
        mtext(paste("Truth = ",truth[j],sep=""),side=3,cex=tsize)
      } else {
        hist(xrw[,i],xlim=c(xmin,xmax),freq=FALSE,main=paste("t = ",cutoff[i]-1,sep=""),xlab="",ylab="",cex.main=msize)
        abline(v=truth[j],col=2,lwd=2)
      }
    }
    dev.off()
  }
}

# pf.scat - function to plot scatterplots of two parameters at specified time points
# Arguments:
# out - 2 by n by tt array of particles (assumed on original scale)
# wts - n by tt matrix of particle weights
# cutoff - k length vector of times at which to plot histograms
# M - number of samples to use to resample particles to get equal weights for constructing histograms
# plabs - vector of length 2 giving labels for parameters
# truth - 2 length vector of true parameter values
# tsize - size of subtitle displaying true parameter values
# mr - number of rows in plot panel window
# mc - number of columns in plot panel window
# mar - corresponds to mar argument in par() for resizing plot windows
# msize - size of titles of upper-left histogram
# labsize - size of axis labels on upper-left histogram
# file - name of outputted pdf file
# width, height - arguments to pdf()
# ... - additional arguments passed to resample()
pf.scat = function(out,wts,cutoff,plabs,truth,file,M=10000,tsize=.65,mr=1,mc=1,mar=c(4,5.2,3,.5)+.1,msize=1.5,labsize=1.5,width=10,height=5,...)
{
  require(smcUtils)

  # Check dimensions
  n = dim(wts)[1]
  tt = dim(wts)[2]
  if(n != dim(out)[2] | tt != dim(out)[3]) stop("dimensions of out and wts don't match")
  nd = length(cutoff)
  if(nd > tt) stop("Too many cutoff points")
  if(dim(out)[1] != 2) stop("Incorrect first dimension of out")

  # Resample particles to get equal weights
  tmps = matrix(NA,nr=M,nc=length(cutoff))
  for(i in 1:length(cutoff))
  {
    tmps[,i] = resample(wts[,cutoff[i]], M, nonuniformity="none",...)$indices
  }

  # Find minimum and maximum values so scatterplots on same scale
  xrw1 = matrix(NA,nr=M,nc=length(cutoff))
  xrw2 = matrix(NA,nr=M,nc=length(cutoff))
  xmin = Inf; xmax = -Inf; ymin = Inf; ymax = -Inf
  for(i in 1:length(cutoff))
  {
    xrw1[,i] = out[1,tmps[,i],cutoff[i]]
    xrw2[,i] = out[2,tmps[,i],cutoff[i]]
    xmin = min(xmin,xrw1[,i])
    xmax = max(xmax,xrw1[,i])
    ymin = min(ymin,xrw2[,i])
    ymax = max(ymax,xrw2[,i])
  }

  # Scatterplots over time
  pdf(file,width=width,height=height)
  par(mfrow=c(mr,mc),mar=mar)
  for(i in 1:length(cutoff))
  {
    if(i == 1)
    {
      plot(xrw1[,i],xrw2[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=plabs[1],ylab=plabs[2],main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize,cex.lab=labsize)
      mtext(paste("Truth = ",truth[1],sep=""),side=1,cex=tsize)
      mtext(paste("Truth = ",truth[2],sep=""),side=2,cex=tsize)
      points(truth[1],truth[2],col=2,pch=20)
    } else {
      plot(xrw1[,i],xrw2[,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize)
      points(truth[1],truth[2],col=2,pch=20)
    }
  }
  dev.off()
}

# Functions to reparameterize theta to [a,b] from the real line and vice versa
theta2u = function(theta,a,b)
{
  etheta = exp(theta)
  return((b*etheta + a) / (1 + etheta))
}
u2theta = function(u,a,b)
{
  U = (u - a) / (b - a)
  return(log(U / (1 - U)))
}