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

# pf.scat - function to resample particles, compute maximum/minimum values, and compute correlations between two parameters at specified time points; returns a list with components xrw1 and xrw2 that are M x length(cutoff) matrices of resampled particles (row = particle, column = time point), xlim and ylim that are each two element vectors that are the minimum and maximum values over all time points of xrw1 and xrw2, respectively, and r that is a length(cutoff) length vector of correlations between xrw1 and xrw2 at each time point
# Arguments:
# out - 2 by n by tt array of particles (assumed on original scale)
# wts - n by tt matrix of particle weights
# cutoff - k length vector of times at which to plot histograms
# M - number of samples to use to resample particles to get equal weights for constructing histograms
# xlim - length 2 vector giving xlim argument in plot(); if any elements NA, are calculated by default
# ylim - length 2 vector giving ylim argument in plot(); if any elements NA, are calculated by default
# ... - additional arguments passed to resample()
pf.scat = function(out, wts, cutoff, M=500, xlim=NA, ylim=NA, ...)
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
    tmps[,i] = resample(wts[,cutoff[i]], M, nonuniformity="none", ...)$indices
  }

  # Find minimum and maximum values so scatterplots on same scale
  xrw1 = matrix(NA,nr=M,nc=length(cutoff))
  xrw2 = matrix(NA,nr=M,nc=length(cutoff))
  if(!all(!is.na(xlim)))
  {
    xlim = c(Inf,-Inf)
    for(i in 1:length(cutoff))
    {
      xrw1[,i] = out[1,tmps[,i],cutoff[i]]
      xlim[1] = min(xlim[1],xrw1[,i])
      xlim[2] = max(xlim[2],xrw1[,i])
    }
  } else {
    for(i in 1:length(cutoff))
    {
      xrw1[,i] = out[1,tmps[,i],cutoff[i]]
    }
  }
  if(!all(!is.na(ylim)))
  {
    ylim = c(Inf,-Inf)
    for(i in 1:length(cutoff))
    {
      xrw2[,i] = out[2,tmps[,i],cutoff[i]]
      ylim[1] = min(ylim[1],xrw2[,i])
      ylim[2] = max(ylim[2],xrw2[,i])
    }
  } else {
    for(i in 1:length(cutoff))
    {
      xrw2[,i] = out[2,tmps[,i],cutoff[i]]
    }
  }

  # Compute correlation between particles at every cutoff point
  r = rep(NA,length(cutoff))
  for(i in 1:length(cutoff))
  {
    r[i] = cov.wt(cbind(out[1,,cutoff[i]],out[2,,cutoff[i]]),wts[,cutoff[i]],cor=TRUE)$cor[1,2]
  }  

  # Return re-weighted particles, correlations, and min/max values of particles
  return(list(xrw1=xrw1,xrw2=xrw2,xlim=xlim,ylim=ylim,r=r))
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