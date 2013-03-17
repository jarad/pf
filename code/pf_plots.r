source("pf_functions.r")

# Set graphics and data path
gpath = "../graphs/"
dpath = "../data/"

# Figure 1 - Compare particle filters over different # particles
# load data frame
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("BF","APF","KD"), resamp = "systematic", prior = "uniform", nonuniformity="ess", threshold=0.8, stringsAsFactors=FALSE)
#mydata = expand.grid(n = 10000, filt = "KD", resamp = "stratified", prior = "normal", nonuniformity="ess", threshold=0.8, stringsAsFactors=FALSE)

# Which filters to plot?
filts = as.character(unique(mydata$filt))
ns = sort(unique(mydata$n), decreasing=FALSE)
resamps = as.character(unique(mydata$resamp))
priors = as.character(unique(mydata$prior))

# Which parameters to plot
p = 1:3

# Which model label
mod = ""

# Which parameters label
ps = ""

# Y-axis limits 
ymins.p = c(0.15,0.08,0.8,0,.75,0)
ymaxs.p = c(0.35,0.16,1.45,1,1.25,.0035)

# pf.plots - function to construct quantile plots for comparing PF methods
pf.plots = function(n, resamp, prior, truth=TRUE, ...)
{
  # Which quantiles to plot (probs = c(.5, .025, .975))
  probs_ind = 2:3
  nq = length(probs_ind)

  # Load state and parameter quantile data
  state.quant = list(); length(state.quant) = length(filts) + truth
  param.quant = list(); length(param.quant) = length(filts) + truth
  index = 0
  if(truth)
  {
    load(paste(dpath,"sim-xy",mod,".rdata",sep=""))
    theta = c(theta,b,varsigma,sigma)
    state.quant[[1]] = array(NA,dim=c(nt,3,1))
    state.quant[[1]][,1,1] = sim$x[1,-1]
    state.quant[[1]][,2,1] = sim$x[2,-1]
    state.quant[[1]][,3,1] = 1 - sim$x[1,-1] - sim$x[2,-1]
    param.quant[[1]] = array(NA,dim=c(nt,length(p),1))
    for(i in 1:length(p)) param.quant[[1]][,i,1] = rep(theta[p[i]],nt)
    index = 1
  }
  for(i in 1:length(filts))
  {
    load(paste(dpath,"PF-quant",mod,"-",filts[i],"-",prior,"-",resamp,"-",n,".rdata",sep=""))
    state.quant[[i+index]] = array(NA,dim=c(nt,3,nq))
    param.quant[[i+index]] = array(NA,dim=c(nt,length(p),nq))
    state.quant[[i+index]][,1:2,] = pf.quant.out$state.quant[-1,,probs_ind]
    state.quant[[i+index]][,3,] = 1 - pf.quant.out$state.quant[-1,1,probs_ind] - pf.quant.out$state.quant[-1,2,probs_ind]
    param.quant[[i+index]] = pf.quant.out$theta.quant[-1,p,probs_ind]
  }

  # Create graph labels
  expr = expression(beta,gamma,nu,b,varsigma,sigma)
  leg = FALSE
  if(n == ns[1])
  {
    mlabs.st = expression(i,s,r)
    xlabs.st = c("Time (days)","","")
    mlabs.p = expr[p]
    xlabs.p = c("Time (days)",rep("",length(p)-1))
    leg = TRUE
  } else {
    mlabs.st = rep("",3)
    xlabs.st = rep("",3)
    mlabs.p = rep("",length(p))
    xlabs.p = rep("",length(p))
  }
  ylabs.st = c(paste("J = ",n,sep=""),"","")
  ylabs.p = c(paste("J = ",n,sep=""),rep("",length(p)-1))
  col = rep(NA,length(filts)+truth)
  index = 0
  if(truth)
  { 
    col[1] = 1
    index = 1
  }
  if("BF" %in% filts)
  { 
    col[index + which(filts == "BF")] = 2
  }
  if("APF" %in% filts)
  { 
    col[index + which(filts == "APF")] = 4
  }
  if("KD" %in% filts)
  { 
    col[index + which(filts == "KD")] = 3
  }
  lty = rep(1,nq)
  labs = filts
  col.leg = col
  lty.leg = rep(1,length(col.leg))
  if(truth) labs = c("Truth",labs)

  # Plot % infected and % susceptible quantiles over time
  ymins=c(0,NA,0)
  ymaxs=c(NA,1,NA)
  file = paste(gpath,"PF-states",mod,"-",prior,"-",resamp,"-",n,".pdf",sep="")
  pf.plot(state.quant,1:nt,1,3,ymins=ymins,ymaxs=ymaxs,mlabs=mlabs.st,xlabs=xlabs.st,ylabs=ylabs.st,col=col,lty=lty,labs=labs,col.leg=col.leg,lty.leg=lty.leg,legend=leg,file=file,width=15,cex.lab=2.2,cex.main=2.4,legendsize=1.5)

  # Plot 95% credible bounds and medians of parameters over time
  ymins = ymins.p[p]
  ymaxs = ymaxs.p[p]
  file = paste(gpath,"PF-params",mod,ps,"-",prior,"-",resamp,"-",n,".pdf",sep="")
  pf.plot(param.quant,1:nt,1,length(p),ymins=ymins,ymaxs=ymaxs,mlabs=mlabs.p,xlabs=xlabs.p,ylabs=ylabs.p,col=col,lty=lty,labs=labs,col.leg=col.leg,lty.leg=lty.leg,file=file,legend=leg,width=5*length(p),cex.lab=2.2,cex.main=2.4,legendsize=1.5)
}

# Construct quantile plots
require(plyr)
mydata2 = expand.grid(n=ns,resamp=resamps,prior=priors,stringsAsFactors=FALSE)
m_ply(mydata2,pf.plots)

# Figure 3 - Compare resampling schemes over different # particles
# load data frame
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = "KD", resamp = c("multinomial","residual","stratified","systematic"), prior = "normal", nonuniformity="ess", threshold=0.8, stringsAsFactors=FALSE)

# Which filters to plot?
filts = as.character(unique(mydata$filt))
ns = sort(unique(mydata$n), decreasing=FALSE)
resamps = as.character(unique(mydata$resamp))
priors = as.character(unique(mydata$prior))

# pf.plots.resamp - function to construct quantile plots, compare resamp schemes
pf.plots.resamp = function(filt, prior, n, ...)
{
  # Which quantiles to plot
  prob_ind = 2:3
  nq = length(prob_ind)

  # Load state and parameter quantile data
  state.quant = list(); length(state.quant) = length(resamps) + 1
  param.quant = list(); length(param.quant) = length(resamps) + 1
  load(paste(dpath,"PF-quant",mod,"-",filt,"-",prior,"-",resamps[1],"-",20000,".rdata",sep=""))
  res.quants = pf.quant.out$theta.quant[-1,p,prob_ind]
  nt = dim(res.quants)[1]
  if(length(resamps) > 1)
  {
    for(i in 2:length(resamps))
    {
      load(paste(dpath,"PF-quant",mod,"-",filt,"-",prior,"-",resamps[i],"-",20000,".rdata",sep=""))
      res.quants = res.quants + pf.quant.out$theta.quant[-1,p,prob_ind]
    }
  }
  param.quant[[1]] = res.quants / length(resamps)
  for(i in 1:length(resamps))
  {
    load(paste(dpath,"PF-quant",mod,"-",filt,"-",prior,"-",resamps[i],"-",n,".rdata",sep=""))
    param.quant[[i+1]] = pf.quant.out$theta.quant[-1,p,prob_ind]
  }

  # Create graph labels
  expr = expression(beta,gamma,nu,b,varsigma,sigma)
  leg = FALSE
  if(n == ns[1])
  {
    mlabs.p = expr[p]
    xlabs.p = c("Time (days)",rep("",length(p)-1))
    leg = TRUE
  } else {
    mlabs.p = rep("",length(p))
    xlabs.p = rep("",length(p))
  }
  ylabs.p = c(paste("J = ",n,sep=""),rep("",length(p)-1))
  col = rep(NA,length(resamps)+1)
  col[1] = 1
  if("multinomial" %in% resamps)
  { 
    col[1 + which(resamps == "multinomial")] = 2
  }
  if("residual" %in% resamps)
  { 
    col[1 + which(resamps == "residual")] = 3
  }
  if("stratified" %in% resamps)
  { 
    col[1 + which(resamps == "stratified")] = 4
  }
  if("systematic" %in% resamps)
  { 
    col[1 + which(resamps == "systematic")] = 5
  }
  lty = rep(1,nq)
  labs = c("True Posterior",resamps)
  col.leg = col
  lty.leg = rep(1,length(col.leg))

  # Plot 95% credible bounds and medians of parameters over time
  ymins = ymins.p[p]
  ymaxs = ymaxs.p[p]
  file = paste(gpath,"PF-params",mod,ps,"-",prior,"-",filt,"-",n,".pdf",sep="")
  pf.plot(param.quant,1:nt,1,length(p),ymins=ymins,ymaxs=ymaxs,pgon=TRUE,mlabs=mlabs.p,xlabs=xlabs.p,ylabs=ylabs.p,col=col,lty=lty,labs=labs,col.leg=col.leg,lty.leg=lty.leg,file=file,legend=leg,width=5*length(p),cex.lab=2.2,cex.main=2.4,legendsize=1.5)
}

# Construct quantile plots
require(plyr)
mydata2 = expand.grid(filt=filts,n=ns,prior=priors,stringsAsFactors=FALSE)
m_ply(mydata2,pf.plots.resamp)

# Figure 2 - Create scatterplots of beta v gamma over time
# load data frame
mydata = expand.grid(n = 10000, filt = "KD", resamp = "systematic", prior = c("normal","uniform"), stringsAsFactors=FALSE)

# pf.scats - function to construct scatterplots of beta v gamma
pf.scats = function(n, filt, resamp, prior, ...)
{
  # Load data
  load(paste(dpath,"sim-xy",mod,".rdata",sep=""))
  load(paste(dpath,"PF",mod,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  if(filt == "KD") out = pf.out$out$theta else out = pf.out$out$state[-(1:2),,]
  ftheta = pf.out$ftheta
  betas = ftheta(out[1,,],1)
  gammas = ftheta(out[2,,],2)

  # Scatterplot
  expr = expression(beta,gamma)
  myout = array(NA,dim=c(2,dim(betas)[1],dim(betas)[2]))
  myout[1,,] = betas
  myout[2,,] = gammas
  file = paste(gpath,"Hist",mod,"-",filt,"-",prior,"-",resamp,"-",n,"-betagamma.pdf",sep="")
  cutoff = seq(16,121,len=8)
  pf.scat(myout,pf.out$out$weight,cutoff,expr[1:2],theta[1:2],file,M=500,mr=2,mc=4,method="stratified")
}

# Construct scatterplots
m_ply(mydata,pf.scats)

###################################################################
# Extra plots - histograms of marginal posteriors and contour plots
###################################################################

# pf.hists - function to construct histograms over time
pf.hists = function(n, filt, resamp, prior, ...)
{
  # Load data
  load(paste(dpath,"sim-xy",mod,".rdata",sep=""))
  theta = c(theta,b,varsigma,sigma)
  load(paste(dpath,"PF",mod,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta

  # Create file labels
  parnames = c("beta","gamma","nu","b","varsigma","sigma")
  expr = expression(beta,gamma,nu,b,varsigma,sigma)
  files = c()
  for(i in p) files = c(files,paste(gpath,"Hist",mod,"-",filt,"-",prior,"-",resamp,"-",n,"-",parnames[i],".pdf",sep=""))

  # Create histograms
  cutoff = seq(16,121,len=8)
  if(filt == "KD") myout = out$theta else myout = out$state[2+p,,]
  pf.hist(myout,out$weight,cutoff,ftheta,expr[p],theta[p],files,mr=2,mc=4,method="stratified")
}

# Construct histograms
m_ply(mydata,pf.hists)

# pf.contours - function to construct contour plots of beta v gamma
pf.contours = function(n, filt, resamp, prior, ...)
{
  # Load data
  load(paste(dpath,"sim-xy",mod,".rdata",sep=""))
  load(paste(dpath,"PF",mod,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  if(filt == "KD") out = pf.out$out$theta else out = pf.out$out$state[-(1:2),,]
  ftheta = pf.out$ftheta
  betas = ftheta(out[1,,],1)
  gammas = ftheta(out[2,,],2)

  # Scatterplot
  expr = expression(beta,gamma)
  myout = array(NA,dim=c(2,dim(betas)[1],dim(betas)[2]))
  myout[1,,] = betas
  myout[2,,] = gammas
  file = paste(gpath,"Hist-contour",mod,"-",filt,"-",prior,"-",resamp,"-",n,"-betagamma.pdf",sep="")
  cutoff = seq(16,121,len=8)
  pf.contour(myout,pf.out$out$weight,cutoff,expr[1:2],theta[1:2],file,mr=2,mc=4,method="stratified")
}

# Construct contour plots
m_ply(mydata,pf.contours)