source("pf_functions.r")

# Set graphics and data path
gpath = "../graphs/"
dpath = "../data/"

# load data frame
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("BF","APF","KD"), resamp = c("multinomial","residual","stratified","systematic"), prior = c("normal","uniform"), nonuniformity="ess", threshold=0.8, stringsAsFactors=FALSE)

# Which filters to plot?
filt = as.character(unique(mydata$filt))
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
ymins.p = c(0.1,0.08,0.8)
ymaxs.p = c(0.55,0.16,1.45)

# pf.plots - function to construct quantile plots
pf.plots = function(n, resamp, prior, label, truth=TRUE, ...)
{
  # Which parameters to plot
  load(paste(dpath,"PF-quant",mod,"-",filt[1],"-",priors[1],"-",resamps[1],"-",ns[1],".rdata",sep=""))
  probs = pf.quant.out$probs
  nq = length(probs)

  # Load state and parameter quantile data
  state.quant = list(); length(state.quant) = length(filt) + truth
  param.quant = list(); length(param.quant) = length(filt) + truth
  index = 0
  if(truth)
  {
    load(paste(dpath,"sim-xy",mod,".rdata",sep=""))
    theta = c(theta,b,varsigma,sigma)
    state.quant[[1]] = array(NA,dim=c(nt,2,1))
    state.quant[[1]][,1,1] = sim$x[1,-1]
    state.quant[[1]][,2,1] = sim$x[2,-1]
    param.quant[[1]] = array(NA,dim=c(nt,length(p),1))
    for(i in 1:length(p)) param.quant[[1]][,i,1] = rep(theta[p[i]],nt)
    index = 1
  }
  for(i in 1:length(filt))
  {
    load(paste(dpath,"PF-quant",mod,"-",filt[i],"-",prior,"-",resamp,"-",n,".rdata",sep=""))
    state.quant[[i+index]] = array(NA,dim=c(nt,2,nq))
    param.quant[[i+index]] = array(NA,dim=c(nt,length(p),nq))
    state.quant[[i+index]] = pf.quant.out$state.quant[-1,,]
    param.quant[[i+index]] = pf.quant.out$theta.quant[-1,p,]
  }

  # Create graph labels
  expr = expression(beta,gamma,nu,b,varsigma,sigma)
  leg = FALSE
  if(label == "n")
  {
    if(n == ns[1])
    {
      mlabs.st = expression(i,s)
      xlabs.st = c("Time (days)","")
      mlabs.p = expr[p]
      xlabs.p = c("Time (days)",rep("",length(p)-1))
      leg = TRUE
    } else {
      mlabs.st = rep("",2)
      xlabs.st = c("","")
      mlabs.p = rep("",length(p))
      xlabs.p = rep("",length(p))
    }
    ylabs.st = c(paste("J = ",n,sep=""),"")
    ylabs.p = c(paste("J = ",n,sep=""),rep("",length(p)-1))
  } else if(label == "resamp") {
    if(resamp == resamps[1])
    {
      mlabs.st = expression(i,s)
      xlabs.st = c("Time (days)","")
      mlabs.p = expr[p]
      xlabs.p = c("Time (days)",rep("",length(p)-1))
      leg = TRUE
    } else {
      mlabs.st = rep("",2)
      xlabs.st = c("","")
      mlabs.p = rep("",length(p))
      xlabs.p = rep("",length(p))
    }
    ylabs.st = c(resamp,"")
    ylabs.p = c(resamp,rep("",length(p)-1))
  } else if(label == "prior") {
    if(prior == priors[1])
    {
      mlabs.st = expression(i,s)
      xlabs.st = c("Time (days)","")
      mlabs.p = expr[p]
      xlabs.p = c("Time (days)",rep("",length(p)-1))
      leg = TRUE
    } else {
      mlabs.st = rep("",2)
      xlabs.st = c("","")
      mlabs.p = rep("",length(p))
      xlabs.p = rep("",length(p))
    }
    ylabs.st = c(paste(prior,"prior",sep=" "),"")
    ylabs.p = c(paste(prior,"prior",sep=" "),rep("",length(p)-1))
  }
  col = rep(NA,length(filt)+truth)
  index = 0
  if(truth)
  { 
    col[1] = 1
    index = 1
  }
  if("BF" %in% filt)
  { 
    col[index + which(filt == "BF")] = 2
  }
  if("APF" %in% filt)
  { 
    col[index + which(filt == "APF")] = 4
  }
  if("KD" %in% filt)
  { 
    col[index + which(filt == "KD")] = 3
  }
  if(nq < 2)
  {
    lty = 1 
    labs = filt
    col.leg = col
    lty.leg = rep(1,length(filt)+truth)
  } else {
    lty = c(1,rep(2,nq-1))
    bound = 100*(max(probs[-1]) - min(probs[-1]))
    labs = c(filt,paste(bound,"% bounds",sep=""))
    col.leg = c(col,1)
    lty.leg = c(rep(1,length(filt)+truth),2)
  }
  if(truth) labs = c("Truth",labs)

  # Plot % infected and % susceptible quantiles over time
  ymins=c(0,NA)
  ymaxs=c(NA,1)
  file = paste(gpath,"PF-states",mod,ps,"-",prior,"-",resamp,"-",n,"-",label,".pdf",sep="")
  pf.plot(state.quant,1:nt,1,2,ymins=ymins,ymaxs=ymaxs,mlabs=mlabs.st,xlabs=xlabs.st,ylabs=ylabs.st,col=col,lty=lty,labs=labs,col.leg=col.leg,lty.leg=lty.leg,legend=leg,file=file)

  # Plot 95% credible bounds and medians of parameters over time
#  ymins=rep(NA,length(p))
#  ymaxs=rep(NA,length(p))
  ymins = ymins.p[p]
  ymaxs = ymaxs.p[p]
  file = paste(gpath,"PF-params",mod,ps,"-",prior,"-",resamp,"-",n,"-",label,".pdf",sep="")
  pf.plot(param.quant,1:nt,1,length(p),ymins=ymins,ymaxs=ymaxs,mlabs=mlabs.p,xlabs=xlabs.p,ylabs=ylabs.p,col=col,lty=lty,labs=labs,col.leg=col.leg,lty.leg=lty.leg,file=file,legend=leg,width=5*length(p),cex.lab=2.2,cex.main=2.4,legendsize=1.5)
}

# Construct quantile plots
require(plyr)
mydata2 = expand.grid(n=ns,resamp=resamps,prior=priors,label=c("n","resamp","prior"),stringsAsFactors=FALSE)
m_ply(mydata2,pf.plots)

#p = 1:6

# pf.hists - function to construct histograms over time
pf.hists = function(n, filt, resamp, prior, ...)
{
  # Load data
  load(paste(dpath,"sim-xy",mod,".rdata",sep=""))
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
  pf.scat(myout,pf.out$out$weight,cutoff,expr[1:2],theta[1:2],file,mr=2,mc=4,method="stratified")
}

# Construct scatterplots
m_ply(mydata,pf.scats)

# pf.scats - function to construct scatterplots of beta v gamma
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

# Clear objects
rm(list=ls(all=TRUE))