source("sir_functions.r")
source("ss_sim.r")

# Initialize .Random.seed
set.seed(sample(1:1000,1))

# Set graphics and data paths, param label
gpath = "../graphs/"
dpath = "../data/"

# Set known parameter values
P = 5000
d = 5
b = c(.25, .27, .23, .29)
varsigma = c(1.07, 1.05, 1.01, .98)
sigma = c(.0012, .0008, .0010, .0011)
dpower = 2

# Set unknown parameter values
theta = c(0.2399, 0.1066, 1.2042)

# Functions to simulate epidemic
revo_sim = function(x){ revo(x, P, d, theta, FALSE)}
robs_sim = function(x){ robs(x, b, varsigma, sigma, dpower)}
rinit_sim = function(){ rinit(10/P)}

# sims - function to simulate many epidemics and store the output
sims = function(nt)
{
  # Simulate epidemic
  sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
  return(sim)
}

# Simulate N times
N = 100
nt = 125
sims_y = array(NA,dim=c(N,4,nt))
for(i in 1:N) sims_y[i,,] = sims(nt)$y

# Which parameters unknown?
p = 1:3
s = rep(0,3); s[p] = 1
sinv = rep(1,3); sinv[p] = 0

# pf.meds - function to calculate quantiles of medians of filtered distributions
pf.meds = function(n, filt, resamp, prior, ...)
{
  # Create function to sample from prior distribution of theta and map theta to original scale
  if(prior == "uniform")
  {
    thetal = c(0.1400, 0.0900, 0.9500)
    thetau = c(0.5000, 0.1430, 1.3000)
    rtheta = function(){ u2theta(runif(3,thetal,thetau),thetal,thetau)}
    ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
  } else {
    theta.mean = c(-1.3296, -2.1764, 0.1055)
    theta.sd = sqrt(c(0.1055, 0.0140, 0.0064))
    rtheta = function(){ rnorm(3,theta.mean,theta.sd)}
    ftheta = function(theta,param=1) exp(theta)
  }

  # Functions to run the particle filters
  if(filt == "BF")
  {
    # Run bootstrap filter
    dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
    revo_bf = function(x){ revo(x, P, d, s*ftheta(x[p+2],p)+sinv*theta)}
    rprior_bf = function()
    { 
      myprior = rprior(sim$y[,ind], rtheta, b, varsigma, sigma, dpower)
      return(c(myprior$x,myprior$theta))
    }
    source("bf.r")
  } else if(filt == "APF"){
    # Run auxiliary particle filter
    dllik_apf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
    pstate_apf = function(x) { revo(x, P, d, s*ftheta(x[p+2],p)+sinv*theta, FALSE)}
    revo_apf = function(x){ revo(x, P, d, s*ftheta(x[p+2],p)+sinv*theta)}
    rprior_apf = function()
    { 
      myprior = rprior(sim$y[,ind], rtheta, b, varsigma, sigma, dpower)
      return(c(myprior$x,myprior$theta))
    }
    source("apf.r")
  } else {
    # Run kernel density particle filter
    dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma, dpower)}
    pstate_kd = function(x, mytheta) { revo(x, P, d, s*ftheta(mytheta,p)+sinv*theta, FALSE)}
    revo_kd = function(x, mytheta){ revo(x, P, d, s*ftheta(mytheta,p)+sinv*theta)}
    rprior_kd = function(){ rprior(sim$y[,ind], rtheta, b, varsigma, sigma, dpower)}
    source("kd_pf.r")
  }

  # run particle filters with n particles for each simulated data set
  np = length(p)
  ns = 2
  no = dim(sims_y)[2]
  nt = dim(sims_y)[3]
  N = dim(sims_y)[1]
  pf.med = array(NA,dim=c(N,nt+1,ns+np))
  require(Hmisc)
  for(j in 1:N)
  {
    # Get simulated data
    y = sims_y[j,,]

    # Get index of first non-empty observation
    empty = TRUE; ind = 1
    while(empty) if(all(is.na(y[,ind]))) ind = ind + 1 else empty = FALSE 
 
    if(filt == "BF")
    {
      # Run bootstrap filter
      out = bf(y, dllik_bf, revo_bf, rprior_bf, n, progress=FALSE, method=resamp, log=F, ...)
    } else if(filt == "APF"){
      # Run auxilliary particle filter
      out = apf(y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=FALSE, method=resamp, log=F, ...)
    } else {
      # Run kernel density particle filter
      out = kd_pf(y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=FALSE, method=resamp, log=F, ...)
    }

    # Calculate median of states and unknown parameters
    for(k in 1:(ns+np))
    {  
      for(i in 1:(nt+1))
      { 
        if(k <= 2)
        {
          if(filt != "KD")
          { 
            pf.med[j,i,k] = wtd.quantile(out$state[k,,i], out$weight[,i], normwt=T, probs=.5)
          } else {
            pf.med[j,i,k] = wtd.quantile(out$state[k,,i], out$weight[,i], normwt=T, probs=.5)
          }
        } else {
          if(filt == "KD")
          {
            pf.med[j,i,k] = wtd.quantile(ftheta(out$state[k,,i],k-ns), out$weight[,i], normwt=T, probs=.5)
          } else {
	    pf.med[j,i,k] = wtd.quantile(ftheta(out$theta[k-ns,,i],k-ns), out$weight[,i], normwt=T, probs=.5)
          }
        }
      }
    }
    cat("\n",j,"\n")
  }

  # Calculate median of the medians and lower/upper quantiles
  pf.m = pf.l = pf.u = matrix(NA,nr=nt+1,nc=ns+np)
  for(k in 1:(ns+np))
  {
    pf.m[,k] = apply(pf.med[,,k],2,median)
    pf.l[,k] = apply(pf.med[,,k],2,function(x) quantile(x,probs=.025))
    pf.u[,k] = apply(pf.med[,,k],2,function(x) quantile(x,probs=.975))
  }

  # Output
  pf.med.quant = array(NA,dim=c(nt+1,ns+np,3)
  pf.med.quant[,,1] = pf.m
  pf.med.quant[,,2] = pf.l
  pf.med.quant[,,3] = pf.u
  return(pf.med.quant)
}

# Run particle filters and save quantiles of medians
mydata = expand.grid(n = 100, filt = c("BF","APF","KD"), resamp = "stratified", prior = c("normal","uniform"), stringsAsFactors = FALSE)
med.quants = dlply(mydata,pf.meds)

# Which filters to plot?
filt = as.character(unique(mydata$filt))
ns = sort(unique(mydata$n), decreasing=FALSE)
resamps = as.character(unique(mydata$resamp))
priors = as.character(unique(mydata$prior))

# pf.plots - function to construct quantile plots
pf.plots = function(n, resamp, prior, label, truth=TRUE, ...)
{
  p = pf.quant.out$p
  probs = c(.5,.025,.975)
  nq = length(probs)
  nt = dim(med.quants)[[1]][1]

  # Load state and parameter quantile data
  state.quant = list(); length(state.quant) = length(filt) + truth
  param.quant = list(); length(param.quant) = length(filt) + truth
  index = 0
  if(truth)
  {
    simx = sims(nt)$x
    state.quant[[1]] = array(NA,dim=c(nt,2,1))
    state.quant[[1]][,1,1] = simx[1,-1]
    state.quant[[1]][,2,1] = simx[2,-1]
    param.quant[[1]] = array(NA,dim=c(nt,length(p),1))
    for(i in 1:length(p)) param.quant[[1]][,i,1] = rep(theta[p[i]],nt)
    index = 1
  }
  mq = med.quants[which(priors == prior & ns == n & resamps == resamp)]
  for(i in 1:length(filt))
  {
    sq = mq[[i]][,1:2,]
    tq = mq[[i]][,3:(2+length(p)),]
    state.quant[[i+index]] = array(NA,dim=c(nt,2,nq))
    param.quant[[i+index]] = array(NA,dim=c(nt,length(p),nq))
    state.quant[[i+index]] = sq[-1,,]
    param.quant[[i+index]] = tq[-1,,]
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
  file = paste(gpath,"PF-states-",prior,"-",resamp,"-",n,"-",label,".pdf",sep="")
  pf.plot(state.quant,1:nt,1,2,ymins=ymins,ymaxs=ymaxs,mlabs=mlabs.st,xlabs=xlabs.st,ylabs=ylabs.st,col=col,lty=lty,labs=labs,col.leg=col.leg,lty.leg=lty.leg,legend=leg,file=file)

  # Plot 95% credible bounds and medians of parameters over time
  ymins=rep(NA,length(p))
  ymaxs=rep(NA,length(p))
  file = paste(gpath,"PF-params-",prior,"-",resamp,"-",n,"-",label,".pdf",sep="")
  pf.plot(param.quant,1:nt,1,length(p),ymins=ymins,ymaxs=ymaxs,mlabs=mlabs.p,xlabs=xlabs.p,ylabs=ylabs.p,col=col,lty=lty,labs=labs,col.leg=col.leg,lty.leg=lty.leg,file=file,legend=leg,width=5*length(p),cex.lab=2.2,cex.main=2.4,legendsize=1.5)
}

# Construct quantile plots
require(plyr)
mydata2 = expand.grid(n=ns,resamp=resamps,prior=priors,label=c("n","resamp","prior"),stringsAsFactors=FALSE)
m_ply(mydata2,pf.plots)