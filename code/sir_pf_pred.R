# Set graphics and data path
gpath = "../gpath/"
dpath = "../dpath/"

# Load simulated data and create revo function
load(paste(dpath,"sim-xy.rdata",sep=""))
pstate_revo = function(x) revo(x, P, d, s*ftheta(x[p+2],p)+sinv*theta, FALSE)

# Function needed to generate predicted SIR curves
ss.pred = function(nt,revo,rinit)
{
  # Find dimension of state
  current.seed = .Random.seed
  ns = length(rinit())
  .Random.seed = current.seed

  # Initialize state matrix
  x = matrix(NA,ns,nt+1)
  x[,1] = rinit()

  # Generate states and observations
  for(i in 1:nt)
  {
    x[,i+1] = revo(x[,i])
  }

  # Output
  return(x)
}

# pf.pred - function to perform prediction of epidemic
pf.pred = function(n, filt, resamp, prior, meth)
{
  # Load data
  load(paste(dpath,"PF-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta
  p = pf.out$p
  tt = dim(out$state)[3]
  nt = tt - 1

  # Initialize output vectors
  require(Hmisc)
  days <- seq(10,35,5)
  pf.li <- pf.i <- pf.ui <- list()
  pf.peak.prop <- pf.peak.calc <- list()
  length(pf.li) = length(pf.i) = length(pf.ui) = length(days)
  length(pf.peak.prop) = length(pf.peak.calc) = length(days)

  # Predict epidemic by propagating particles after being filtered for different numbers of days
  pb = txtProgressBar(0,n*length(days),style=3)
  for(j in 1:length(days))
  {
    # Grab particles from each filter at time days[j]
    out.theta = array(NA,dim=c(length(p),n,tt)
    if(filt != "KD"){
      for(j in p) out.theta[j,,] = ftheta(out$state[j+2,,],j+2)
    } else {
      for(j in p) out.theta[j,,] = ftheta(out$theta[j,,],j)
    }
    mystate = rbind(out$state[1:2,,days[j]+1],out.theta[,,days[j]+1])

    # Propagate particles forward to end of epidemic
    pred = array(NA,dim=c(2+length(p),n,nt-days[j]))
    for(i in 1:n)
    {
      setTxtProgressBar(pb,i + n*(j-1))
      pred[,i,] = ss.pred(nt-(days[j]+1),pstate_revo,function() return(mystate[,i]))
    }

    # Calculate .025 and .975 quantiles of predicted particle samples
    pf.li[[j]] = pf.i[[j]] = pf.ui[[j]] = rep(NA,nt-days[j])
    for(i in 1:(nt-days[j]))
    {
      pf.li[[j]][i] = wtd.quantile(pred[1,,i], out$weight[,days[j]+1], normwt=T, probs=.025)
      pf.i[[j]][i] = wtd.quantile(pred[1,,i], out$weight[,days[j]+1], normwt=T, probs=.5)
      pf.ui[[j]][i] = wtd.quantile(pred[1,,i], out$weight[,days[j]+1], normwt=T, probs=.975)
    }

    # Get peak times and intensities of particles
    # By propagation
    if(meth = "prop" | meth == "both")
    {
      get.peak.time = function(x){which(x == max(x))}
      peak_intensity = apply(pred[1,,],1,max)
      peak_time = apply(pred[1,,],1,get.peak.time) - 1 + days[j] 

      # Calculate .025, .5, and .975 quantiles of peak times and intensities
      pf.peak.prop[[j]] = matrix(NA,nr=2,nc=3)
      pf.peak.prop[[j]][1,] = wtd.quantile(peak_time[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
      pf.peak.prop[[j]][2,] = wtd.quantile(peak_intensity[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
    }
    # By calculation
    if(meth == "calc" | meth == "both")
    {
      calc.peak.time = function(x)
      {
        y = rep(NA,6); ind = 3
        y[1:2] = x[1:2]
        for(j in 1:3)
        {
          if(j %in% p)
          {
            y[j+2] = x[ind]
            ind = ind + 1 
          } else {
            y[j+2] = theta[j]
          }
        }
        y[6] = x[length(x)]
        if(y[6] < 1e-300) y[6] = 1e-300
        return(log(1/y[6])/(y[4]*(y[3]/y[4] - 1)) + 0.4/y[4])
      }
      calc.peak.intensity = function(x)
      {
        y = rep(NA,6); ind = 3
        y[1:2] = x[1:2]
        for(j in 1:3)
        {
          if(j %in% p)
          {
            y[j+2] = x[ind]
            ind = ind + 1 
          } else {
            y[j+2] = theta[j]
          }
        }
        return(1 - x[4]/x[3] + log(x[4]/x[3])/(x[3]/x[4]))
      }
      peak_intensity_calc = apply(rbind(mystate[1:2,],mystate[-(1:2),],out$state[1,,1]),2,calc.peak.intensity)
      peak_time_calc = apply(rbind(mystate[1:2,],mystate[-(1:2),],out$state[1,,1]),2,calc.peak.time)

      pf.peak.calc[[j]] = matrix(NA,nr=2,nc=3)
      pf.peak.calc[[j]][1,] = wtd.quantile(peak_time_calc[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
      pf.peak.calc[[j]][2,] = wtd.quantile(peak_intensity_calc[,1],out$weight[,days[j]+1],normwt=T,probs=c(.025,.5,.975))
    }
   
    # Save data
    save.image(paste(dpath,"PF-pred-",meth,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  }
}

# Perform prediction
require(plyr)
mydata = expand.grid(n = c(100,1000,10000), filt = c("BF","APF","KD"), resamp = c("multinomial","residual","stratified","systematic"), prior = c("normal","uniform"), meth = "both", stringsAsFactors = FALSE)
m_ply(mydata, pf.pred)

# Plot predicted % population infected, peak times/intensities
pred.plots(n, filt, resamp, prior, meth)
{
  # Load data
  load(paste(dpath,"PF-pred-",meth,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))

  # Create color depending on filt
  if(filt == "BF"){ col = 2
  } else if(filt == "APF"){ col = 4
  } else{ col = 3}
  
  # Construct plot of predicted % pop. infected
  pdf(paste(gpath,"PF-pred-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""),width=10,height=8)
  par(mfrow=c(2,3)) # dimensions of plot should depend on days
  ymax = max(sapply(pf.ui,max),max(sim$x[1,]))
  for(j in 1:length(days))
  {
    if(j == 1)
    {
      plot(1:nt,sim$x[1,-1],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab="% Population",main="Predicted %Pop Infected")
      mtext(paste("k = ",days[j]," Days",sep=""),side=3)
      abline(v=days[j]+1)
    } else {
      plot(1:nt,sim$x[1,-1],type="l",ylim=c(0,ymax),xlab="",ylab="")
      mtext(paste("k = ",days[j]," Days",sep=""),side=3)
      abline(v=days[j]+1)
    }
    lines((days[j]+1):nt,pf.i[[j]],col=2)
    lines((days[j]+1):nt,pf.li[[j]],lty=2,col=2)
    lines((days[j]+1):nt,pf.ui[[j]],lty=2,col=2)
    if(j == 1){ legend("topright",legend=c("Truth","Median","95% CI"),lty=c(1,1,2),col=c(1,2,2))}
  }
  dev.off()

  # Plot predicted peak times and intensities
  tpeak = which(sim$x[1,]==max(sim$x[1,]))-1
  ipeak = max(sim$x[1,])
  if(meth == "prop" |  meth == "both")
  {
    times = sapply(pf.peak.prop,function(x) x[1,2])
    times.l = sapply(pf.peak.prop,function(x) x[1,1])
    times.u = sapply(pf.peak.prop,function(x) x[1,3])
    ints = sapply(pf.peak.prop,function(x) x[2,2])
    ints.l = sapply(pf.peak.prop,function(x) x[2,1])
    ints.u = sapply(pf.peak.prop,function(x) x[2,3])
  }
  if(meth == "calc" | meth == "both")
  {
    timesc = sapply(pf.peak.calc,function(x) x[1,2])
    timesc.l = sapply(pf.peak.calc,function(x) x[1,1])
    timesc.u = sapply(pf.peak.calc,function(x) x[1,3])
    intsc = sapply(pf.peak.calc,function(x) x[2,2])
    intsc.l = sapply(pf.peak.calc,function(x) x[2,1])
    intsc.u = sapply(pf.peak.calc,function(x) x[2,3])
  }
  if(meth == "prop")
  {
    ymin = min(times.l,tpeak)
    ymax = max(times.u,tpeak)
    pdf(paste(gpath,"PF-peaks-",meth,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""),width=10,height=5)
    par(mfrow=c(1,2))
    plot(days,times,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),xlab="Days",ylab="Peak (Days)",main="Predicted Time of Epidemic Peak")
    segments(days,times.l,days,times.u)
    abline(h=tpeak,lty=2)
    legend("top",legend=c("Predicted","95% Cred Int","Truth"),lty=c(NA,1,2),pch=c(1,NA,NA))
    ymin = min(ints.l,ipeak)
    ymax = max(ints.u,ipeak)
    plot(days,ints,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),xlab="Days",ylab="Proportion Infected",main="Predicted Intensity of Epidemic Peak")
    segments(days,ints.l,days,ints.u)
    abline(h=ipeak,lty=2)
    dev.off()
  } else if(meth == "calc") {
    ymin = min(timesc.l,tpeak)
    ymax = max(timesc.u,tpeak)
    pdf(paste(gpath,"PF-peaks-",meth,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""),width=10,height=5)
    par(mfrow=c(1,2))
    plot(days,timesc,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),xlab="Days",ylab="Peak (Days)",main="Predicted Time of Epidemic Peak")
    segments(days,timesc.l,days,timesc.u)
    abline(h=tpeak,lty=2)
    legend("top",legend=c("Predicted","95% Cred Int","Truth"),lty=c(NA,1,2),pch=c(1,NA,NA))
    ymin = min(intsc.l,ipeak)
    ymax = max(intsc.u,ipeak)
    plot(days,intsc,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),xlab="Days",ylab="Proportion Infected",main="Predicted Intensity of Epidemic Peak")
    segments(days,intsc.l,days,intsc.u)
    abline(h=ipeak,lty=2)
    dev.off()
  } else { 
    ymin = min(times.l,timesc.l,tpeak)
    ymax = max(times.u,timesc.u,tpeak)
    ymin = min(times.l,tpeak)
    ymax = max(times.u,tpeak)
    pdf(paste(gpath,"PF-peaks-",meth,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""),width=10,height=5)
    par(mfrow=c(1,2))
    plot(days,times,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),xlab="Days",ylab="Peak (Days)",main="Predicted Time of Epidemic Peak")
    points(days+0.5,timesc,pch=2)
    segments(days,times.l,days,times.u)
    segments(days+.5,timesc.l,days+.5,timesc.u)
    abline(h=tpeak,lty=2)
    legend("top",legend=c("Prop","Calc","95% Cred Int","Truth"),lty=c(NA,NA,1,2),pch=c(1,2,NA,NA))
    ymin = min(ints.l,intsc.l,ipeak)
    ymax = max(ints.u,intsc.u,ipeak)
    plot(days,ints,ylim=c(ymin,ymax),xlim=c(min(days),max(days)+1),xlab="Days",ylab="Proportion Infected",main="Predicted Intensity of Epidemic Peak")
    points(days+0.5,intsc,pch=2)
    segments(days,ints.l,days,ints.u)
    segments(days+.5,intsc.l,days+.5,intsc.u)
    abline(h=ipeak,lty=2)
    dev.off()
  }
}

# Construct plots
m_ply(mydata, pred.plots)