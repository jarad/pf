source("pf_functions.r")

# Set graphics and data path
gpath = "../graphs/"
dpath = "../data/"

# Initialize .Random.seed
set.seed(sample(1:1000,1))

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
pf.pred = function(n, filt, resamp, prior, cutoff)
{
  # Load data
  load(paste(dpath,"sim-xy.rdata",sep=""))
  load(paste(dpath,"PF-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta
  tt = nt + 1

  # Grab particles from filter at time cutoff
  out.theta = array(NA,dim=c(3,n,tt))
  if(!("theta" %in% names(out))){
    for(j in 1:3) out.theta[j,,] = ftheta(out$state[j+2,,],j)
  } else {
    for(j in 1:3) out.theta[j,,] = ftheta(out$theta[j,,],j)
  }
  mystate = rbind(out$state[1:2,,cutoff+1],out.theta[,,cutoff+1])
  pstate_revo = function(x) revo(x, P, d, x[3:5], FALSE)

  # Propagate particles forward to end of epidemic
  pred = array(NA,dim=c(2,n,tt-cutoff))
  pb = txtProgressBar(0,n,style=3)
  for(i in 1:n)
  {
    setTxtProgressBar(pb,i)
    pred[,i,] = ss.pred(nt-cutoff,pstate_revo,function() return(mystate[,i]))[1:2,]
  }

  # Save predicted particles with associated weights
  wts = matrix(NA,nr=n,nc=tt-cutoff)
  for(i in 1:(tt-cutoff)) wts[,i] = out$weight[,cutoff+1]
  pred = list(pred=pred,wts=wts)
  save(pred,file=paste(dpath,"PF-pred-",cutoff,"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
  return(pred)
}

# Perform prediction for different cutoff day values
require(plyr)
days = seq(10,35,5)
mydata = expand.grid(n = 10000, filt = "KD", resamp = "stratified", prior = "normal", cutoff = days, stringsAsFactors = FALSE)
pp = mlply(mydata,pf.pred)
#pp = list(); length(pp) = length(days)
#for(i in 1:length(days))
#{ 
#  load(paste(dpath,"PF-pred-",days[i],"-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
#  pp[[i]] = pred
#}

# Calculate .025 and .975 quantiles of predicted particle samples
probs = c(0.025,0.975)
quant = list(); length(quant) = length(days)
for(i in 1:length(days)) quant[[i]] = pf.quantile(pp[[i]]$pred, pp[[i]]$wts, function(x,param=1) x, probs)
save(quant,file=paste(dpath,"PF-pred-quant-KD-normal-stratified-10000.rdata",sep=""))

# Construct plot of predicted % pop. infected
load(paste(dpath,"sim-xy.rdata",sep=""))
ymax = max(sapply(quant,function(x) max(x[,1,2])),max(sim$x[1,]))
msize = 3
labsize = 2.5
axsize = 1.5
pdf(file=paste(gpath,"PF-pred-KD-normal-stratified-10000.pdf",sep=""),width=10,height=8)
par(mfrow=c(2,3),mar=c(6,8,2,1)+.1,mgp=c(4,1,0)) # dimensions of plot should depend on days
for(j in 1:length(days))
{
  if(j == 1)
  {
    plot(1:nt,sim$x[1,-1],type="l",ylim=c(0,ymax),xlab="Time (days)",ylab="% Population Infected",main=paste("t = ",days[j]," Days",sep=""),cex.main=msize,cex.lab=labsize,cex.axis=axsize)
    abline(v=days[j],lty=2)
  } else {
    plot(1:nt,sim$x[1,-1],type="l",ylim=c(0,ymax),xlab="",ylab="",main=paste("t = ",days[j]," Days",sep=""),cex.main=msize,cex.axis=axsize)
    abline(v=days[j],lty=2)
  }
  lines((days[j]):nt,quant[[j]][,1,1],col=2)
  lines((days[j]):nt,quant[[j]][,1,2],col=2)
  if(j == 1){ legend("topright",legend=c("Truth","95% CI"),lty=c(1,1),col=c(1,2),cex=1.5)}
}
dev.off()
