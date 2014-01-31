source("sir_functions.r")
source("ss_sim.r")

# Set .Random.seed
set.seed(60)

# Set data and graphics path
gpath = "../graphs/"
dpath = "../data/"

# Set known parameter values
P = 5000
b = c(.25, .27, .23, .29)
varsigma = c(1.07, 1.05, 1.01, .98)
sigma = c(.0012, .0008, .0010, .0011)
eta = rep(0,4)

nt = 125
n.sims = 40
sims = list()
mysims = list()
for(i in 1:n.sims)
{
  # Set unknown parameter values and initial proportion of population infected
  rtheta <- function()
  {
    theta <- rep(NA, 3)
    log.params <- find.mu.sigma(c(1.5, .09, .95), c(3, .143, 1.3))
    theta[2:3] <- exp(rnorm(2, log.params[[1]][2:3], log.params[[2]][2:3]))
    theta[1] <- theta[2]*exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
    return(theta)
  }
  init <- rprior(rtheta)
  
  # Simulate epidemic
  revo_sim = function(x){ revo(x, init$theta, P)}
  robs_sim = function(x){ robs(x, b, varsigma, sigma, eta)}
  rinit_sim = function(){ return(init$x)}
  sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
  mysims[[i]] = list(sim=sim,true.params=list(theta=init$theta,b=b,sigma=sigma,varsigma=varsigma,eta=eta,P=P))
}
save(mysims,file=paste(dpath,"sim-orig.rdata",sep=""))

# Write data to .csv files
for(i in 1:n.sims)
{
  epid.data = data.frame(seq(0,125,1),cbind(mysims[[i]]$sim$x[1,],mysims[[i]]$sim$x[2,],1-mysims[[i]]$sim$x[2,]-mysims[[i]]$sim$x[1,]),t(cbind(rep(NA,dim(mysims[[i]]$sim$y)[1]),mysims[[i]]$sim$y)))
  streams = rep(NA, dim(mysims[[i]]$sim$y)[1])
  for(j in 1:dim(mysims[[i]]$sim$y)[1]) streams[j] = paste("Stream ",j,sep="")
  names(epid.data) = c("Day","s","i","r",streams)
  write.csv(epid.data,file=paste(dpath,"simdata-orig-",i,".csv",sep=""),row.names=FALSE)
}

# Plot the data
# Epidemic curves
pdf(paste(gpath,"sim-orig-epid.pdf",sep=""))
par(mar=c(5,7,4,1)+.1)
if(n.sims > 1)
{
  for(i in 2:n.sims)
  {
    if(i == 2)
    {
      plot(0:nt,mysims[[i]]$sim$x[2,],type="l",ylim=c(0,1),col="gray80",ylab="% Population",xlab="Time (days)",main="True Epidemic Curves",cex.lab=2,cex.main=2,cex.axis=1.6)
#    lines(0:nt,1 - mysims[[i]]$sim$x[1,] - mysims[[i]]$sim$x[2,],col=4)
#    lines(0:nt,mysims[[i]]$sim$x[1,])
    } else {
      lines(0:nt,mysims[[i]]$sim$x[2,],col="gray80")
    }
  }
  lines(0:nt,1 - mysims[[1]]$sim$x[1,] - mysims[[1]]$sim$x[2,],col=4)
  lines(0:nt,mysims[[1]]$sim$x[1,])
  lines(0:nt,mysims[[1]]$sim$x[2,],col=2)
} else {
  plot(0:nt,1 - mysims[[1]]$sim$x[1,] - mysims[[1]]$sim$x[2,],type="l",ylim=c(0,1),col=4,ylab="% Population",xlab="Time (days)",main="True Epidemic Curves",cex.lab=2,cex.main=2,cex.axis=1.6)
  lines(0:nt,mysims[[1]]$sim$x[1,])
  lines(0:nt,mysims[[1]]$sim$x[2,],col=2)
}
legend("topright",legend=c("Susceptible","Infected","Recovered"),lty=rep(1,3),col=c(1,2,4),cex=1.5)
dev.off()

# Epidemic curves for each simulation
for(j in 1:n.sims)
{
  pdf(paste(gpath,"sim-orig-epid-",j,".pdf",sep=""))
  par(mar=c(5,7,4,1)+.1)
  plot(0:nt,1 - mysims[[j]]$sim$x[1,] - mysims[[j]]$sim$x[2,],type="l",ylim=c(0,1),col=4,ylab="% Population",xlab="Time (days)",main="True Epidemic Curves",cex.lab=2,cex.main=2,cex.axis=1.6)
  lines(0:nt,mysims[[j]]$sim$x[1,])
  lines(0:nt,mysims[[j]]$sim$x[2,],col=2)
  legend("topright",legend=c("Susceptible","Infected","Recovered"),lty=rep(1,3),col=c(1,2,4),cex=1.5)
  dev.off()
}

# Observed data (z) over all simulations
pdf(paste(gpath,"sim-orig-z.pdf",sep=""))
par(mar=c(5,5,4,1)+.1)
miny = min(sapply(mysims, function(x) min(x$sim$y,na.rm=T)))
maxy = max(sapply(mysims, function(x) max(x$sim$y,na.rm=T)))
no = dim(mysims[[1]]$sim$y)[1]
x = which(!is.na(mysims[[1]]$sim$y[1,]))
z = mysims[[1]]$sim$y[1,x]
plot(x,z,ylim=c(miny,maxy),xlim=c(0,nt),type="l",xlab="Time (Days)",ylab=expression(paste("Observed data (",y[t],")",sep="")),main="Syndromic Data",cex.lab=2,cex.main=2,cex.axis=1.6)
if(no > 1)
{
  for(i in 2:no)
  {
     x = which(!is.na(mysims[[1]]$sim$y[i,]))
     z = mysims[[1]]$sim$y[i,x]
     lines(x,z,col=i)
  }
}
if(n.sims > 1)
{
  for(j in 2:n.sims)
  {
    no = dim(mysims[[j]]$sim$y)[1]
    for(i in 1:no)
    {
      x = which(!is.na(mysims[[j]]$sim$y[i,]))
      z = mysims[[j]]$sim$y[i,x]
      lines(x,z,col=i)
    }
  }
}
legend("topright",legend=seq(1,length(streams),1),pch=rep(1,length(streams)),col=seq(1,length(streams),1),title="Stream",cex=1.5)
dev.off()

# Observed data for each individual simulation
for(j in 1:n.sims)
{
  no = dim(mysims[[j]]$sim$y)[1]
  pdf(paste(gpath,"sim-orig-z-",j,".pdf",sep=""))
  par(mar=c(5,5,4,1)+.1)
  x = which(!is.na(mysims[[j]]$sim$y[1,]))
  z = mysims[[j]]$sim$y[1,x]
  plot(x,z,ylim=c(miny,maxy),xlim=c(0,nt),xlab="Time (Days)",ylab=expression(paste("Observed data (",y[t],")",sep="")),main="Syndromic Data",cex.lab=2,cex.main=2,cex.axis=1.6)
  if(no > 1)
  {
    for(i in 2:no)
    {
       x = which(!is.na(mysims[[j]]$sim$y[i,]))
       z = mysims[[j]]$sim$y[i,x]
       points(x,z,col=i)
    }
  }
  legend("topright",legend=seq(1,length(streams),1),pch=rep(1,length(streams)),col=seq(1,length(streams),1),title="Stream",cex=1.5)
  dev.off()
}