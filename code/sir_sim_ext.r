source("sir_functions.r")
source("ss_sim.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data from original model
load(paste(dpath,"sim-orig.rdata",sep=""))
require(plyr)
n.sims = length(mysims)
nt = dim(mysims[[1]]$sim$y)[2]

# Set .Random.seed
set.seed(60)

for(i in 1:n.sims)
{
  # Set parameters in observation equation
  log.params = find.mu.sigma(c(.1, .85, .0005), c(.4, 1.15, .0015))
  b = exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
  varsigma = exp(rnorm(1, log.params[[1]][2], log.params[[2]][2]))
  sigma = exp(rnorm(1, log.params[[1]][3], log.params[[2]][3]))
  eta = rnorm(1, 2.5, 1)
  mysims[[i]]$true.params$b = b
  mysims[[i]]$true.params$varsigma = varsigma
  mysims[[i]]$true.params$sigma = sigma
  mysims[[i]]$true.params$eta = eta
  
  # Generate new observations
  robs_sim = function(x){ robs(x, b, varsigma, sigma, eta)}
  mysims[[i]]$sim$y = matrix(NA,nr=1,nc=nt)
  for(j in 1:nt) mysims[[i]]$sim$y[1,j] = robs_sim(mysims[[i]]$sim$x[,j+1])
}

# Save data
save(mysims,file=paste(dpath,"sim-ext.rdata",sep=""))

# Write data to .csv files
for(i in 1:n.sims)
{
  epid.data = data.frame(seq(0,125,1),cbind(mysims[[i]]$sim$x[1,],mysims[[i]]$sim$x[2,],1-mysims[[i]]$sim$x[2,]-mysims[[i]]$sim$x[1,]),t(cbind(rep(NA,dim(mysims[[i]]$sim$y)[1]),mysims[[i]]$sim$y)))
  streams = rep(NA, dim(mysims[[i]]$sim$y)[1])
  for(j in 1:dim(mysims[[i]]$sim$y)[1]) streams[j] = paste("Stream ",j,sep="")
  names(epid.data) = c("Day","s","i","r",streams)
  write.csv(epid.data,file=paste(dpath,"simdata-ext-",i,".csv",sep=""),row.names=FALSE)
}

# Observed data (z) over all simulations
pdf(paste(gpath,"sim-ext-z.pdf",sep=""))
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
  pdf(paste(gpath,"sim-ext-z-",j,".pdf",sep=""))
  par(mar=c(5,5,4,1)+.1)
  x = which(!is.na(mysims[[j]]$sim$y[1,]))
  z = mysims[[j]]$sim$y[1,x]
  miny = min(mysims[[j]]$sim$y, na.rm=T)
  maxy = max(mysims[[j]]$sim$y, na.rm=T)  
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