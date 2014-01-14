source("sir_functions.r")
source("ss_sim.r")

# Set .Random.seed
set.seed(60)

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data from original model
load(paste(dpath,"sim-orig.rdata",sep=""))
require(plyr)
sim = list()
length(sim) = 2
names(sim) = c("x","y")
sim$x <- mysims[[1]]$sim$x
theta = mysims[[1]]$true.params$theta
P = mysims[[1]]$true.params$P
nt = dim(sim$x)[2] - 1

# Set parameters in observation equation
b = .25
varsigma = 1
sigma = 0.001
eta = 2

# Generate new observations
robs_sim = function(x){ robs(x, b, varsigma, sigma, eta)}
sim$y = matrix(NA,nr=1,nc=nt)
for(i in 1:nt) sim$y[1,i] = robs_sim(sim$x[,i+1])

# Save data
mysim = list(sim=sim,true.params=list(theta=theta,b=b,sigma=sigma,varsigma=varsigma,eta=eta,P=P))
save(mysim,file=paste(dpath,"sim-ext.rdata",sep=""))

# Create table of data for .csv file
epid.data = data.frame(seq(0,125,1),cbind(sim$x[2,],sim$x[1,],1-sim$x[2,]-sim$x[1,],c(NA,sim$y[1,])))
names(epid.data) = c("Day","s","i","r","Stream 1")
write.csv(epid.data,file=paste(dpath,"simdata-ext.csv",sep=""),row.names=FALSE)

# Plot the data
no = dim(sim$y)[1]
pdf(paste(gpath,"sim-ext-z.pdf",sep=""))
par(mar=c(5,5,4,1)+.1)
x = which(!is.na(sim$y[1,]))
z = sim$y[1,x]
plot(x,z,ylim=c(min(sim$y,na.rm=T),max(sim$y,na.rm=T)),xlim=c(0,nt),xlab="Time (Days)",ylab=expression(paste("Observed data (",y[t],")",sep="")),main="Syndromic Data",cex.lab=2,cex.main=2,cex.axis=1.6)
if(no > 1)
{
  for(i in 2:no)
  {
     x = which(!is.na(sim$y[i,]))
     z = sim$y[i,x]
     points(x,z,col=i)
  }
}
legend("topright",legend=1,pch=1,col=1,title="Stream",cex=1.5)
dev.off()