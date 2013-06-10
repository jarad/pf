source("sir_functions.r")
source("ss_sim.r")

# Initialize .Random.seed
rnorm(1)

# Set data and graphics path
gpath = "../graphs/"
dpath = "../data/"

# Set known parameter values
P = 5000
d = 1
b = c(.25, .27, .23, .29)
varsigma = c(1.07, 1.05, 1.01, .98)
sigma = c(.0012, .0008, .0010, .0011)
mu = rep(0,4)

# Set unknown parameter values
theta = c(0.2399, 0.1066, 1.2042)

# Simulate epidemic
revo_sim = function(x){ revo(x, P, d, theta)}
robs_sim = function(x){ robs(x, b, varsigma, sigma, mu)}
rinit_sim = function(){ rinit(10/P)}
nt = 125
sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
save.image(paste(dpath,"sim-orig.rdata",sep=""))
#load(paste(dpath,"sim-orig.rdata",sep=""))

# Plot the data
# Epidemic curves
pdf(paste(gpath,"sim-orig-epid.pdf",sep=""))
par(mar=c(5,7,4,1)+.1)
plot(0:nt,1 - sim$x[1,] - sim$x[2,],type="l",ylim=c(0,1),col=4,ylab="% Population",xlab="Time (days)",main="True Epidemic Curves",cex.lab=2,cex.main=2,cex.axis=1.6)
lines(0:nt,sim$x[2,])
lines(0:nt,sim$x[1,],col=2)
legend("topright",legend=c("Susceptible","Infected","Recovered"),lty=rep(1,3),col=c(1,2,4),cex=1.5)
dev.off()

# Observed data (z)
no = dim(sim$y)[1]
pdf(paste(gpath,"sim-orig-z.pdf",sep=""))
par(mar=c(5,5,4,1)+.1)
x = which(!is.na(sim$y[1,]))
z = sim$y[1,x]
plot(x,z,ylim=c(min(sim$y,na.rm=T),max(sim$y,na.rm=T)),xlim=c(0,nt),xlab="Time (Days)",ylab=expression(paste("Observed data (",z,")",sep="")),main="Syndromic Data",cex.lab=2,cex.main=2,cex.axis=1.6)
if(no > 1)
{
  for(i in 2:no)
  {
     x = which(!is.na(sim$y[i,]))
     z = sim$y[i,x]
     points(x,z,col=i)
  }
}
legend("topright",legend=seq(1,4,1),pch=rep(1,4),col=seq(1,4,1),title="Stream",cex=1.5)
dev.off()

# Clear objects
rm(list=ls(all=TRUE))