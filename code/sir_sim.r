source("sir_functions.r")
source("ss_sim.r")

# Initialize .Random.seed
set.seed(sample(1:1000,1))

# Set data and graphics path
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

# Simulate epidemic
revo_sim = function(x){ revo(x, P, d, theta)}
robs_sim = function(x){ robs(x, b, varsigma, sigma, dpower)}
rinit_sim = function(){ rinit(10/P)}
nt = 125
sim = ss.sim(nt, revo_sim, robs_sim, rinit_sim)
save.image(paste(dpath,"sim-xy.rdata",sep=""))
#load(paste(dpath,"sim-xy.rdata",sep=""))

# Plot the data
# Epidemic curves
pdf(paste(gpath,"sim-x.pdf",sep=""))
par(mar=c(5,7,4,1)+.1)
plot(0:nt,1 - sim$x[1,] - sim$x[2,],ylim=c(0,1),col=4,ylab="% Population",xlab="Time (days)",main="True Epidemic Curves",cex.lab=2,cex.main=2,cex.axis=1.6)
points(0:nt,sim$x[2,])
points(0:nt,sim$x[1,],col=2)
legend("topright",legend=c("Susceptible","Infected","Recovered"),pch=rep(1,3),col=c(1,2,4),cex=1.5)
dev.off()

# log scale (y)
no = dim(sim$y)[1]
pdf(paste(gpath,"sim-y-log.pdf",sep=""))
par(mar=c(5,5,4,1)+.1)
x = which(!is.na(sim$y[1,]))
y = sim$y[1,x]
plot(x,y,ylim=c(min(sim$y,na.rm=T),max(sim$y,na.rm=T)),xlim=c(0,nt),xlab="Time (Days)",ylab=expression(y),cex.lab=1.5)
if(no > 1)
{
  for(i in 2:no)
  {
     x = which(!is.na(sim$y[i,]))
     y = sim$y[i,x]
     points(x,y,col=i)
  }
}
dev.off()

# original scale (z)
pdf(paste(gpath,"sim-y.pdf",sep=""))
par(mar=c(5,7,4,1)+.1)
x = which(!is.na(sim$y[1,]))
y = sim$y[1,x]
plot(x,exp(y),ylim=c(min(exp(sim$y),na.rm=T),max(exp(sim$y),na.rm=T)),xlim=c(0,nt),xlab="Time (Days)",ylab=expression(z),main="Simulated Data",cex.lab=2,cex.main=2,cex.axis=1.6)
if(no > 1)
{
  for(i in 2:no)
  {
     x = which(!is.na(sim$y[i,]))
     y = sim$y[i,x]
     points(x,exp(y),col=i)
  }
}
legend("topright",legend=seq(1,no,1),pch=rep(1,no),col=seq(1,no,1),title="Syndrome")
dev.off()

# Clear objects
rm(list=ls(all=TRUE))