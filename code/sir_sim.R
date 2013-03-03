source("sir_functions.R")
source("ss_sim.R")

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

# Plot the data
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
pdf(paste(gpath,"sim-y.pdf",sep=""))
par(mar=c(5,5,4,1)+.1)
x = which(!is.na(sim$y[1,]))
y = sim$y[1,x]
plot(x,exp(y),ylim=c(min(exp(sim$y),na.rm=T),max(exp(sim$y),na.rm=T)),xlim=c(0,nt),xlab="Time (Days)",ylab=expression(e^y),cex.lab=1.5)
if(no > 1)
{
  for(i in 2:no)
  {
     x = which(!is.na(sim$y[i,]))
     y = sim$y[i,x]
     points(x,exp(y),col=i)
  }
}
dev.off()

# Clear objects
rm(list=ls(all=TRUE))