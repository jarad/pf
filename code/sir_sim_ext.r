source("sir_functions.r")
source("ss_sim.r")

# Initialize .Random.seed
rnorm(1)

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data from original model
load(paste(dpath,"sim-orig.rdata",sep=""))
sim.orig = mysim
sim = list(x=sim.orig$sim$x)
theta = sim.orig$true.params$theta
P = sim.orig$true.params$P
nt = dim(sim$x)[2] - 1

# Set parameters in observation equation
b = .25
varsigma = 1
sigma = 0.001
eta = 2

# Generate new observations
robs_sim = function(x){ robs(x, b, varsigma, sigma, eta)}
z = matrix(NA,nr=1,nc=nt)
for(i in 1:nt) z[1,i] = robs_sim(sim$x[,i+1])

# Add new observations to list
sim$y = z

# Save data
mysim = list(sim=sim,true.params=list(theta=theta,b=b,sigma=sigma,varsigma=varsigma,eta=eta,P=P))
save(mysim,file=paste(dpath,"sim-ext.rdata",sep=""))

# Create table of data for .csv file
epid.data = data.frame(seq(0,125,1),cbind(sim$x[2,],sim$x[1,],1-sim$x[2,]-sim$x[1,],c(NA,sim$y[1,])))
names(epid.data) = c("Day","s","i","r","Stream 1")
write.csv(epid.data,file=paste(dpath,"simdata-ext.csv",sep=""),row.names=FALSE)

# Export epid.data as latex xtable
names(epid.data) = c("Day","$s$","$i$","$r$","$y_{1,t}$")
require(xtable)
caption = "Simulated epidemic and syndromic data"
label = "tab:data"
align = "|c|c|ccc|c|"
digits = c(0,0,rep(6,4))
print(xtable(epid.data,caption,label,align,digits),type="latex",file="../latex/simdata-ext.txt",include.rownames=FALSE)

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

# Clear objects
rm(list=ls(all=TRUE))