# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data from original model
load(paste(dpath,"sim-xy.rdata",sep=""))

# Create table of data for .csv file
epid.data = data.frame(seq(0,125,1),cbind(sim$x[2,],sim$x[1,],1-sim$x[2,]-sim$x[1,]),t(cbind(rep(NA,4),exp(sim$y))))

# Alter parameters in observation equation
b = .25
varsigma = 1
sigma = 0.001

# Generate new observations
y = matrix(NA,nr=1,nc=nt)
for(i in 1:nt) y[1,i] = robs_sim(sim$x[1,i+1])

# Replace old observations with new ones
sim$y = y

# Save data
save.image(paste(dpath,"sim-xy-ext.rdata",sep=""))
#load(paste(dpath,"sim-xy-ext.rdata",sep=""))

# Add to .csv file of data and save
epid.data = data.frame(epid.data,c(NA,exp(sim$y[1,])))
names(epid.data) = c("Day","s","i","r","Syndrome 1","Syndrome 2","Syndrome 3","Syndrome 4","Syndrome 1 (Extended Analysis)")
write.csv(epid.data,file=paste(dpath,"simdata.csv",sep=""),row.names=FALSE)

# Export epid.data as latex xtable
epid.data = epid.data[,-1]
names(epid.data) = c("$s$","$i$","$r$","$z_{1,t}$","$z_{1,t}$","$z_{1,t}$","$z_{1,t}$","$z_{1,t}$ (Extended Analysis)")
rownames(epid.data) = seq(0,125,1)
require(xtable)
caption = "Simulated epidemic and syndromic data"
label = "tab:data"
align = "|c|ccc|cccc|c|"
digits = 6
print(xtable(epid.data,caption,label,align,digits),type="latex",file="../latex/simdata.txt")

# Plot the data
no = dim(sim$y)[1]
pdf(paste(gpath,"sim-y-log-ext.pdf",sep=""))
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
pdf(paste(gpath,"sim-y-ext.pdf",sep=""))
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