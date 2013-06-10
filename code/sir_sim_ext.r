# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data from original model
load(paste(dpath,"sim-orig.rdata",sep=""))

# Create table of data for .csv file
epid.data = data.frame(seq(0,125,1),cbind(sim$x[2,],sim$x[1,],1-sim$x[2,]-sim$x[1,]),t(cbind(rep(NA,4),sim$y)))

# Alter parameters in observation equation
b = .25
varsigma = 1
sigma = 0.001
mu = 2

# Generate new observations
z = matrix(NA,nr=1,nc=nt)
for(i in 1:nt) z[1,i] = robs_sim(sim$x[1,i+1])

# Replace old observations with new ones
sim$y = z

# Save data
save.image(paste(dpath,"sim-ext.rdata",sep=""))
#load(paste(dpath,"sim-ext.rdata",sep=""))

# Add to .csv file of data and save
epid.data = data.frame(epid.data,c(NA,exp(sim$y[1,])))
names(epid.data) = c("Day","s","i","r","Stream 1","Stream 2","Stream 3","Stream 4","Stream 1 (Extended Analysis)")
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
pdf(paste(gpath,"sim-ext-z.pdf",sep=""))
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
legend("topright",legend=1,pch=1,col=1,title="Stream",cex=1.5)
dev.off()

# Clear objects
rm(list=ls(all=TRUE))