load("../Data/sir.pf.test-6P.rdata") #Use same simulated data

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

  return(x)
}

# Get 95% Credible intervals of predicted epidemic curve after running PFs for different numbers of days
days <- seq(10,35,5)
bf.li <- bf.ui <- list()
apf.li <- apf.ui <- list()
kd.li <- kd.ui <- list()
length(bf.li) = length(bf.ui) = length(days)
length(apf.li) = length(apf.ui) = length(days)
length(kd.li) = length(kd.ui) = length(days)
for(j in 1:length(days))
{
  mystate = out$state[,,days[j]+1]
  mystate2 = out2$state[,,days[j]+1]
  mybeta = theta2u(out3$theta[1,,days[j]+1],betal,betau)
  mygamma = theta2u(out3$theta[2,,days[j]+1],gammal,gammau)
  mynu = theta2u(out3$theta[3,,days[j]+1],nul,nuu)
  myb = theta2u(out3$theta[4,,days[j]+1],bl,bu)
  myvarsigma = theta2u(out3$theta[5,,days[j]+1],varsigmal,varsigmau)
  mysigma = theta2u(out3$theta[6,,days[j]+1],sigmal,sigmau)
  mystate3 = rbind(out3$state[,,days[j]+1],mybeta,mygamma,mynu,myb,myvarsigma,mysigma)
  pred = pred2 = pred3 = matrix(NA,n,nt-days[j])
  pb = txtProgressBar(0,n,style=3)
  for(i in 1:n)
  {
    setTxtProgressBar(pb,i)
    pred[i,] = ss.pred(nt-(days[j]+1),pstate,function() return(mystate[,i]))[1,]
    pred2[i,] = ss.pred(nt-(days[j]+1),pstate,function() return(mystate2[,i]))[1,]
    pred3[i,] = ss.pred(nt-(days[j]+1),pstate,function() return(mystate3[,i]))[1,]
  }
  bf.li[[j]] = bf.ui[[j]] = apf.li[[j]] = apf.ui[[j]] = kd.li[[j]] = kd.ui[[j]] = rep(NA,nt-days[j])
  for(i in 1:(nt-days[j]))
  {
    bf.li[[j]][i] = wtd.quantile(pred[,i], out$weight[,days[j]+1], normwt=T, probs=.025)
    bf.ui[[j]][i] = wtd.quantile(pred[,i], out$weight[,days[j]+1], normwt=T, probs=.975)
    apf.li[[j]][i] = wtd.quantile(pred2[,i], out2$weight[,days[j]+1], normwt=T, probs=.025)
    apf.ui[[j]][i] = wtd.quantile(pred2[,i], out2$weight[,days[j]+1], normwt=T, probs=.975)
    kd.li[[j]][i] = wtd.quantile(pred3[,i], out3$weight[,days[j]+1], normwt=T, probs=.025)
    kd.ui[[j]][i] = wtd.quantile(pred3[,i], out3$weight[,days[j]+1], normwt=T, probs=.975)
  }
  print(j)
}

# Plot true % population infected
# Bootstrap filter
plot(sim$x[1,],type="l",ylim=c(0,.43),xlab="Time (days)",ylab="% Population",
	main="95% Credible Intervals of Predicted %Pop Infected")
mtext(paste("n = ",n," particles",sep=""),side=3)
for(j in 1:length(days))
{
  lines((days[j]+1):nt,bf.li[[j]],lty=2,col=j)
  lines((days[j]+1):nt,bf.ui[[j]],lty=2,col=j)
}
legend("topright",legend=c("Truth",days),lty=c(1,rep(2,length(days))),col=c(1,1:length(days)),title="Days")

# Auxilliary particle filter
windows()
plot(sim$x[1,],type="l",ylim=c(0,.43),xlab="Time (days)",ylab="% Population",
	main="95% Credible Intervals of Predicted %Pop Infected")
mtext(paste("n = ",n," particles",sep=""),side=3)
for(j in 1:length(days))
{
  lines((days[j]+1):nt,apf.li[[j]],lty=2,col=j)
  lines((days[j]+1):nt,apf.ui[[j]],lty=2,col=j)
}
legend("topright",legend=c("Truth",days),lty=c(1,rep(2,length(days))),col=c(1,1:length(days)),title="Days")

# Kernel density particle filter
windows()
plot(sim$x[1,],type="l",ylim=c(0,.43),xlab="Time (days)",ylab="% Population",
	main="95% Credible Intervals of Predicted %Pop Infected")
mtext(paste("n = ",n," particles",sep=""),side=3)
for(j in 1:length(days))
{
  lines((days[j]+1):nt,kd.li[[j]],lty=2,col=j)
  lines((days[j]+1):nt,kd.ui[[j]],lty=2,col=j)
}
legend("topright",legend=c("Truth",days),lty=c(1,rep(2,length(days))),col=c(1,1:length(days)),title="Days")

# Histograms of unknown parameters at day 25
# bootstrap filter
require(plotrix)
cutoff = 26
msize = labsize = 1.5
pdf(paste("../Graphs/PF/Hist-BF-6P-",n,"-day",cutoff-1,".pdf",sep=""),width=10,height=8)
par(mfrow=c(2,3))
weighted.hist(out$state[3,,cutoff],out$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=beta)),side=3)
weighted.hist(out$state[4,,cutoff],out$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=gamma)),side=3)
weighted.hist(out$state[5,,cutoff],out$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=nu)),side=3)
weighted.hist(out$state[6,,cutoff],out$weight[,cutoff],
	xlab="b",main=expression(paste("Histogram of ",b,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste("b"," = ",aa,sep=""),list(aa=b)),side=3)
weighted.hist(out$state[7,,cutoff],out$weight[,cutoff],
	xlab=expression(varsigma),main=expression(paste("Histogram of ",varsigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=varsigma)),side=3)
weighted.hist(out$state[8,,cutoff],out$weight[,cutoff],
	xlab=expression(sigma),main=expression(paste("Histogram of ",sigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=sigma)),side=3)
dev.off()

# auxillary particle filter
pdf(paste("../Graphs/PF/Hist-APF-6P-",n,"-day",cutoff-1,".pdf",sep=""),width=10,height=8)
par(mfrow=c(2,3))
weighted.hist(out2$state[3,,cutoff],out2$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=beta)),side=3)
weighted.hist(out2$state[4,,cutoff],out2$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=gamma)),side=3)
weighted.hist(out2$state[5,,cutoff],out2$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=nu)),side=3)
weighted.hist(out2$state[6,,cutoff],out2$weight[,cutoff],
	xlab="b",main=expression(paste("Histogram of ",b,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste("b"," = ",aa,sep=""),list(aa=b)),side=3)
weighted.hist(out2$state[7,,cutoff],out2$weight[,cutoff],
	xlab=expression(varsigma),main=expression(paste("Histogram of ",varsigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=varsigma)),side=3)
weighted.hist(out2$state[8,,cutoff],out2$weight[,cutoff],
	xlab=expression(sigma),main=expression(paste("Histogram of ",sigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=sigma)),side=3)
dev.off()

# kernel density particle filter
pdf(paste("../Graphs/PF/Hist-KD-6P-",n,"-day",cutoff-1,".pdf",sep=""),width=10,height=8)
par(mfrow=c(2,3))
weighted.hist(theta2u(out3$theta[1,,cutoff],betal,betau),out3$weight[,cutoff],
	xlab=expression(beta),main="Histogram of Contact Rate",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(beta," = ",aa,sep=""),list(aa=beta)),side=3)
weighted.hist(theta2u(out3$theta[2,,cutoff],gammal,gammau),out3$weight[,cutoff],
	xlab=expression(gamma),main="Histogram of Recovery Time",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(gamma," = ",aa,sep=""),list(aa=gamma)),side=3)
weighted.hist(theta2u(out3$theta[3,,cutoff],nul,nuu),out3$weight[,cutoff],
	xlab=expression(nu),main="Histogram of Mixing Intensity",
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(nu," = ",aa,sep=""),list(aa=nu)),side=3)
weighted.hist(theta2u(out3$theta[4,,cutoff],bl,bu),out3$weight[,cutoff],
	xlab="b",main=expression(paste("Histogram of ",b,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste("b"," = ",aa,sep=""),list(aa=b)),side=3)
weighted.hist(theta2u(out3$theta[5,,cutoff],varsigmal,varsigmau),out3$weight[,cutoff],
	xlab=expression(varsigma),main=expression(paste("Histogram of ",varsigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(varsigma," = ",aa,sep=""),list(aa=varsigma)),side=3)
weighted.hist(theta2u(out3$theta[6,,cutoff],sigmal,sigmau),out3$weight[,cutoff],
	xlab=expression(sigma),main=expression(paste("Histogram of ",sigma,sep="")),
	cex.main=msize,cex.lab=labsize)
mtext(substitute(paste(sigma," = ",aa,sep=""),list(aa=sigma)),side=3)
dev.off()

