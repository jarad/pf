# Set data and graphics path
dpath = "../data/"

## Perform Kolmogorov-Smirnov tests

# Load "true" posterior distributions - KD, stratified at 20000 particles
load(paste(dpath,"PF-KD-normal-stratified-20000.rdata",sep=""))
y = pf.out$ftheta(pf.out$out$theta)
require(smcUtils)
tt = dim(y)[3]; nt = tt - 1
M = 1000
y1 = array(NA,dim=c(dim(y)[1],M,dim(y)[3]))
for(i in 1:tt)
{
  tmp = resample(pf.out$out$weight[,i], M, method = "stratified")
  y1[,,i] = y[,tmp$indices,i]
}

# Function to calculate p-values of ks tests
ks <- function(n = 10000, filt = "KD", resamp = "stratified", prior = "normal", nonunif, thresh)
{
  require(smcUtils) 

  # Load data
  load(paste(dpath,"PF-",filt,"-",prior,"-",resamp,"-",n,"-",nonunif,"-",100*thresh,".rdata",sep=""))
  x = pf.out$ftheta(pf.out$out$theta)
  x1 = array(NA,dim=c(dim(x)[1],M,dim(x)[3]))
  tt = dim(x)[3]; nt = tt - 1
  for(i in 1:tt)
  {
    tmp = resample(pf.out$out$weight[,i], M, method = resamp)
    x1[,,i] = x[,tmp$indices,i]
  }

  # Calculate p-values over time
  pvals = matrix(NA,nr=nt,nc=dim(x)[1])
  for(j in 1:dim(x)[1])
  {
    for(i in 2:tt)
    {
      pvals[i-1,j] = ks.test(x1[j,,i], y1[j,,i])$p.value
    }
  }

  print(paste(nonunif, thresh, sep=", "))
  return(pvals)
}

require(plyr)
n = 10000
mydata = expand.grid(nonunif = c("ess","cov","entropy"), thresh = seq(.05,.95,.05), stringsAsFactors = FALSE)
mypvals = maply(mydata, ks)
save(mypvals, file = paste(dpath,"PF-KD-normal-stratified-",n,"-ks.rdata",sep=""))
