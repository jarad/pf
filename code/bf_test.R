library(smcUtils)
demo(kalman)


# Functions written so they work in all 3 particle filters
revo = function(x, theta=log(sqrt(W))) rnorm(1,x,exp(theta))
dllik = function(y, x, theta=NULL) dnorm(y,x,sqrt(V),log=T)
rprior = function(.stateonly=TRUE) 
{
  if (.stateonly) 
  {
    return(rnorm(1, m0, sqrt(C0)))
  } else {
    list(x    = rnorm(1, m0, sqrt(C0)),   
         theta= log(rgamma(1,100*sqrt(W),100)))
  }
}
pstate = function(x, theta=NULL) x

# Run bootstrap filter 
# for unknown state only
source("bf.r")
out = bf(y, dllik, revo, rprior, J, 
         method="stratified",nonuniformity="ess",threshold=0.8,log=F)

# Run auxiliary particle filter 
# for unknown state only
source("apf.r")
out2 = apf(y, dllik, pstate, revo, rprior, J)

# Run Liu-West kernel density particle filter
# for unknown state and parameter
source("kd_pf.r")
rprior_kd = function() rprior(FALSE)
out3 = kd_pf(y, dllik, pstate, revo, rprior_kd, J)


# Comparison
tt = N+1
pbs = c(.025,.975)
qnorm(pbs,y[tt-1])

m[N]+c(-2,2)*sqrt(M[N])       # Truth via Kalman filter
pf.m[N]+c(-2,2)*sqrt(pf.v[N]) # Truth via demo bootstrap filter
wtd.quantile(x[N,], ws[N,], normwt=T, probs=pbs)
wtd.quantile(out$state[,,tt], out$weight[,tt], normwt=T, probs=pbs)
wtd.quantile(out2$state[,,tt], out2$weight[,tt], normwt=T, probs=pbs)
wtd.quantile(out3$state[,,tt], out3$weight[,tt], normwt=T, probs=pbs)
