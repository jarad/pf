source("sir_pmcmc_functions.r")

# Run pfilter
np = 100
params = matrix(NA, nr=3, nc=np)
for(i in 1:np) params[,i] = rprior(i)
rownames(params) = c('beta','gamma','nu')
sir.pf = pfilter(sir.pomp, params, np, filter.mean = TRUE, save.states=TRUE, save.params=TRUE, verbose=TRUE)

# Plot 95% credible intervals for beta
lbeta = sapply(sir.pf$saved.params, function(x) quantile(x[1,],.025))
ubeta = sapply(sir.pf$saved.params, function(x) quantile(x[1,],.975))
ymin = min(lbeta, mysims[[1]]$true.params$theta[1])
ymax = max(ubeta, mysims[[1]]$true.params$theta[1])
plot(1:dim(y)[2], lbeta, ylim=c(ymin,ymax), type="l", col = 2)
lines(1:dim(y)[2], ubeta, col = 2)
abline(h=mysims[[1]]$true.params$theta[1], col = "gray70", lwd = 3)

# Liu and West filter
rprocess.kd <- function(xstart, times, params, ...) rprocess(xstart, times, params, P, log.par = TRUE, ...)
skeleton.kd <- function(x, t, params, ...) skeleton(x, t, params, P, log.par = TRUE, ...)
rprior.kd <- function(params, ...) rprior(params, log = TRUE, ...)
dprior.kd <- function(params, log=FALSE, ...) dprior(params, log=FALSE, log.par = TRUE, ...)
sir.pomp.kd <- pomp(data = y, times = 1:dim(y)[2], t0 = 0, P=P, b=b, varsigma=varsigma, sigma=sigma, eta = eta, rprocess=rprocess.kd, rmeasure=rmeasure, dmeasure=dmeasure, skeleton = skeleton.kd, skeleton.type = "map", initializer=initializer, rprior=rprior.kd, dprior=dprior.kd, parameter.transform = par.trans, parameter.inv.transform = par.inv.trans)

np = 100
par.init = rprior.kd(1)
#par.init = matrix(NA, nr=3, nc=np)
#for(i in 1:np) par.init[,i] = log(rprior(i))
sir.pf.kd = bsmc(sir.pomp.kd, params = par.init, Np = np, est = c('beta','gamma','nu'), smooth = 0.1, verbose=TRUE, transform = TRUE)

# Calculate 95% credible intervals for theta
(ltheta = apply(sir.pf.kd$post, 1, function(x) quantile(x, .025)))
(utheta = apply(sir.pf.kd$post, 1, function(x) quantile(x, .975)))
mysims[[1]]$true.params$theta