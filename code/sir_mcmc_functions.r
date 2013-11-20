source("sir_functions.r")

sir_mcmc <- function(y, psi, mcmc.details, steps, progress=TRUE) {
  nt = dim(y)[2]
  L = dim(y)[1]
  
  # Set up initial values
  tmp = rprior(function() exp(rnorm(2, c(-1.3296, -2.1764), c(.3248, .1183))))
  theta = tmp$theta
  x = matrix(NA, 2, nt + 1)
  x[,1] = tmp$x
  beta.x <- exp(rnorm(nt, -1.3296, .3248))
  gamma.x <- exp(rnorm(nt, -2.1764, .1183))
  for(i in 1:nt) x[,i+1] = revo(x[,i], c(beta.x[i], gamma.x[i], 1), psi$P)

  # Deal with missing arguments
  if (missing(mcmc.details)) {
    n.thin <- 1   # save every n.thin(th) iteration
    n.sims <- 10
    n.burn <- 0  
  } else {
    attach(mcmc.details)
  }
  if (missing(steps)) {
    steps=c('x','theta')
  } 

  # save structures
  n.iter <- n.burn+n.sims%/%n.thin
  keep.x   <- array(NA, c(n.sims, 2, nt + 1))
  keep.theta    <- matrix(NA, n.sims, 2)

  # Run mcmc
  if(progress) pb = txtProgressBar(0,nt,style=3)
  for (i in 1:n.sims) {
    if (progress) setTxtProgressBar(pb,i)
    
    if ('x'  %in% steps) x = sample.x(y, x, c(theta,1), psi)
    if ('theta'   %in%steps) theta = sample.theta(y, x, c(theta,1), psi)[1:2]

    # Only save every n.thin iteration
    if (ii <- save.iteration(i,n.burn,n.thin)) {
      keep.x[ii,,] = x
      keep.theta[ii,] = theta
    }
  }

  return(list(x=keep.x, theta=keep.theta))
}

rmove <- function(y, x, theta, psi)
{
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)

  # Move states
  states = sample.x(y, x, theta, psi)
  
  # Move parameters
  params = sample.theta(y, x, theta, psi)
  
  return(list(state=states, theta=params[1:2]))
}

sample.x <- function(y, x, theta, psi)
{
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)  
  x.prop = rx(x, theta, psi)
  logMH = dljoint(y, x.prop, theta, psi) + dlx(x, x.prop, theta, psi) - dljoint(y, x, theta, psi) - dlx(x.prop, x, theta, psi)
  if(log(runif(1)) < logMH) x = x.prop
  return(x)
}

sample.theta <- function(y, x, theta, psi, tuning = c(0.003, 0.001))
{
  # Sample beta
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  beta.prop = rnorm(1, theta[1], tuning[1])
  while(!(beta.prop > 0)) beta.prop = rnorm(1, theta[1], tuning[1])
  logMH = dljoint(y, x, c(beta.prop, theta[2:3]), psi) + dnorm(theta[1], beta.prop, tuning[1], log=TRUE) - dljoint(y, x, theta, psi) -  dnorm(beta.prop, theta[1], tuning[1], log=TRUE)
  if(log(runif(1)) < logMH) theta[1] = beta.prop
  
  # Sample gamma
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  gamma.prop = rnorm(1, theta[2], tuning[2])
  while(!(gamma.prop > 0)) gamma.prop = rnorm(1, theta[2], tuning[2])
  logMH = dljoint(y, x, c(theta[1], gamma.prop, theta[3]), psi) + dnorm(theta[2], gamma.prop, tuning[2], log=TRUE) - dljoint(y, x, theta, psi) -  dnorm(gamma.prop, theta[2], tuning[2], log=TRUE)
  if(log(runif(1)) < logMH) theta[2] = gamma.prop
  
  return(theta)
}

save.iteration <- function(current.iter, n.burn, n.thin) {
  if (current.iter<n.burn | ((current.iter-n.burn)%%n.thin)!=0) {
    return(0)
  } else {
    return((current.iter-n.burn)/n.thin)
  }
}

#### Utility functions ####

rx <- function(x, theta, psi, sdfac = 1)
{
  s <- rnorm(dim(x)[2], x[1,], sdfac*sqrt(theta[1] / psi$P^2))
  i <- rnorm(dim(x)[2], x[2,] + -1*(s - x[1,]), sdfac*sqrt(theta[2] / psi$P^2))
  while(!all(apply(cbind(s,i), 1, in.Omega)))
  {
    s <- rnorm(dim(x)[2], x[1,], sdfac*sqrt(theta[1] / psi$P^2))
    i <- rnorm(dim(x)[2], x[2,] + -1*(s - x[1,]), sdfac*sqrt(theta[2] / psi$P^2))
  }
  return(rbind(s,i))
}

dlx <- function(x, x.curr, theta, psi, sdfac=1)
{
  nt = dim(x)[2] - 1
  log.density = sum(dnorm(x[,1], c(x.curr[1,1], x.curr[2,1] + -1*(x[1,1] - x.curr[1,1])), sdfac*sqrt(theta[1:2] / psi$P^2), log=TRUE))
  for(i in 1:nt) log.density = log.density + sum(dnorm(x[,i+1], c(x.curr[1,i+1], x.curr[2,i+1] + -1*(x[1,i+1] - x.curr[1,i+1])), sdfac*sqrt(theta[1:2] / psi$P^2), log=TRUE))
  return(log.density)
}

dljoint <- function(y, x, theta, psi)
{
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  log.density = dlprior(x[,1], theta)
  for(i in 1:nt) log.density = log.density + dllik(y[,i], x[,i+1], psi$b, psi$varsigma, psi$sigma, psi$eta) + dlevo(x[,i+1], x[,i], theta, 
  psi$P)
  return(log.density)
}

in.Omega <- function(x) all(x >= 0) & sum(x) <= 1
