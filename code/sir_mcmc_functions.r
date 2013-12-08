source("sir_functions.r")

# sir_mcmc - function to run MCMC for SIR model described in 'sir_functions.r'
# Arguments:
# y - L by T data matrix with L = number of data streams and T = number of time points. y may contain NA's for missing data.
# psi - list contain known parameters P, b, varsigma, sigma, and eta
# initial - optionally, a list specifying initial values of unknown states x (2 by T+1 matrix) parameters theta (length-3 vector)
# tuning - optionally, a list containing initial values of tuning parameters tuning.x (2 by T+1 matrix) for sampling states and tuning.theta (length-3 vector) for sampling parameters
# mcmc.details - optionally, a list containing scalars n.thin, n.sims, and n.burn for the how often to save iterations, total number of iterations, and burn in period, respectively, and boolean tune to specify whether tuning parameters should be adjusted during burn-in period
# steps - optionally, a character vector specifying which states/unknown parameters to sample from in the Gibbs sampler. May contain any of 'x', 'beta', 'gamma', or 'nu'.
# progress - boolean, should progress bar be displayed?
sir_mcmc <- function(y, psi, initial, tuning, mcmc.details, steps, progress=TRUE) {
  nt = dim(y)[2] # How many time points?

  # Deal with missing arguments
  if(missing(initial)) {
    # Set up initial values
    log.params = find.mu.sigma(c(.14, .09, .95), c(.5, .143, 1.3))
    tmp = rprior(function() exp(rnorm(3, log.params[[1]], log.params[[2]])))
    theta = tmp$theta
    x = matrix(NA, 2, nt + 1)
    x[,1] = tmp$x
    beta.x <- exp(rnorm(nt, log.params[[1]][1], log.params[[2]][1]))
    gamma.x <- exp(rnorm(nt, log.params[[1]][2], log.params[[2]][2]))
    nu.x <- exp(rnorm(nt, log.params[[1]][3], log.params[[2]][3]))
    for(i in 1:nt) x[,i+1] = revo(x[,i], c(beta.x[i], gamma.x[i], nu.x[i]), psi$P)
  } else {
    x = initial$x
    theta = initial$theta
  }
  if(missing(tuning))
  {
    tuning.x = matrix(0.001, nr=2, nc = nt+1)
    tuning.theta = c(0.01, 0.001, 0.01)
  } else {
    tuning.x = tuning$tuning.x
    tuning.theta = tuning$tuning.theta
  }
  if (missing(mcmc.details)) {
    n.thin <- 1   # save every n.thin(th) iteration
    n.sims <- 10
    n.burn <- 0
    tune = FALSE
  } else {
    n.thin = mcmc.details$n.thin
    n.sims = mcmc.details$n.sims
    n.burn = mcmc.details$n.burn
    tune = mcmc.details$tune
  }
  if (missing(steps)) {
    steps=c('x','beta','gamma','nu')
  } 

  # save structures
  n.iter <- (n.sims - n.burn)%/%n.thin
  keep.x   <- array(NA, c(n.iter, 2, nt + 1))
  keep.theta <- matrix(NA, n.iter, 3)
  accept.x <- matrix(NA, n.sims, nt + 1)
  accept.beta <- accept.gamma <- accept.nu <- rep(NA, n.sims)
  tune.states.all <- array(NA, c(2, nt+1, n.sims))
  tune.params.all <- matrix(NA, nr=n.sims, nc=3)

  # Run mcmc
  if(progress) pb = txtProgressBar(0,n.sims,style=3)
  for (i in 1:n.sims) {
    if (progress) setTxtProgressBar(pb,i)
    if(i <= n.burn & tune) tune.burn = TRUE else tune.burn = FALSE
    
    if('x'  %in% steps)
    {
      samp.x = sample.x(y, x, theta, psi, tuning.x, tune.burn)
      x = samp.x$x
      tuning.x <- tune.states.all[,,i] <- samp.x$tuning
      accept.x[i,] = samp.x$accept
    }
    if('beta' %in% steps)
    {
      samp.beta = sample.beta(y, x, theta, psi, tuning.theta[1], tune.burn)
      theta[1] = samp.beta$beta
      tuning.theta[1] <- tune.params.all[i,1] <- samp.beta$tuning
      accept.beta[i] = samp.beta$accept
    }
    if('gamma' %in% steps)
    {
      samp.gamma = sample.gamma(y, x, theta, psi, tuning.theta[2], tune.burn)
      theta[2] = samp.gamma$gamma
      tuning.theta[2] <- tune.params.all[i,2] <- samp.gamma$tuning
      accept.gamma[i] = samp.gamma$accept
    }
    if('nu' %in% steps)
    {
      samp.nu = sample.nu(y, x, theta, psi, tuning.theta[3], tune.burn)
      theta[3] = samp.nu$nu
      tuning.theta[3] <- tune.params.all[i,3] <- samp.nu$tuning
      accept.nu[i] = samp.nu$accept
    }

    # Only save every n.thin iteration
    if (ii <- save.iteration(i,n.burn,n.thin)) {
      keep.x[ii,,] = x
      keep.theta[ii,] = theta
    }
  }

  return(list(x=keep.x, theta=keep.theta, tuning=list(states=tune.states.all,params=tune.params.all), accept.x=accept.x, accept.theta=cbind(accept.beta,accept.gamma,accept.nu)))
}

#rmove <- function(y, x, theta, psi, tuning)
#{
#  nt = dim(y)[2]
#  stopifnot(nt == dim(x)[2] - 1)
#
#  # Move states
#  x = sample.x(y, x, theta, psi, tuning[1:2])
#  
#  # Move parameters
#  theta[1] = sample.beta(y, x, theta, psi, tuning[3])
#  theta[2] = sample.gamma(y, x, theta, psi, tuning[4])
#  theta[3] = sample.nu(y, x, theta, psi, tuning[5])
#  
#  return(list(state=x, theta=theta))
#}

# sample.x - function to sample from the full conditional of the state vector x at each time point separately, starting from initial state x[,0]
sample.x <- function(y, x, theta, psi, tuning, tune = FALSE)
{
  # Check dimensions
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  accept = rep(FALSE, nt + 1)
  
  # resample prior state
  i0 = rx(x[,1], tuning[,1])[2]
  x.prop = c(1-i0,i0)
  logMH = dlevo(x[,2], x.prop, theta, psi$P) + dlx0(x.prop) - dlevo(x[,2], x[,1], theta, psi$P) - dlx0(x[,1])
  if(log(runif(1)) < logMH & in.Omega(x.prop))
  {
    x[,1] = x.prop
    accept[1] = TRUE
    if(tune) tuning[,1] = tuning[,1]*1.1
  } else {
    if(tune) tuning[,1] = tuning[,1] / 1.1 
  }
  
  # resample middle states
  for(i in 2:nt)
  {
    x.prop = rx(x[,i], tuning[,i])
    logMH = dllik(y[,i-1], x.prop, psi$b, psi$varsigma, psi$sigma, psi$eta) + dlevo(x[,i+1], x.prop, theta, psi$P) + dlevo(x.prop, x[,i-1], theta, psi$P) - dllik(y[,i-1], x[,i], psi$b, psi$varsigma, psi$sigma, psi$eta) - dlevo(x[,i+1], x[,i], theta, psi$P) - dlevo(x[,i], x[,i-1], theta, psi$P)
    if(log(runif(1)) < logMH & in.Omega(x.prop))
    {
      x[,i] = x.prop
      accept[i] = TRUE
      if(tune) tuning[,i] = tuning[,i]*1.1
    } else {
      if(tune) tuning[,i] = tuning[,i] / 1.1 
    }
  }
  
  # resample last state
  x.prop = rx(x[,nt+1], tuning[,nt+1])
  logMH = dllik(y[,nt], x.prop, psi$b, psi$varsigma, psi$sigma, psi$eta) + dlevo(x.prop, x[,nt], theta, psi$P) - dllik(y[,nt], x[,nt+1], psi$b, psi$varsigma, psi$sigma, psi$eta) - dlevo(x[,nt+1], x[,nt], theta, psi$P)
  if(log(runif(1)) < logMH & in.Omega(x.prop))
  {
    x[,nt+1] = x.prop
    accept[nt + 1] = TRUE
    if(tune) tuning[,nt+1] = tuning[,nt+1]*1.1
  } else {
    if(tune) tuning[,nt+1] = tuning[,nt+1] / 1.1
  }
  return(list(x=x,accept=accept,tuning=tuning))
}

sample.beta <- function(y, x, theta, psi, tuning, tune = FALSE)
{
  # Check dimensions
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  accept = FALSE
  
  # Propose beta
  beta.curr = theta[1]
  beta.prop = rnorm(1, beta.curr, tuning)
  
  # Calculate MH ratio and iterate
  dlevo.prop = dlevo.curr = 0
  for(i in 1:nt)
  {
    dlevo.prop = dlevo.prop + dlevo(x[,i+1], x[,i], c(beta.prop,theta[2:3]), psi$P)
    dlevo.curr = dlevo.curr + dlevo(x[,i+1], x[,i], theta, psi$P)
  }
  logMH = dlevo.prop + dlbeta0(beta.prop) - dlevo.curr - dlbeta0(beta.curr)
  if(log(runif(1)) < logMH & beta.prop > 0)
  { 
    beta.curr = beta.prop
    accept = TRUE
    if(tune) tuning = tuning*1.1
  } else {
    if(tune) tuning = tuning / 1.1
  }
  
  return(list(beta=beta.curr,accept=accept,tuning=tuning))
}

sample.gamma <- function(y, x, theta, psi, tuning, tune = FALSE)
{
  # Check dimensions
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  accept = FALSE
  
  # Propose beta
  gamma.curr = theta[2]
  gamma.prop = rnorm(1, gamma.curr, tuning)
  
  # Calculate MH ratio and iterate
  dlevo.prop = dlevo.curr = 0
  for(i in 1:nt)
  {
    dlevo.prop = dlevo.prop + dlevo(x[,i+1], x[,i], c(theta[1],gamma.prop,theta[3]), psi$P)
    dlevo.curr = dlevo.curr + dlevo(x[,i+1], x[,i], theta, psi$P)
  }
  logMH = dlevo.prop + dlgamma0(gamma.prop) - dlevo.curr - dlgamma0(gamma.curr)
  if(log(runif(1)) < logMH & gamma.prop > 0)
  {
    gamma.curr = gamma.prop
    accept = TRUE
    if(tune) tuning = tuning*1.1
  } else {
    if(tune) tuning = tuning / 1.1
  }

  return(list(gamma=gamma.curr,accept=accept,tuning=tuning))
}

sample.nu <- function(y, x, theta, psi, tuning, tune = FALSE)
{
  # Check dimensions
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  accept = FALSE
  
  # Propose beta
  nu.curr = theta[3]
  nu.prop = rnorm(1, nu.curr, tuning)
  
  # Calculate MH ratio and iterate
  dlevo.prop = dlevo.curr = 0
  for(i in 1:nt)
  {
    dlevo.prop = dlevo.prop + dlevo(x[,i+1], x[,i], c(theta[1:2],nu.prop), psi$P)
    dlevo.curr = dlevo.curr + dlevo(x[,i+1], x[,i], theta, psi$P)
  }
  logMH = dlevo.prop + dlnu0(nu.prop) - dlevo.curr - dlnu0(nu.curr)
  if(log(runif(1)) < logMH & nu.prop > 0)
  {
    nu.curr = nu.prop
    accept = TRUE
    if(tune) tuning = tuning*1.1
  } else {
    if(tune) tuning = tuning / 1.1
  }
  
  return(list(nu=nu.curr,accept=accept,tuning=tuning))
}

save.iteration <- function(current.iter, n.burn, n.thin) {
  if (current.iter<n.burn | ((current.iter-n.burn)%%n.thin)!=0) {
    return(0)
  } else {
    return((current.iter-n.burn)/n.thin)
  }
}

#### Utility functions ####

# rx - samples a new state x from the proposal density of x given a tuning parameter
rx <- function(x, tuning = c(0.001, 0.001)) rnorm(2, x, tuning)

# in.Omega - returns TRUE if all elements of x are nonnegative and sum to less than or equal to 1
in.Omega <- function(x) all(x >= 0) & sum(x) <= 1
