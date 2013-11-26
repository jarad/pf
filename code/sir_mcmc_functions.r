source("sir_functions.r")

sir_mcmc <- function(y, psi, initial, tuning, mcmc.details, steps, progress=TRUE) {
  nt = dim(y)[2]

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
    attach(initial)
  }
  if(missing(tuning))
  {
    tuning.x = c(0.001, 0.001)
    tuning.beta = 0.01
    tuning.gamma = 0.001
    tuning.nu = 0.01
  } else {
    attach(tuning)
  }
  if (missing(mcmc.details)) {
    n.thin <- 1   # save every n.thin(th) iteration
    n.sims <- 10
    n.burn <- 0
    tune = FALSE
  } else {
    attach(mcmc.details)
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
  tune.all <- matrix(NA, n.sims, 5)

  # Run mcmc
  if(progress) pb = txtProgressBar(0,n.sims,style=3)
  for (i in 1:n.sims) {
    if (progress) setTxtProgressBar(pb,i)
    if(i <= n.burn & tune) burn = TRUE else burn = FALSE
    
    if('x'  %in% steps)
    {
      samp.x = sample.x(y, x, theta, psi, tuning.x, burn)
      x = samp.x$x
      tuning.x <- tune.all[i,1:2] <- samp.x$tuning
      accept.x[i,] = samp.x$accept
    }
    if('beta' %in% steps)
    {
      samp.beta = sample.beta(y, x, theta, psi, tuning.beta, burn)
      theta[1] = samp.beta$beta
      tuning.beta <- tune.all[i,3] <- samp.beta$tuning
      accept.beta[i] = samp.beta$accept
    }
    if('gamma' %in% steps)
    {
      samp.gamma = sample.gamma(y, x, theta, psi, tuning.gamma, burn)
      theta[2] = samp.gamma$gamma
      tuning.gamma <- tune.all[i,4] <- samp.gamma$tuning
      accept.gamma[i] = samp.gamma$accept
    }
    if('nu' %in% steps)
    {
      samp.nu = sample.nu(y, x, theta, psi, tuning.nu, burn)
      theta[3] = samp.nu$nu
      tuning.nu <- tune.all[i,5] <- samp.nu$tuning
      accept.nu[i] = samp.nu$accept
    }

    # Only save every n.thin iteration
    if (ii <- save.iteration(i,n.burn,n.thin)) {
      keep.x[ii,,] = x
      keep.theta[ii,] = theta
    }
  }

  return(list(x=keep.x, theta=keep.theta, tuning=tune.all, accept.x=accept.x, accept.theta=cbind(accept.beta,accept.gamma,accept.nu)))
}

rmove <- function(y, x, theta, psi, tuning)
{
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)

  # Move states
  x = sample.x(y, x, theta, psi, tuning[1:2])
  
  # Move parameters
  theta[1] = sample.beta(y, x, theta, psi, tuning[3])
  theta[2] = sample.gamma(y, x, theta, psi, tuning[4])
  theta[3] = sample.nu(y, x, theta, psi, tuning[5])
  
  return(list(state=x, theta=theta))
}

sample.x <- function(y, x, theta, psi, tuning, burn = FALSE)
{
  # Check dimensions
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  accept = rep(FALSE, nt + 1)
  
  # resample prior state
  x.prop = rx(x[,1], tuning)
  logMH = dlevo(x[,2], x.prop, theta, psi$P) + dlx0(x.prop) - dlevo(x[,2], x[,1], theta, psi$P) - dlx0(x[,1])
  if(log(runif(1)) < logMH)
  {
    x[,1] = x.prop
    accept[1] = TRUE
    if(burn) tuning = tuning*1.1
  } else {
    if(burn) tuning = tuning / 1.1 
  }
  
  # resample middle states
  for(i in 2:nt)
  {
    x.prop = rx(x[,i], tuning)
    logMH = dllik(y[,i-1], x.prop, psi$b, psi$varsigma, psi$sigma, psi$eta) + dlevo(x[,i+1], x.prop, theta, psi$P) + dlevo(x.prop, x[,i-1], theta, psi$P) - dllik(y[,i-1], x[,i], psi$b, psi$varsigma, psi$sigma, psi$eta) - dlevo(x[,i+1], x[,i], theta, psi$P) - dlevo(x[,i], x[,i-1], theta, psi$P)
    if(log(runif(1)) < logMH)
    {
      x[,i] = x.prop
      accept[i] = TRUE
      if(burn) tuning = tuning*1.1
    } else {
      if(burn) tuning = tuning / 1.1 
    }
  }
  
  # resample last state
  x.prop = rx(x[,nt+1], tuning)
  logMH = dllik(y[,nt], x.prop, psi$b, psi$varsigma, psi$sigma, psi$eta) + dlevo(x.prop, x[,nt], theta, psi$P) - dllik(y[,nt], x[,nt+1], psi$b, psi$varsigma, psi$sigma, psi$eta) - dlevo(x[,nt+1], x[,nt], theta, psi$P)
  if(log(runif(1)) < logMH)
  {
    x[,nt+1] = x.prop
    accept[nt + 1] = TRUE
    if(burn) tuning = tuning*1.1
  } else {
    if(burn) tuning = tuning / 1.1
  }
  return(list(x=x,accept=accept,tuning=tuning))
}

sample.beta <- function(y, x, theta, psi, tuning, burn = FALSE)
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
  if(log(runif(1)) < logMH)
  { 
    beta.curr = beta.prop
    accept = TRUE
    if(burn) tuning = tuning*1.1
  } else {
    if(burn) tuning = tuning / 1.1
  }
  
  return(list(beta=beta.curr,accept=accept,tuning=tuning))
}

sample.gamma <- function(y, x, theta, psi, tuning, burn = FALSE)
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
  if(log(runif(1)) < logMH)
  {
    gamma.curr = gamma.prop
    accept = TRUE
    if(burn) tuning = tuning*1.1
  } else {
    if(burn) tuning = tuning / 1.1
  }

  return(list(gamma=gamma.curr,accept=accept,tuning=tuning))
}

sample.nu <- function(y, x, theta, psi, tuning, burn = FALSE)
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
  if(log(runif(1)) < logMH)
  {
    nu.curr = nu.prop
    accept = TRUE
    if(burn) tuning = tuning*1.1
  } else {
    if(burn) tuning = tuning / 1.1
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
