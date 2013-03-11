# sir_functions.R - functions needed for simulating data, running particle filters, calculating quantiles of filtered distributions, and constructing plots for SIR model of a disease epidemic
#
# revo - function that propagates state forward given previous state x and parameters theta
# Arguments:
# x - 2-element vector, current state (i,s)
# P - scalar, population size
# d - scalar, positive integer indicating number of times to propagate state within 1 time unit before returning
# theta - 3-element vector, current parameter (beta,gamma,nu)
# random - boolean, if TRUE function returns a sample from predictive distribution of the next state and if FALSE function returns the mean of next state
revo = function(x,P,d=1,theta,random=TRUE)
{
  d = floor(d)
  if(d < 1) stop("d must be a positive integer")
  tau = 1/d
  x1var = (theta[1] + theta[2])*tau^2 / P^2
  x2var = theta[1]*tau^2 / P^2
  for(i in 1:d)
  { 
    is = -1
    ss = -1
    entered = FALSE # in case of numerical instability
    while(!(is >= 0 & is <= 1 & ss >= 0 & ss <= 1))
    {
      is = x[1] + tau*x[1]*(theta[1]*x[2]^theta[3] - theta[2]) + random*rnorm(1,0,sqrt(x1var))
      ss = x[2] - tau*theta[1]*x[1]*x[2]^theta[3] + random*rnorm(1,0,sqrt(x2var))
      if(entered == TRUE & !random)
      {
        if(is < 0) is = 0
        if(ss < 0) ss = 0
        if(is > 1) is = 1
        if(ss > 1) ss = 1
      }
      entered = TRUE
    }
    x[1] = is
    x[2] = ss
  }
  return(x)
}

# robs - function that produces a simulated observation y given current state x
# Arguments:
# x - 2-element vector, current state (i,s)
# b - vector of length L, the number of syndromes, with elements values of parameters b_l
# varsigma - vector of length L, the number of syndromes, with elements values of parameters varsigma_l
# sigma - vector of length L, number of syndromes, with elements values of parameters sigma_l
# dpower - scalar, power of b_l*i^varsigma_l in the denominator of the observation variance
robs = function(x,b,varsigma,sigma,dpower=2)
{
  L = length(b)
  if(!(L == length(varsigma) & L == length(sigma))) stop("b, varsigma, and sigma must all have same length")
  l = sample(1:L,sample(0:L,1))
  y = rep(NA,L)
  if(length(l) > 0) y[l] = rnorm(length(l),varsigma[l]*log(b[l]*x[1]),sigma[l]/sqrt(b[l]*x[1]^varsigma[l])^dpower)
  return(y)
}

# rinit - function to initialize values of initial state
# Arguments:
# i0 - scalar between 0 and 1, initial proportion of infected individuals in population
rinit = function(i0=.001)
{
  if(i0 < 0 | i0 > 1) stop("i0 must be between 0 and 1")
  return(c(i0,1-i0))
}

# dllik - function to return the log of the likelihood function given current observation y, current state x, and parameter theta
# Arguments:
# y - vector, current observation of length L, the number of syndromes; may have empty elements
# x - vector, current state (i,s)
# b - vector of length L, values of b_l's
# varsigma - vector of length L, values of varsigma_l's
# sigma - vector of length L, values of sigma_l's
# dpower - scalar, power of b_l*i^varsigma_l in the denominator of the observation variance
dllik = function(y,x,b,varsigma,sigma,dpower=2)
{
  L = length(y)
  if(!(L == length(b) & L == length(varsigma) & L == length(sigma))) stop("b, varsigma, and sigma must all have same length")
  h = log(b*x[1]^varsigma)
  s = sigma / sqrt(b*x[1]^varsigma)^dpower
  return(sum(dnorm(y,h,s,log=T),na.rm=TRUE))
}

# rprior - function that samples from the prior distribution of the initial states and unknown parameters; returns a list with vector elements initial state vector x and parameter values theta
# Arguments:
# rtheta - a function that takes no arguments and returns sampled values from the prior distribution of beta, gamma, and nu; sampled values are assumed to be already mapped to the real line
# obsparam - boolean, if TRUE b_j's, varsigma_j's and sigma_j's assumed unknown and prior samples are included in returned theta
# rthetaPlus - a function that takes no arguments and returns a list with elements b, varsigma, and sigma, each vectors containing sampled values from the prior distribution of b, varsigma and sigma; sampled values are assumed to be already mapped to the real line; ignored if obsparam = FALSE
rprior = function(rtheta,obsparam=FALSE,rthetaPlus=NULL)
{
  .current.seed = .Random.seed
  theta0 = rtheta()
  i0 = -1
  while(i0 < 0 | i0 > 1)
  {
    .Random.seed = .current.seed
    i0 = rnorm(1,0.002,0.00255)
    .Random.seed = .current.seed
    set.seed(sample(1:1000,1))
  }
  s0 = 1 - i0
  if(obsparam)
  {
    addparam = rthetaPlus()
    theta0 = c(theta0,addparam$b,addparam$varsigma,addparam$sigma)
  }
  return(list(x=c(i0,s0),theta=theta0))
}

# Functions to reparameterize theta to [a,b] from the real line and vice versa
theta2u = function(theta,a,b)
{
  etheta = exp(theta)
  return((b*etheta + a) / (1 + etheta))
}
u2theta = function(u,a,b)
{
  U = (u - a) / (b - a)
  return(log(U / (1 - U)))
}






