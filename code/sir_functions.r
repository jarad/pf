# sir_functions.R - functions needed for simulating data and running particle filters for SIR model of a disease epidemic
#
# robs - function that produces a simulated observation y given current state x
# Arguments:
# x - 2-element vector, current state (s,i)
# b - vector of length L, the number of syndromes, with elements values of parameters b_l
# varsigma - vector of length L, the number of syndromes, with elements values of parameters varsigma_l
# sigma - vector of length L, number of syndromes, with elements values of parameters sigma_l
# eta - vector of length L, number of syndromes, with elements values of parameters eta_l
robs = function(x,b,varsigma,sigma,eta)
{
  L = length(b)
  if(!(L == length(varsigma) & L == length(sigma) & L == length(eta))) stop("b, varsigma, sigma, and eta must all have same length")
  l = sample(1:L,sample(0:L,1))
  y = rep(NA,L)
  if(length(l) > 0) y[l] = rlnorm(length(l),b[l]*x[2]^varsigma[l]+eta[l],sigma[l])
  return(y)
}

# rinit - function to initialize values of initial state
# Arguments:
# i0 - scalar between 0 and 1, initial proportion of infected individuals in population
rinit = function(i0=.002)
{
  if(i0 < 0 | i0 > 1) stop("i0 must be between 0 and 1")
  return(c(1-i0,i0))
}

# revo - function that propagates state forward given previous state x and parameters theta
# Arguments:
# x - 2-element vector, current state (s,i)
# theta - 3-element vector, current parameter (beta,gamma,nu)
# P - scalar, population size
# random - boolean, if TRUE function returns a sample from predictive distribution of the next state and if FALSE function returns the mean of next state
revo = function(x,theta,P,random=TRUE,sdfac=1)
{
  tmp = c(NA, NA)
  tmp[1] = x[1] - theta[1]*x[2]*x[1]^theta[3] + random*rnorm(1, 0, sdfac*sqrt(theta[1] / P^2))
  tmp[2] = x[2] + theta[1]*x[2]*x[1]^theta[3] - x[2]*theta[2] + random*rnorm(1, -1*(tmp[1] - (x[1] - theta[1]*x[2]*x[1]^theta[3])), sdfac*sqrt(theta[2] / P^2))
  while(!(all(tmp >= 0) & sum(tmp) <= 1))
  {
    tmp[1] = x[1] - theta[1]*x[2]*x[1]^theta[3] + random*rnorm(1, 0, sdfac*sqrt(theta[1] / P^2))
    tmp[2] = x[2] + theta[1]*x[2]*x[1]^theta[3] - x[2]*theta[2] + random*rnorm(1, -1*(tmp[1] - (x[1] - theta[1]*x[2]*x[1]^theta[3])), sdfac*sqrt(theta[2] / P^2))
  }
  x[1:2] = tmp
  return(x)
}

# dlevo - function to evaluate the logarithm of the state transition density
# Arguments:
# Arguments:
# x - 2-element vector of state (s,i) at which to evaluate the density
# x.curr - 2-element vector, current state (s,i)
# theta - 3-element vector, current parameter (beta,gamma,nu)
# P - scalar, population size
dlevo <- function(x,x.curr,theta,P,sdfac=1)
{
  if(all(x >= 0) & sum(x) <= 1)
  {
    dnorm(x[1], x.curr[1] - theta[1]*x.curr[2]*x.curr[1]^theta[3], sdfac*sqrt(theta[1]/P^2), log=TRUE) + dnorm(x[2], x.curr[2]*(1 - theta[2]) + x.curr[1] - x[1], sdfac*sqrt(theta[2] / P^2), log=TRUE)
  } else { return(-Inf) }
}

# dllik - function to return the log of the likelihood function given current observation y, current state x, and parameter theta
# Arguments:
# y - vector, current observation of length L, the number of syndromes; may have empty elements
# x - vector, current state (s,i)
# b - vector of length L, values of b_l's
# varsigma - vector of length L, values of varsigma_l's
# sigma - vector of length L, values of sigma_l's
# eta - vector of length L, values of eta_l's
dllik = function(y,x,b,varsigma,sigma,eta)
{
  L = length(y)
  if(!(L == length(b) & L == length(varsigma) & L == length(sigma) & L == length(eta))) stop("b, varsigma, sigma, and eta must all have same length")
  h = b*x[2]^varsigma+eta
  return(sum(dlnorm(y,h,sigma,log=T),na.rm=TRUE))
}

# rprior - function that samples from the prior distribution of the initial states and unknown parameters; returns a list with vector elements initial state vector x and parameter values theta
# Arguments:
# rtheta - a function that takes integer argument particle number and returns a vector of sampled values from the prior distribution of the unknown parameters; sampled values are assumed to be already mapped to the real line
rprior = function(rtheta)
{
  theta0 = rtheta()
  i0 = rnorm(1,0.002,0.0005)
  while(i0 < 0 | i0 > 1) i0 = rnorm(1,0.002,0.0005)
  s0 = 1 - i0
  return(list(x=c(s0,i0),theta=theta0))
}

# dlx0 - function that evaluates the log density of the prior state (s,i)
dlx0 <- function(x)
{
  if(all(x >= 0) & sum(x) <= 1)
  {
    dnorm(x[2], .002, .0005, log=TRUE)
  } else { return(-Inf) }
}

# dlbeta0 - function that evaluates the log density of the prior on beta
dlbeta0 <- function(beta)
{
  log.params = find.mu.sigma(.14, .50)
  return(dlnorm(beta, log.params$mu, log.params$sigma, log=TRUE))
}

# dlgamma0 - function that evaluates the log density of the prior on gamma
dlgamma0 <- function(gamma)
{
  log.params = find.mu.sigma(.09, .143)
  return(dlnorm(gamma, log.params$mu, log.params$sigma, log=TRUE))
}

# dlnu0 - function that evaluates the log density of the prior on nu
dlnu0 <- function(nu)
{
  log.params = find.mu.sigma(0.95, 1.3)
  dlnorm(nu, log.params$mu, log.params$sigma, log=TRUE)
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

# Functions to find mean and sd on log scale so that (1-alpha)*100% of random draws fall b/w a and b
find.mu.sigma <- function(a, b, alpha = 0.05)
{
  z = qnorm(1 - alpha/2)
  sigma = log(b/a) / (2*z)
  mu = log(b) - z*sigma
  return(list(mu=mu,sigma=sigma))
}
