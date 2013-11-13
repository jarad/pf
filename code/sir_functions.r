# sir_functions.R - functions needed for simulating data and running particle filters for SIR model of a disease epidemic
#
# revo - function that propagates state forward given previous state x and parameters theta
# Arguments:
# x - 2-element vector, current state (i,s)
# theta - 3-element vector, current parameter (beta,gamma,nu)
# P - scalar, population size
# d - scalar, positive integer indicating number of times to propagate state within 1 time unit before returning
# random - boolean, if TRUE function returns a sample from predictive distribution of the next state and if FALSE function returns the mean of next state
revo = function(x,theta,P,random=TRUE,d=1)
{
  d = floor(d)
  if(d < 1) stop("d must be a positive integer")
  tau = 1/d
  x1var = (theta[1] + theta[2])*tau^2 / P^2
  x2var = theta[1]*tau^2 / P^2
  x1x2cov = -theta[1]*tau^2 / P^2
  vr = cbind(c(x1var,x1x2cov),c(x1x2cov,x2var))
  cl = t(chol(vr))
  for(i in 1:d)
  { 
    tmp = c(-1,-1)
    while(!(all(tmp >= 0) & sum(tmp) <= 1))
    {
      tmp[1] = x[1] + tau*(x[1]*theta[1]*x[2]^theta[3] - x[1]*theta[2])
      tmp[2] = x[2] - tau*theta[1]*x[1]*x[2]^theta[3]
      if(random) tmp = tmp + cl%*%rnorm(2)
    }
    x[1:2] = tmp
  }
  return(x)
}

rmove <- function(y, x, theta, b, varsigma, sigma, eta, P)
{
  require(rjags)

  # Set data values
  N = dim(y)[2]
  stopifnot(N == dim(x)[2] - 1)
  L = length(b)
  stopifnot(L == length(varsigma) & L == length(sigma) & L == length(eta) & L == dim(y)[1])

  # Create JAGS model
  d = list(y=y,N=N,L=L,b=b,varsigma=varsigma,sigma=sigma,eta=eta,P=P)
  x[2,1] = NA
  inits = list(list(beta=theta[1],gamma=theta[2],nu=theta[3],x=x))
  mod = jags.model("jags-model.txt", data=d, n.chains=1, n.adapt=0, inits=inits, quiet=TRUE)
  
  # Generate new samples
  samps = coda.samples(mod, c("beta","gamma","nu","x"), n.iter=1)
  
  # Return moved particle
  return(list(state=matrix(samps[[1]][4:(2*(N+1)+3)],nr=2,byrow=FALSE),theta=samps[[1]][1:3]))
}

# robs - function that produces a simulated observation y given current state x
# Arguments:
# x - 2-element vector, current state (i,s)
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
  if(length(l) > 0) y[l] = rlnorm(length(l),b[l]*x[1]^varsigma[l]+eta[l],sigma[l])
  return(y)
}

# rinit - function to initialize values of initial state
# Arguments:
# i0 - scalar between 0 and 1, initial proportion of infected individuals in population
rinit = function(i0=.002)
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
# eta - vector of length L, values of eta_l's
dllik = function(y,x,b,varsigma,sigma,eta)
{
  L = length(y)
  if(!(L == length(b) & L == length(varsigma) & L == length(sigma) & L == length(eta))) stop("b, varsigma, sigma, and eta must all have same length")
  h = b*x[1]^varsigma+eta
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






