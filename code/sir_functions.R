# Functions for simulating data; revo also needed for particle filters
# revo - function that propagates state forward given previous state x and, optionally, parameters theta; intended to be revo argument in ss.sim, bf, apf, or kd_pf functions
# Arguments:
# x - vector, current state
# P - scalar, population size
# d - scalar, positive integer indicating number of times to run state equation in 1 day before predicting the next state; larger d gives better approximation of differential equations specifying disease dynamics
# theta - optional vector, current values of parameters
# thetal - optional vector, prior lower bounds on parameters; must be given if theta given
# thetau - optional vector, prior upper bounds on parameters; must be given if theta given
# stateonly - boolean, if TRUE theta is ignored; if FALSE, theta must be given
# random - boolean, if TRUE function returns a sample from predictive distribution; if FALSE function returns a point estimate of next state
revo = function(x,P,d=1,theta=NULL,thetal=NULL,thetau=NULL,stateonly=TRUE,random=TRUE)
{
  tau = 1/d
  if(stateonly)
  {
    x1var = (x[3] + x[4])*tau^2 / P^2
    x2var = x[3]*tau^2 / P^2
    for(i in 1:d)
    { 
      is = -1
      ss = -1
      while(!(is >= 0 & is <= 1 & ss >= 0 & ss <= 1))
      {
        is = x[1] + tau*x[1]*(x[3]*x[2]^x[5] - x[4]) + random*rnorm(1,0,sqrt(x1var))
        ss = x[2] - tau*x[3]*x[1]*x[2]^x[5] + random*rnorm(1,0,sqrt(x2var))
      }
      x[1] = is
      x[2] = ss
    }
  } else {
    theta[1] = theta2u(theta[1],thetal[1],thetau[1])
    theta[2] = theta2u(theta[2],thetal[2],thetau[2])
    theta[3] = theta2u(theta[3],thetal[3],thetau[3])

    x1var = (theta[1] + theta[2])*tau^2 / P^2
    x2var = theta[1]*tau^2 / P^2
    for(i in 1:d)
    { 
      is = -1
      ss = -1
      while(!(is >= 0 & is <= 1 & ss >= 0 & ss <= 1))
      {
        is = x[1] + tau*x[1]*(theta[1]*x[2]^theta[3] - theta[2]) + random*rnorm(1,0,sqrt(x1var))
        ss = x[2] - tau*theta[1]*x[1]*x[2]^theta[3] + random*rnorm(1,0,sqrt(x2var))
      }
      x[1] = is
      x[2] = ss
    }
  }
  return(x)
}

# robs - function that produces a simulated observation y given current state x; intended to be robs argument in ss.sim function
# Arguments:
# x - vector, current state
# b - vector of length J = number of syndromes with elements values of parameters b_j
# varsigma - vector of length J = number of syndromes with elements values of parameters varsigma_j
# sigma - vector of length J = number of syndromes with elements values of parameters sigma_j
robs = function(x,b,varsigma,sigma)
{
  j = length(b)
  if(!(j == length(varsigma) & j == length(sigma))) stop("b, varsigma, and sigma must all have same length")
  J = sample(1:j,1)
  y = rep(NA,j)
  y[J] = rnorm(1,varsigma[J]*log(b[J]*x[1]),sigma[J])
  return(y)
}

# rinit - function to initialize values of initial states and unknown parameters; intended to be rinit function in ss.sim
# Arguments:
# i0 - scalar between 0 and 1, initial proportion of infected individuals in population
# theta - vector of true values of unknown parameters
rinit = function(i0,theta)
{
  if(i0 < 0 | i0 > 1) stop("i0 must be between 0 and 1")
  return(c(i0,1-i0,theta))
}

# Functions to reparameterize theta to and from the real line
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

# dllik - function to return the log of the likelihood function given current observation y, current state x, and optionally parameters theta; intended to be dllik argument in bf, apf, or kd_pf functions
# Arguments:
# y - vector, current observation of length J = number of syndromes with all elements NA except for element j = syndrome from which current observation arrives
# x - vector, current state
# b - vector of length J, values of b_j's if known; ignored if addparam or stateonly FALSE
# varsigma - vector of length J, values of varsigma_j's if known; ignored if addparam or stateonly FALSE
# sigma - vector of length J, values of sigma_j's if known; ignored if addparam or stateonly FALSE
# theta - optional vector, current values of parameters
# thetal - optional vector, prior lower bounds on parameters; must be given if theta given
# thetau - optional vector, prior upper bounds on parameters; must be given if theta given
# stateonly - boolean, if TRUE theta is ignored; if FALSE theta must be given
# addparam - boolean, if TRUE b, varsigma, and sigma assumed unknown and included in either x or theta; if FALSE b, varsigma, and sigma must be given
dllik = function(y,x,b,varsigma,sigma,theta=NULL,thetal=NULL,thetau=NULL,stateonly=TRUE,addparam=FALSE)
{
  J = which(!is.na(y))
  if (!(length(J) == 1)) stop("y must be a vector with all but 1 element NA")
  if(stateonly)
  {
    if(addparam)
    {
      h = x[7]*log(x[6]*x[1])
      s = x[8]
    } else {
      h = varsigma[J]*log(b[J]*x[1])
      s = sigma[J]
    }
  } else {
    theta[4] = theta2u(theta[4],thetal[4],thetau[4])
    theta[5] = theta2u(theta[5],thetal[5],thetau[5])
    theta[6] = theta2u(theta[6],thetal[6],thetau[6])
    h = theta[5]*log(theta[4]*x[1])
    s = theta[6]
  }
  return(dnorm(y[J],h,s,log=T))
}

# rprior - function that samples from prior distributions of initial states and unknown parameters and returns a state vector or a list with initial state vector x and parameter values theta; intended to be rprior argument in functions bf, apf, or kd_pf
# Arguments:
# sim - list with an element called y that is a J x nt matrix of observations, where J is the number of syndromes, nt is the number of total time points at which observations were given, and each column y is defined as the argument y in function dllik (above)
# thetal - vector, prior lower bounds on unknown parameters
# thetau - vector, prior upper bounds on unknown parameters
# b - vector of length J, values of b_j's if known; ignored if addparam FALSE
# varsigma - vector of length J, values of varsigma_j's if known; ignored if addparam FALSE
# sigma - vector of length J, values of sigma_j's if known; ignored if addparam FALSE
# stateonly - boolean, if TRUE only a state vector is returned; if FALSE, a list with state vector x and parameters theta is returned
# addparam - boolean, if TRUE initial values of unknown b_j's, varsigma_j's and sigma_j's are included in returned x or theta (depending on value of stateonly)
rprior = function(sim,thetal,thetau,b=NULL,varsigma=NULL,sigma=NULL,stateonly=TRUE,addparam=FALSE)
{
  J = which(!is.na(sim$y[,1]))
  if (!(length(J) == 1)) stop("y must be a vector with all but 1 element NA")
  
  beta0 = runif(1,thetal[1],thetau[1])
  gamma0 = runif(1,thetal[2],thetau[2])
  nu0 = runif(1,thetal[3],thetau[3])
  if(addparam)
  {
    b0 = runif(1,thetal[4],thetau[4])
    varsigma0 = runif(1,thetal[5],thetau[5])
    sigma0 = runif(1,thetal[6],thetau[6])
    tauj = sim$y[J,1] - log(b0) - 1
    while(tauj < sim$y[J,1] - log(b0)) tauj = rnorm(1,0,sigma0)
    i0 = exp(sim$y[J,1] - tauj)/b0
    s0 = 1 - i0
  } else {
    tauj = sim$y[J,1] - log(b[J]) - 1
    while(tauj < sim$y[J,1] - log(b[J])) tauj = rnorm(1,0,sigma[J])
    i0 = exp(sim$y[J,1] - tauj)/b[J]
    s0 = 1 - i0
  }
  if(stateonly & addparam){
    return(c(i0,s0,beta0,gamma0,nu0,b0,varsigma0,sigma0))
  } else if(stateonly) {
    return(c(i0,s0,beta0,gamma0,nu0))
  } else if(addparam) {
    return(list(x=c(i0,s0),theta=u2theta(c(beta0,gamma0,nu0,b0,varsigma0,sigma0),thetal,thetau)))
  } else {
    return(list(x=c(i0,s0),theta=u2theta(c(beta0,gamma0,nu0),thetal[1:3],thetau[1:3])))
  }
}
