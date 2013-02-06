# 'ss.sim' - a function to simulate data from a state space model of an epidemic outbreak
#
# param nt a scalar, number of days for which to generate observations
# param revo a function to generate next state given previous state x
# param robs a function to generate next observation given current state x
# param rinit a function to generate initial state vector
# returns a list containing an ns x nt + 1 matrix of states (x) and an no x nt matrix of generated observations (y)
#
ss.sim = function(nt,revo,robs,rinit)
{
  # Find dimension of state
  current.seed = .Random.seed
  ns = length(rinit())
  .Random.seed = current.seed

  # Find dimension of observations
  current.seed = .Random.seed
  no = length(robs(rinit()))
  .Random.seed = current.seed

  # Initialize state matrix, observation vector
  x = matrix(NA,ns,nt+1)
  x[,1] = rinit()
  y = matrix(NA,no,nt)

  # Generate states and observations
  pb = txtProgressBar(0,nt,style=3)
  for(i in 1:nt)
  {
    setTxtProgressBar(pb,i)
    x[,i+1] = revo(x[,i])
    y[,i] = robs(x[,i+1])
  }

  return(list(x=x,y=y))
}