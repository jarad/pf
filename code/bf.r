#' A function to perform the bootstrap filter
#'
#' @param y a matrix with no, dimension of each observation, rows and  nt, number of observation time points, columns
#' @param dllik a function to evaluate the logarithm of the likelihood with arguments y and x for data and state respectively
#' @param revo a function to propagate the state with argument the state x
#' @param rprior a function to sample for the prior for the state
#' @param n the number of particles
#' @param progress a boolean to display progress bar if TRUE
#' @param seed a boolean to set the seed before sampling prior state
#' @param ... arguments passed on to resample
#' @return a list containing an n x (nt+1) matrix of normalized particles weights, a ns x n x (nt+1) array of state draws, and an n x nt parent matrix
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @export bf
#' @references Gordon, N. J.; Salmond, D. J. and Smith, A. F. M. (1993). 
#'    "Novel approach to nonlinear/non-Gaussian Bayesian state estimation". 
#'    IEEE Proceedings F on Radar and Signal Processing 140 (2): 107–113. 
#' @seealso \code{\link{resample}}
#'
bf = function(y, dllik, revo, rprior, n, progress = TRUE, seed = TRUE, ...)
{
  require(smcUtils)

  if (!is.matrix(y)) y = matrix(y, 1)
  no = nrow(y) # not currently used
  nt = ncol(y)

  # Find dimension of state
  current.seed = .Random.seed
  ns = length(rprior())
  .Random.seed = current.seed

  # Set up initial state
  state = array(NA, dim=c(ns,n,nt+1))
  for (j in 1:n)
  {
    if(seed) set.seed(j)
    state[,j,1] = rprior()
  }

  # Initialize weights
  weight = matrix(NA, n, nt+1)
  weight[,1] = 1/n
  
  # Initialize parent
  parent = matrix(NA, n, nt+1)

  # Run bootstrap filter
  if(progress) pb = txtProgressBar(0,nt,style=3)
  for (i in 1:nt) 
  {
    if(progress) setTxtProgressBar(pb,i)
    if (i==1) # No resampling
    {
      kk = 1:n
      tmp = list(weights=rep(1/n,n))
    } else 
    {
      tmp = resample(weight[,i],...)
      kk = tmp$indices
      tmp$weights
    }

    for (j in 1:n) 
    {
      state[,j,i+1] = revo(state[,kk[j],i])
      weight[j,i+1] = log(tmp$weights[kk[j]]) + 
                      dllik(y[,i], state[,j,i+1])
    }

    weight[,i+1] = renormalize(weight[,i+1], log=T)

    parent[,i+1] = kk
  }

  return(list(state=state, weight=weight, parent=parent))
}