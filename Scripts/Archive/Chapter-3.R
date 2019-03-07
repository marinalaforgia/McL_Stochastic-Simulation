#### Chapter 3: Stochastic Population Models ####
rm(list = ls())

#### Read in Data ####


#### Define Terms ####


#### Write Function to evaluate population numbers across years ####

Annual <- function(A, tf, N0, p, q, df, dg, afg, agf, c) {
  n=dim(A)[3]
  P = cbind(c(1-p, q), c(p, 1-q)) # Create matrix of probabilities of switching between good and bad years
  env = 1 # set initial environmental state
  N = matrix(NA,tf,2) # create matrix of zeros for storing population size for forbs and grasses
  N[1,] = N0 # set initial value
  vec <- rep(0,tf) # empty vector to keep track of environmental state
  vec[1] <- env
  for(t in 1:(tf-1)) { # loop through each time step up to nsteps
    N[t+1,1] = (1-A[2,2,env])*(1-df)*N[t,1] + A[1,2,env]*A[2,2,env]*N[t,1]/(c+A[2,2,env]*N[t,1]+afg*A[2,1,env]*N[t,2]) # forb 
    N[t+1,2] = (1-A[2,1,env])*(1-dg)*N[t,2] + A[1,1,env]*A[2,1,env]*N[t,2]/(c+A[2,1,env]*N[t,2]+agf*A[2,2,env]*N[t,1]) # grass
    env=sample(1:n,size=1,prob=P[env,])
    vec[t+1] <- env #env 1 is bad, env 2 is good
  }
  return(list(vec = vec,N = N)) # return population size and env state vectors
}