# Recreate figure from Levine & Rees 2004

# two species competition model with environmental stochasticity; all models start with two individuals for each species
rm(list=ls())

# 1. Write Function to evaluate population numbers across years
# 2. Define matrices of lambda/g in good/bad/intermediate year
# 3. Put matrices together in an array
# 4. Define Parameter Values
# 5. Plot density across years with
# a. equal probability of switching between years (graph for germinated and ungerminated fraction, and graph for just germinated fraction)
# b. constantly favorable environment 
# c. constantly unfavorable environment
# d. intermediate environment 

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

# Define matrix of lambdas and germination rates together
A.bad.hd <- matrix(c(5,5,0.9,0.1), 2, byrow = T) # row 1 - lambda, row 2, germination, col 1 grass, col 2 - high dormant forb
A.bad.ld <- matrix(c(5,5,0.9,0.5), 2, byrow = T) # row 1 - lambda, row 2, germination, col 1 grass, col 2 - low dormant forb
A.good.hd <- matrix(c(30, 30,0.9,0.6), 2, byrow = T)
A.good.ld <- matrix(c(30, 30,0.9,0.8), 2, byrow = T)
A.inter <- matrix(c(17.5, 17.5, 0.9, 0.7), 2, byrow = T)

# Put matrices together in array
A.ld = array(0,dim = c(2,2,2)) #2 2x2 matrices
A.ld[,,1] = A.bad.ld
A.ld[,,2] = A.good.ld

A.hd = array(0,dim = c(2,2,2)) #2 2x2 matrices
A.hd[,,1] = A.bad.hd
A.hd[,,2] = A.good.hd

A2 = array(0,dim = c(2,2,2)) #Probably a more efficient way to do this to look at a non-varying environment but I couldnt figure it out
A2[,,1] = A.inter
A2[,,2] = A.inter

# Define Parameter Values
df = 0.1
dg = 0.7
afg = 2
agf = 0.5
c = 1
N0 = c(2,2)
tf = 10

par(mfrow=c(2,1))
# equal probability of good and bad years (germ and ungerm fraction); low dormancy forb
eq <- Annual(tf = 100, N0 = N0, p = 0.5, q = 0.5, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A.ld)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Germinated and Ungerminated Fraction")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

# same graph only germinated fraction 
N.germ <- eq[[2]]
N.germ[,2] <- eq[[2]][,2] - eq[[2]][,2]*.1 # grass germinated fraction
N.germ[,1] <- ifelse(eq[[1]]>1, eq[[2]][,1] - eq[[2]][,1]*.3, eq[[2]][,1] - eq[[2]][,1]*.9) # forb germinated fraction in good year and bad year
matplot(N.germ,type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Germinated Fraction Only")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

par(mfrow=c(1,3))
# constantly favorable year
eq <- Annual(tf = 100, N0 = N0,  p = 1, q = 0, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A.ld) 
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Constantly Favorable Environment")
legend("topleft", legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

# constantly unfavorable year
eq <- Annual(tf = 100, N0 = N0,  p = 0, q = 1, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A.ld) 
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Constantly Unfavorable Environment")
legend("topleft", legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

# constantly intermediate year
eq <- Annual(tf = 100, N0 = N0,  p = 0, q = 1, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A2)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Constantly Intermediate Environment")
legend("topleft", legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

par(mfrow=c(2,1))
# equal probability of good and bad years (germ and ungerm fraction); high dormancy forb
eq <- Annual(tf = 100, N0 = N0, p = 0.5, q = 0.5, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A.hd)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Germinated and Ungerminated Fraction")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

# same graph only germinated fraction 
N.germ <- eq[[2]]
N.germ[,2] <- eq[[2]][,2] - eq[[2]][,2]*.1 # grass germinated fraction
N.germ[,1] <- ifelse(eq[[1]]>1, eq[[2]][,1] - eq[[2]][,1]*.3, eq[[2]][,1] - eq[[2]][,1]*.9) # forb germinated fraction in good year and bad year
matplot(N.germ,type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Germinated Fraction Only")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

par(mfrow=c(1,3))
# constantly favorable year
eq <- Annual(tf = 100, N0 = N0,  p = 1, q = 0, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A.hd) 
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Constantly Favorable Environment")
legend("topleft", legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

# constantly unfavorable year
eq <- Annual(tf = 100, N0 = N0,  p = 0, q = 1, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A.hd) 
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Constantly Unfavorable Environment")
legend("topleft", legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

# constantly intermediate year
eq <- Annual(tf = 100, N0 = N0,  p = 0, q = 1, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A2)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Constantly Intermediate Environment")
legend("topleft", legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n")

