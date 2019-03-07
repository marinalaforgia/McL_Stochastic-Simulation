# Changing germination rates based on previous year
# Recreate figure from Levine & Rees 2004
#Germination rates are a function of yield the PREVIOUS year
# changing competition coefficients

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

Annual <- function(A, tf, N0, p, q, df, dg, c) {
  n=dim(A)[3]
  P = cbind(c(1-p, q), c(p, 1-q)) # Create matrix of probabilities of switching between good and bad years
  env = 1 # set initial environmental state
  N = matrix(NA,tf,2) # create matrix of zeros for storing population size for forbs and grasses
  N[1,] = N0 # set initial value
  vec <- rep(0,tf) # empty vector to keep track of environmental state
  vec[1] <- env
  for(t in 1:(tf-1)) { # loop through each time step up to nsteps
    N[t+1,1] = (1-A[2,2,env])*(1-df)*N[t,1] + A[1,2,env]*A[2,2,env]*N[t,1]/(c+A[2,2,env]*N[t,1]+A[3,1,env]*A[2,1,env]*N[t,2]) # forb 
    N[t+1,2] = (1-A[2,1,env])*(1-dg)*N[t,2] + A[1,1,env]*A[2,1,env]*N[t,2]/(c+A[2,1,env]*N[t,2]+A[3,2,env]*A[2,1,env]*N[t,2]) # grass
    env=sample(1:n,size=1,prob=P[env,])
    vec[t+1] <- env #env 1 is bad, env 2 is good
  }
  return(list(vec = vec,N = N)) # return population size and env state vectors
}

# Define matrix of lambdas and germination rates together
A.bad.SA <- matrix(c(5,5,1,0.1,1,0), 3, byrow = T) # row 1 - lambda, row 2, germination, row 3, competition, col 1 grass, col 2 - high dormant forb
A.bad.ST <- matrix(c(5,10,1,0.5,0.5,0), 3, byrow = T) # row 1 - lambda, row 2, germination, col 1 grass, col 2 - low dormant forb
A.good.SA <- matrix(c(30, 30,1,0.4,1.5,0), 3, byrow = T) # think LACA
A.good.ST <- matrix(c(30, 20,1,0.5,3,0), 3, byrow = T) # think AGHE
#A.inter <- matrix(c(17.5, 17.5, 0.9, 0.7), 3, byrow = T)

# Put matrices together in array
A.ST = array(0,dim = c(3,2,2)) #2 2x2 matrices
A.ST[,,1] = A.bad.ST
A.ST[,,2] = A.good.ST

A.SA = array(0,dim = c(3,2,2)) #2 2x2 matrices
A.SA[,,1] = A.bad.SA
A.SA[,,2] = A.good.SA

#A2 = array(0,dim = c(2,2,2)) #Probably a more efficient way to do this to look at a non-varying environment but I couldnt figure it out
#A2[,,1] = A.inter
#A2[,,2] = A.inter

# Define Parameter Values SA
df.ST = 0.3
df.SA = 0.2
dg = 1
c = 1
N0 = c(2,2)
tf = 10

par(mfrow=c(2,1), mai = c(.5,1,.5,1))
##### PLOTS FOR PRESENTATION
#good years too frequent ST
eq <- Annual(tf = 100, N0 = N0, p = 0.6, q = 0.4, df = df.ST, dg = dg, c = c, A = A.ST)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Stress Tolerant Forb, Good Year Probability = 0.6")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n",y.intersp = 1.5)

#bad years too frequent ST
eq <- Annual(tf = 100, N0 = N0, p = 0.4, q = 0.6, df = df.ST, dg = dg, c = c, A = A.ST)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Stress Tolerant Forb, Good Year Probability = 0.4")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n",y.intersp = 1.5)

#good years too frequent SA
eq <- Annual(tf = 100, N0 = N0, p = 0.6, q = 0.4, df = df.SA, dg = dg, c = c, A = A.SA)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Stress Avoiding Forb, Good Year Probability = 0.6")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n",y.intersp = 1.5)

#bad years too frequent SA
eq <- Annual(tf = 100, N0 = N0, p = 0.4, q = 0.6, df = df.SA, dg = dg, c = c, A = A.SA)
# plot
matplot(eq[[2]],type="l",ylim = c(0,35), xlab = "Years", ylab = "Density", main="Stress Avoiding Forb, Good Year Probability = 0.4")
points(eq[[1]]+33, pch = ifelse(eq[[1]]+33 > 34, "+", "-"))
legend(0,33, legend = c("Grass","Forb"), col = c("red","black"), lty = c(2,1),lwd = 1, bty = "n",y.intersp = 1.5)

#############

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


# equal probability of good and bad years (germ and ungerm fraction); high dormancy forb
eq <- Annual(tf = 100, N0 = N0, p = 0.2, q = 0.8, df = df, dg = dg, afg = afg, agf = agf, c = c, A = A.hd)
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