nyear<-100  #number of years for simulation
rho<- 0  #temporal autocorrelation
Pfav<-0.5  #probability of a favorable year
gFu<-.1  #germination (g) of the forb (F) in an unfavorable year (u)
gFf<-.7  #germination (g) of the forb (F) in a favorable year (f)
gGu<-.9  #germination (g) of the grass (G) in an unfavorable year (u)
gGf<-.9  #germination (g) of the grass (G) in a favorable year (f)
dF<-.1  #annual death rate (d) of forb (F) seeds in the soil
dG<-.7  #annual death rate (d) of grass (G) seeds in the soil
lambda.u<-5  #fecundity for both species in an unfavorable year (u) 
lambda.f<-30  #fecundity for both species in a favorable year (f)
alphaFG<-2  #competition coefficient for grass (G) effects on forbs (F)
alphaGF<-.5  #competition coefficient for forb (f) effects on grass (G) 
c<-1  #constant related to how competition affects density

q<-c()  # Vector of 1s and 2s randomly assigned based on Pfav and rho (1 is a unfav year, 2 is a fav year)
q[1]<-1+rbinom(1,1,Pfav)
for (y in  2:nyear){
  if (q[y-1]==1) q[y]<-1+rbinom(1,1,(Pfav*(1-rho))) else q[y]<-1+rbinom(1,1,(1-(1-Pfav)*(1-rho)))
}


S<-matrix(NA, nrow=nyear, ncol=2)  #matrix of (seed) densities where the forb is in column 1 and the
#grass is in column 2; each row is a year
S[1,]<-2  #begins the simulation with two individuals of each species in year 1
for (t in 2:nyear){  #loops through the model applying different parameters if it is a favorable or 
  #unfavorable year, following equations (2) and (3) of the paper
  
  if(q[t-1]==1) S[t,1]<- (1-gFu)*(1-dF)*S[t-1,1]+gFu*lambda.u*S[t-1,1]/(c+ gFu*S[t-1,1]+ alphaFG*gGu*S[t-1,2])
  if(q[t-1]==1) S[t,2]<- (1-gGu)*(1-dG)*S[t-1,2]+gGu*lambda.u*S[t-1,2]/(c+ gGu*S[t-1,2]+ alphaGF*gFu*S[t-1,1])
  
  if(q[t-1]==2) S[t,1]<- (1-gFf)*(1-dF)*S[t-1,1]+gFf*lambda.f*S[t-1,1]/(c+ gFf*S[t-1,1]+ alphaFG*gGf*S[t-1,2])
  if(q[t-1]==2) S[t,2]<- (1-gGf)*(1-dG)*S[t-1,2]+gGf*lambda.f*S[t-1,2]/(c+ gGf*S[t-1,2]+ alphaGF*gFf*S[t-1,1])
  
}
par(mfrow=c(2,1))
matplot(1:nyear, S, col=1, type="l")  #Plots grass and forb dynamics through time