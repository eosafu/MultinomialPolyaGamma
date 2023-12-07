# Code to perform simulation
# including fixed and spatial effect.

library(mvtnorm)
library(magic)
library(BayesLogit)
sorce("utilityfunctions.R")

#%%%%%%%%%%%%%%%%%%%%
# Simulation begins
#%%%%%%%%%%%%%%%%%%%
# Sample syntetic data
p=3   # For fixed effect Including intercept
K=3   # Number of buckets
N=6 # Total number of items

### Nigeria map adjacency matrix to build spatial model

A <- read.csv("https://raw.githubusercontent.com/eosafu/Model/master/Nig_adj_Matrix.txt",sep=",",header = F) # Nigeria adjacency matrix
A = A[,-1]

#### Get Precision matrix
m<-apply(A,1,sum)	
kappa<-.99	  	            # Spatial dependence parameter ~ 1 for intrinsic CAR
Q<- (diag(m))-kappa*(as.matrix(A))

## Fix spatial effect 
tau.phi<-tau.phiTrue <- 1
Q.phi  <- adiag(tau.phi*Q,1*Q)

m.phi.1  <- sample(seq(0.01,1,length.out=200), dim(Q)[1],replace = TRUE) 
m.phi.2  <- sample(seq(-1,-0.01,length.out=200), dim(Q)[1],replace = TRUE) 


# Number of sample sizes
SampleSizes = 200

# Replication within sample sizes to calculate adequacy measures (MSE, MRE and CPO)
numbIT=30

# Cpo storage initialization
cpo    = matrix(0,length(SampleSizes),numbIT)


for (kk in 1:length(SampleSizes)) {
  
  set.seed(1456)
  
  n =  SampleSizes[kk]
  FixRE  =matrix(0,(K-1)*p,numbIT)
  RNDRE  =matrix(0,(K-1)*q,numbIT)
  RNtauRE= vector("numeric",numbIT)
  
  #Initialize matrix for Mean Square Error
  FixMSE  =matrix(0,(K-1)*p,numbIT)
  RNDMSE  =matrix(0,(K-1)*q,numbIT)
  RNtauMSE=vector("numeric",numbIT)

  
  ## Sample spatial covariate
  q = dim(Q)[1] # Size of spatial effect
  # Fix true parameters
  
  ## Fix true betas
  
  beta10=  2.5
  beta11=  -1.5
  beta12=  -0.5
  
  beta20=  -0.5
  beta21=  2.0
  beta22=  -1.0
  Beta    <- BetaTrue <- matrix(c(beta10,beta11,beta12,beta20,beta21,beta22))
  
  # Fix true spatial effects

  PhiTrue    <- rmvnorm(1,c(m.phi.1,m.phi.2), solve(Q.phi))
  phi.tru <- PhiTrue 
  tau.phi <- tau.phiTrue <- 1 
  
  # Get spatial covariate
  id.phi  <- c(c(1:q),sample(1:q,n-q,replace =T))
  
  # For loop to generate synthetic data
  
  for (index in 1:numbIT) {
    
    x0 <- rep(1,n)       # covariates
    x1 <- rnorm(n, 0,1)
    x2 <- rnorm(n, 0,1)
    
    XX <- as.matrix(cbind(x0,x1,x2)) #Design matrix
    
    
    Y <-  matrix(0, n,K)          ## To hold outcome samples
    
    for (i in 1:n) {
      
      X  <- matrix(0,(K-1)*p,(K-1)) ## Splited Design matrix for fixed
      Z  <- matrix(0,(K-1)*r,(K-1)) ## Splited Design matrix for random
      
      X[1:p,1]     <- XX[i,]        # o obtain X_i
      X[(1:p)+p,2] <- XX[i,]        # covariate are the same for all the classes
      
      Z[id.phi[i],1] <- 1
      Z[(id.phi[i]+q),2] <- 1
      
      eta <-  t(X)%*%Beta + t(Z)%*% matrix(phi.tru[1,])    # Linear Predictor
      
      y <- Tild_pi <- exp(eta)/(1+exp(eta))                # Logit inverse to calculate \pi tilda
      
      y[1] <- rbinom(1,N,Tild_pi[1])
      y[2] <- rbinom(1,N-y[1],Tild_pi[2])
      y <- c(y[1],y[2],N-(y[1]+y[2]))
      
      # Response matrix
      Y[i,] <- y
      
    }

#%%%%%%%%%%%%%%%%%
# Estimation step
#%%%%%%%%%%%%%%%%

# Hyperparameters
    a1 <- 5
    bb1 <- 2
    Sigma <- 1000*diag(rep(1,(K-1)*p))
# Initialization
    Beta     <- matrix(0,(K-1)*p,1)
    Phi      <- matrix(0,(K-1)*q,1)
    
    phi.tau  <- 1
    Num      <- 3000
    burning <- 1000
    BETA     <- matrix(0,(K-1)*p,Num)
    PHI      <- matrix(0,(K-1)*q,Num)
    PHI.tau  <-vector("numeric",Num)
    
#Log likelihood storage
    
    LIK <-  matrix(0,n,Num)
    
#%%%%%%%%%%%%%%%%%
## Estimation begins
#%%%%%%%%%%%%%%%%%
     for (k in 1:Num) {
       # Update model parameters
      Aux    <-  SamplePostBeta(Beta = Beta,Phi = Phi, phi.tau = phi.tau)
      
   # Extract samples from the posterior distribution
      
      BETA[,k]  <-  Beta  <- Aux[[1]]
      PHI[,k] <- Phi <- Aux[[2]]
      Phi[1:q,1]  <- Phi[1:q,1] -mean(Phi[1:q,1])
      Phi[(q+1):(2*q),1] <- Phi[(q+1):(2*q),1]-mean(Phi[(q+1):(2*q),1])
      PHI.tau[k]  <- phi.tau  <- Aux[[3]]
      
  # calculate likelihood
      
      Likelihood = NULL
      for (i in 1:n) {
        
        X  <- matrix(0,(K-1)*p,(K-1)) ## Splited Design matrix for fixed
        Z  <- matrix(0,(K-1)*r,(K-1)) ## Splited Design matrix for random
        
        X[1:p,1]     <- XX[i,]   # o obtain X_i
        X[(1:p)+p,2] <- XX[i,]  # covariate are the same for all the classes
        
        Z[id.phi[i],1] <- 1
        Z[(id.phi[i]+37),2] <- 1
        
        eta <-  t(X)%*%Beta + t(Z)%*% Phi #matrix(phi.tru[1,])    # Linear Predictor
        
        y <- Tild_pi <- exp(eta)/(1+exp(eta)) # Logit inverse to calculate \pi tilda
        
        y[1] <- rbinom(1,N,Tild_pi[1])
        y[2] <- rbinom(1,N-y[1],Tild_pi[2])
        y <- c(y[1],y[2],N-(y[1]+y[2]))
        t_p1 <- Tild_pi[1]               # Compute \pi 
        t_p2 <- Tild_pi[2]*(1-t_p1)
        t_p3 <- 1-t_p1-t_p2
        Likelihood <- c(Likelihood,dmultinom(Y[i,],N,c(t_p1,t_p2,t_p3)))
        
      }
      LIK[,k] <- Likelihood
      if(k%%500==0){
        cat(k,"\n")
      }
      
    }
    CPO <-  mean(log(1/rowMeans(1/LIK))[sample(1:n,200,replace = T)])
   
     ### Model Adequacy

    #Relative error
    cpo[kk,index] <- CPO
    FixRE[,index]=rowMeans(BETA[,burning:Num ])/BetaTrue
    RNDRE[,index]=rowMeans(PHI[,burning:Num ])/PhiTrue 
    RNtauRE[index]=mean(PHI.tau[burning:Num ])/tau.phiTrue 
    #Mean Square Error
    FixMSE[,index]=(rowMeans(BETA[,burning:Num ])-BetaTrue)^2
    RNDMSE[,index]=(rowMeans(PHI[,burning:Num ])-PhiTrue )^2
    RNtauMSE[index]=(mean(PHI.tau[burning:Num ])-tau.phiTrue)^2

  }
  propRE=rbind(FixRE,RNDRE,RNtauRE)
  propMSE=rbind(FixMSE,RNDMSE,RNtauMSE)
 
}


#### spatial plot

est.comp1 <- rowMeans(PHI[,burning:3000])[1:q]
est.comp2 <- rowMeans(PHI[,burning:3000])[(q+1):(2*q)]
true.comp1 <- phi.tru[1:q]
true.comp2 <- phi.tru[(q+1):(2*q)]

data.frame(true.comp1,est.comp1,true.comp2,round(est.comp2,3))

###########
############
require(rgdal)
require(ggplot2)
fn <- file.path(getwd(), "nigeria-lgas.zip", fsep = "\\")
utils::unzip(fn, exdir = getwd())
shp <- read_sf(dsn ="new_lga_nigeria_2003.shp")

shp.map <- fortify(shp, region = "STATE")

shp.map$estimate=0

ss  <- read.csv("https://raw.githubusercontent.com/eosafu/Poly--GammaInference/5a230ab3e24457bfa18ff682d80d669401995f50/NIGBND.csv")
state2 <-  ss$state2
cloc2 <-   ss$cloc2

est_data <- data.frame(region=state2,est=est.comp1[cloc2])

s <- function(region,est){
  shp.map$estimate[which(shp.map$group==as.character(region))] <- est
  return(shp.map)
}

# Mean 

for (i in 1:37) {
  shp.map <-s(est_data[i,1],est_data[i,2]) 
}
shp.map$estimate[which(shp.map$group=="Lagos.2")] <- est_data[25,2]
shp.map$estimate[which(shp.map$group=="Lake.1")] <- est_data[25,2]

ggplot(shp.map, aes(x = long, y = lat, group = group, fill = count)) +
  geom_polygon(colour = "black", size = 0.5, aes(group = group)) +
  theme_minimal()+ scale_fill_gradient(name="Risk",high='green', low='red')
