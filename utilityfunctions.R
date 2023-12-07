

####################
# Required functions
####################

# Function to draw samples from a polya gamma function

vectorizedRpg <- function(h,z){ ## Vectorized function to facilitate samples 
  #from ploya Gamma
  omega <- NULL
  for (i in 1:length(h)) {
    omega <- c(omega,ifelse(h[i]>0,rpg.devroye(1,h[i],z[i]),0))
  }
  return(diag(omega))
}

# Function to draw posterior samples

SamplePostBeta <- function(Beta,Phi,phi.tau){  ## Posterior sample function
  #Beta -> Fixed Effect
  #Phi -> Spatial effect
  #phi.tau -> spatial precision
  
  p = dim(Beta)[1]/2
  q = dim(Phi)[1]/2
  
  A=matrix(0,2*p,2*p)    # Auxiliary matrix to compute Sum of matrices
  b=matrix(0,K-1,1)
  C=matrix(0,2*p,1)
  
  for (i in 1:n) {
    
    X  <- matrix(0,(K-1)*p,(K-1)) ## Splited Design matrix for fixed
    Z  <- matrix(0,(K-1)*r,(K-1)) ## Splited Design matrix for random
    
    
    N1 <- sum(Y[i,])    ## Obtain N1
    N2 <-  N-cumsum(Y[i,])[1] ## N2 for the ith observation
    
    X[1:p,1]     <-XX[i,]   ## Extract the corresponding design mat
    X[(1:p)+p,2] <-XX[i,]
    
    Z[id.phi[i],1] <- 1
    Z[(id.phi[i]+37),2] <- 1
    
    eta   <- t(X)%*%Beta +t(Z)%*% Phi  ## Vector of predictor
    
    
    Omega <- vectorizedRpg(c(N1,N2),eta) ## Sample Polya gamma latent variable vectors
    
    ### Linear Effect
    A1 <- X%*%Omega%*%t(X)
    b1 <- (matrix(Y[i,1:2])-matrix(c(N1,N2))/2) # Posterior quantities
    C1 <- X%*%b1-X%*%Omega%*%t(Z)%*%Phi
    
    A <- A+A1  # Sum 
    C <- C+C1 
    
  }             # Exit loop 
  #######
  ####### Posterior draws for fixed effect
  
  PostSigmaBeta <- solve(A+ solve(Sigma))
  PostMeanBeta  <- PostSigmaBeta%*%C
  post.beta     <- t(rmvnorm(1,PostMeanBeta,PostSigmaBeta)) # sample from mvnorm
  
  
  ###### Posterior draws for spatial Effect
  
  phi.A = matrix(0,2*q,2*q)
  phi.b = matrix(0,K-1,1)
  phi.C = matrix(0,2*q,1)
  
  for (i in 1:n) {
    
    X  <- matrix(0,(K-1)*p,(K-1)) ## Splited Design matrix for fixed effect
    Z  <- matrix(0,(K-1)*r,(K-1)) ## Splited Design matrix for spatial efect
    
    ######################## The following chunk should be adjusted for a different K
    # For K=3. 
    N1 <- sum(Y[i,])           ## Obtain N1
    N2 <-  N-cumsum(Y[i,])[1]  ## N2 for the ith observation
    #N3 <-  N-cumsum(Y[i,])[2]  ## N3 for the ith observation
    #N4 <-  N-cumsum(Y[i,])[3]  ## N4 and so on and so forth
    ########################
    X[1:p,1]     <-XX[i,]      ## Extract the corresponding design matrix
    X[(1:p)+p,2] <-XX[i,]
    
    Z[id.phi[i],1] <- 1
    Z[(id.phi[i]+q),2] <- 1
    
    eta   <- t(X)%*%Beta + t(Z)%*% Phi  ## Compute Vector of predictor
    
    
    Omega <- vectorizedRpg(c(N1,N2),eta) ## Sample Polya gamma latent vector
    
    
    phi.A1 <- Z%*%Omega%*%t(Z)
    phi.b1 <- (matrix(Y[i,1:2])-matrix(c(N1,N2))/2) # Posterior quantities
    phi.C1 <- Z%*%phi.b1-Z%*%Omega%*%t(X)%*%post.beta
    
    phi.A <- phi.A+phi.A1  # Sum 
    phi.C <- phi.C+phi.C1 
    
  }      
  
  PostSigmaPhi <- solve(phi.A+ phi.tau*Q.phi)
  PostMeanPhi  <- PostSigmaPhi%*%phi.C
  post.phi     <- t(rmvnorm(1,PostMeanPhi,PostSigmaPhi)) # sample from mvnorm
  phi.tau      <- rgamma(1,a1+q,bb1+0.5*t(post.phi)%*%Q.phi%*%post.phi)
  
  return(list(Beta=(post.beta),
              Phi=(post.phi),
              phi.tau=phi.tau
              
  ))
}


