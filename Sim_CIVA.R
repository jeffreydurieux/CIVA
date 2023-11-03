# Thu Oct 26 14:15:37 2023 ------------------------------

# What: function to simulate CIVA data

#V = 1000; Ti = 2; Nr = 1; Q = 2; R = 2; E = 0.01; D = 3
Simulate_CIVA <- function(Nr = 40, V = 1000, Ti = 100, Q = 2, R = 2, D = 3, E = 0.1){
  ###
  # N = subjects
  # V = voxels
  # Ti = dimensions
  # Q = cluster specific components
  # R = clusters
  # E = Gaussian noise
  
  # notes:
  # datasets D are fixed at three for now
  
  # generate cluster specific components (FC patterns)
  # use rmvl from laplacedemon package
  # double exponential = b
  # Dimensions: R lists of V X Q matrices
  
  
  ### data structure:
  # Sr are multiple arrays of size Q x V x D
  # Ai are multiple arrays of size Ti x Q x D
  # X are multiple arrays of size Ti x V x D
  # Xe are multiple arrays of size Ti x V x D with scaled gaussian noise
  
  #### helper functions
  ### helper functions for generating gaussian
  addError<-function(datablock,error)
  {
    
    errorM<-replicate(ncol(datablock),rnorm(nrow(datablock)))
    errorM<-SSequal(errorM,datablock)
    errorlevel<-error/(1-error)
    
    res<-datablock + (errorM * sqrt(errorlevel))
    return(res)
  }
  
  SSequal<-function(m1,m2)
  {
    res<-(m1/sqrt(sum(m1^2)) * sqrt(sum(m2^2))) #R
    return(res)
  }
  
  
  
  ##### end helper functions
  
  
  Qr <- list()
  for(i in 1:R){
    
    Stemp <- array(NA, c(Q, V, D))
    
    for (j in 1:Q) {
      Stemp[j, , ] <- LaplacesDemon::rmvl(V, rep(0, D), diag(D))
    }
    Qr[[i]] <- Stemp
  }
  
  Ais <- list()
  for(i in 1:R){
    Atemp <- list()
    for (j in 1:Nr) {
      Atemp[[j]] <- array(rnorm(Ti * Q * D), c(Ti, Q, D))
    }
    Ais[[i]] <- Atemp
  }
  
  Xlist <- list()
  for(i in 1:R){
     
    Xjlist <- list()
    for(j in 1:Nr){
      Xtemp <- array(NaN, c(Ti, V, D))  
      for (d in 1:D) {
        Xtemp[, , d] <- Ais[[i]][[j]][, , d] %*% Qr[[i]][ , , d]
      }  
     Xjlist[[j]] <- Xtemp
    }
    Xlist[[i]] <- Xjlist
  }
 
  X <- unlist(Xlist, recursive = F)
  
  Xe <- list()
  for(i in 1:length(X)){
    Xtemp <- X[[i]]
    for(d in 1:D){
      Xtemp[ , ,d] <- addError(Xtemp[ , ,d], error = E)
    }
   Xe[[i]] <- Xtemp
  }
  
  Out <- list() 
  Out$P <- rep(1:R, each = Nr)
  Out$Qr <- Qr
  Out$Ais <- Ais
  Out$X <- X
  Out$Xe <- Xe
  return(Out)
}

# set.seed(2407)
# testdata <- Simulate_CIVA(Nr = 10, V = 1000, Ti = 2, Q = 2,R = 2,D = 3, E = .001)
# 
# X <- abind::abind(testdata$Xe[1:2], along = 1)
# A <- abind::abind(testdata$Ais[[1]][1:2], along = 1)
# Al <- testdata$Ais[[1]][1:2]
# 
# X <- testdata$Xe[[1]]
# set.seed(2407)
# tmp <- proc.time()
# res <- fastIVA(X, source_density = 'laplace_diag')
# time <- proc.time() - tmp
# time
# tt <- coef(res, which.dataset = 1)
# avg_ISI(tt,testdata$Ais[[1]][[1]])
# 
