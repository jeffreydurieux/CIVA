# Fri Nov 10 14:26:23 2023 ------------------------------
# JD  
# CIVA sim 3 

# use two clusters

# both clusters have med visual and DMN intact in ses1

# one cluster have intact med visual and DMN intact in ses2
# one cluster have intact med visual and slightly affected DMN in ses2

# one cluster have intact med visual and DMN intact in ses3
# one cluster have intact med visual and severe affected DMN in ses2

# hypo: clustering will be driven completly by disrupted DMN

setwd('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/')
library(mclust)
library(gtools)
library(ivaBSS)
library(multiway)
library(abind)
library(NMFN)
source('EstAih.R')
source('fastIVA_custom.R')
source('Sim_CIVA.R')
source('whiten_custom.R')
source('clusf.R')
source('invW_Asplit.R')
source('ConcData.R')
source('ExtractIVA.R')
source('Lir.R')
source('Reclus.R')
source('Xscale.R')
source('SearchEmptyClusters.R')
source('CIVA.R')

### helper functions #####
FindOptimalPermutSingle <- function( Sest , Strue, verbose = FALSE)
{
  # code to search the optimal permutation of estimated ICA components for
  # comparing it with simulated components
  # Author(s): Tom F. Wilderjans and minor adjustments by Jeffrey Durieux
  
  #JD: code from Tom, adjusted for matrix vs matrix comparison
  #Sest, Strue (nVoxels x nSources)
  
  #library(gtools)
  N_sources = dim(Sest)[2]
  
  
  AllPerms = permutations( n = N_sources , r = N_sources , v = 1:N_sources )
  nPerms = dim(AllPerms)[1]
  
  #Find best permutation
  BestRecov = -9999
  BestPerm = -9999
  for( permtel in 1:nPerms )
  {
    if(verbose == TRUE)
    {
      if( (permtel%%50) == 0)
      {
        print( paste( "perm: " , permtel , "/" , nPerms ) )
      }
    }
    
    
    tp = AllPerms[permtel,]
    tempRecovBlock = matrix( -9999 , 1 , 1 )
    
    tempRecovBlock[1] = mean( abs( diag( congru(Strue ,
                                                Sest[, tp] ) ) ) )
    # niet nodig als het goed is
    tempRecov = mean(tempRecovBlock)
    
    if( permtel==1 )
    {
      BestRecov = tempRecov
      BestRecovBlock = tempRecovBlock
      BestPerm = tp
    }
    else
    {
      if( (tempRecov-BestRecov)>.0000000001 )
      {
        BestRecov = tempRecov
        BestRecovBlock = tempRecovBlock
        BestPerm = tp
      }
    }
    rm(tp,tempRecov,tempRecovBlock)
  }
  Out = list()
  Out$BestRecov = BestRecov
  Out$BestRecovBlock = BestRecovBlock
  Out$BestPerm = BestPerm
  Out$TuckerMatrix = congru(Strue , Sest[, BestPerm] )
  return(Out)
}

FindOptimalClusPermut <- function(Pest, Ptrue){
  # find optimal cluster permutation of estimated clustering
  # compared to simulated clustering
  # Author(s): Tom F. Wilderjans and minor adjustments by Jeffrey Durieux
  clus <- length(unique(Pest))
  
  AllPerms = gtools::permutations( n = clus , r = clus)
  nPerms = dim(AllPerms)[1]
  
  BestRecov = -9999
  BestPerm = -9999
  for( permtel in 1:nPerms )
  {
    if( (permtel%%50) == 0)
    {
      print( paste( "perm: " , permtel , "/" , nPerms ) )
    }
    
    tp = AllPerms[permtel,]
    tempRecovBlock = matrix( -9999 , 1 , 1 )
    
    tab <- table(Ptrue, Pest)
    
    tempRecovBlock[1] = sum( diag( tab[,tp] ) )
    
    tempRecov = mean(tempRecovBlock)
    
    if( permtel==1 )
    {
      BestRecov = tempRecov
      BestRecovBlock = tempRecovBlock
      BestPerm = tp
    }
    else
    {
      if( (tempRecov-BestRecov)>.0000000001 )
      {
        BestRecov = tempRecov
        BestRecovBlock = tempRecovBlock
        BestPerm = tp
      }
    }
    rm(tp,tempRecov,tempRecovBlock)
  }
  
  Out = list()
  Out$BestRecov = BestRecov
  Out$BestRecovBlock = BestRecovBlock
  Out$BestPerm = BestPerm
  return(Out)
  
}

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


Nr <- 20 #c(50, 75)
E <- c(.1)#, .3, .60)
Ti <- 50
rep <- 1:20

grid <- expand.grid(Nr=Nr, E=E, rep = rep, Ti = Ti)
setwd('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/Simulation3/Maps/')

files <- dir()
files <- mixedsort(files)
files

id <- 1
ffs1 <- paste('Maps_S1_',id,'.Rdata', sep = '')
ffs2 <- paste('Maps_S2_',id,'.Rdata', sep = '')

ids1 <- which(files==ffs1)
ids2 <- which(files==ffs2)

Qr <- list()
Qr[[1]] <- get(load(files[ids1]))
Qr[[2]] <- get(load(files[ids2]))

Ais <- list()
for(i in 1:2){
  Atemp <- list()
  for (j in 1:Nr) {
    Atemp[[j]] <- array(rnorm(Ti * 2 * 3), c(Ti, 2, 3))
  }
  Ais[[i]] <- Atemp
}

Xlist <- list()
for(i in 1:2){
  
  Xjlist <- list()
  for(j in 1:Nr){
    Xtemp <- array(NaN, c(Ti, 2553, 3))  
    for (d in 1:3) {
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
  for(d in 1:3){
    Xtemp[ , ,d] <- addError(Xtemp[ , ,d], error = E)
  }
  Xe[[i]] <- Xtemp
}


X <- Xe
X <- lapply(X, xscale) 


ptm <- proc.time()
civa <- CIVA(X = X, nc = 2, R = 2, starts = 5, 
             parallel = FALSE,verbose = TRUE)
time <- proc.time() - ptm
time

Pests <- sapply(seq_along(civa), function(x) civa[[x]]$P)

simP <- c(rep(1,20), rep(2,20))
ARI_starts <- apply(Pests, MARGIN = 2, adjustedRandIndex, y = simP)
#ARI_starts
#cat(length( which(ARI_starts==1) ) /30 * 100,'%')

losses <- sapply(seq_along(civa), function(x) civa[[x]]$LossStart)
losses
Optimal <- which.min(losses)

OptStart <- civa[[Optimal]]
# evaluation metrics here


###### write to nifti ######
OptStart$P
per_p <- FindOptimalClusPermut(Pest = OptStart$P, Ptrue = simP)

FindOptimalPermutSingle(Sest = t(OptStart$IVA$Qr[[1]][,,1]), 
                          Strue = t(Qr[[ per_p$BestPerm[1] ]][,,1]) )
FindOptimalPermutSingle(Sest = t(OptStart$IVA$Qr[[1]][,,2]), 
                        Strue = t(Qr[[ per_p$BestPerm[1] ]][,,2]) )
FindOptimalPermutSingle(Sest = t(OptStart$IVA$Qr[[1]][,,3]), 
                        Strue = t(Qr[[ per_p$BestPerm[1] ]][,,3]) )

FindOptimalPermutSingle(Sest = t(OptStart$IVA$Qr[[2]][,,1]), 
                        Strue = t(Qr[[ per_p$BestPerm[2] ]][,,1]) )
FindOptimalPermutSingle(Sest = t(OptStart$IVA$Qr[[2]][,,2]), 
                        Strue = t(Qr[[ per_p$BestPerm[2] ]][,,2]) )
FindOptimalPermutSingle(Sest = t(OptStart$IVA$Qr[[2]][,,3]), 
                        Strue = t(Qr[[ per_p$BestPerm[2] ]][,,3]) )




# library(CICA)
# #n <- readNifti('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Data/23x_all_beckmann_8_and_csf_and_wm.nii.gz')
# mask <- readNifti('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Data/brainmask_23_28_23.nii.gz')
# mask <- matrix(mask)
# maskid <- which(mask == 1)
# 
# EmptyS1_ses1 <- matrix(data = 0, nrow = dim(mask)[1], ncol = 2)
# EmptyS1_ses2 <- matrix(data = 0, nrow = dim(mask)[1], ncol = 2)
# EmptyS1_ses3 <- matrix(data = 0, nrow = dim(mask)[1], ncol = 2)
# 
# EmptyS2_ses1 <- matrix(data = 0, nrow = dim(mask)[1], ncol = 2)
# EmptyS2_ses2 <- matrix(data = 0, nrow = dim(mask)[1], ncol = 2)
# EmptyS2_ses3 <- matrix(data = 0, nrow = dim(mask)[1], ncol = 2)
# 
# OptStart$P
# EmptyS1_ses1[maskid,] <- abs(t(OptStart$IVA$Qr[[1]][,,3]))
# #EmptyS1_ses1[maskid,] <- t(Qr[[1]][,,1])
# EmptyS2_ses1[maskid,] <- abs(t(OptStart$IVA$Qr[[2]][,,3]))
# 
# 
# x <- list()
# x$Sr[[1]] <- EmptyS1_ses1
# x$Sr[[2]] <- EmptyS2_ses1
# class(x) <- 'CICA'
# plot(x)
# 
# congru(t(OptStart$IVA$Qr[[1]][,,1]),t(OptStart$IVA$Qr[[1]][,,2]))
# congru(t(OptStart$IVA$Qr[[1]][,,1]),t(OptStart$IVA$Qr[[1]][,,3]))
# congru(t(OptStart$IVA$Qr[[1]][,,2]),t(OptStart$IVA$Qr[[1]][,,3]))
# 
# congru(t(OptStart$IVA$Qr[[2]][,,1]),t(OptStart$IVA$Qr[[2]][,,2]))
# congru(t(OptStart$IVA$Qr[[2]][,,1]),t(OptStart$IVA$Qr[[2]][,,3]))
# congru(t(OptStart$IVA$Qr[[2]][,,2]),t(OptStart$IVA$Qr[[2]][,,3]))
# 
# 
# congru(t(OptStart$IVA$Qr[[2]][,,1]),t(Qr[[1]][,,1]))
# congru(t(OptStart$IVA$Qr[[2]][,,2]),t(Qr[[1]][,,2]))
# congru(t(OptStart$IVA$Qr[[2]][,,3]),t(Qr[[1]][,,3]))
# 
# congru(t(OptStart$IVA$Qr[[1]][,,1]),t(Qr[[2]][,,1]))
# congru(t(OptStart$IVA$Qr[[1]][,,2]),t(Qr[[2]][,,2]))
# congru(t(OptStart$IVA$Qr[[1]][,,3]),t(Qr[[2]][,,3]))