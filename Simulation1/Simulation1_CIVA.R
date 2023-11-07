# Mon Nov  6 08:43:06 2023 ------------------------------
# JD

# What: simulation script for CIVA
# simulation run on ALICE

#change path
ext <- '/exports/fsw/durieuxj/P4/Sim1/'  

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

#### Design ####

Q <- c(2,4,6)
R <- c(2, 4)
Nr <- 20 #c(50, 75)
Err <- c(.1, .3, .60)
Ti <- 50
rep <- 1:20

grid <- expand.grid(R=R,Q=Q, Nr=Nr, Err=Err, rep = rep, Ti = Ti)

args <- commandArgs(TRUE)
#args <- as.numeric(args)

splits <- split(1:360, ceiling(seq_along(1:360)/5)) # list of 72

sp <- args[1]

rows <- splits[[sp]]

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




#### main sim code here

grid <- grid[rows,]

for(sim in 1:nrow(grid)){
  
  #cat('\n')
  #cat('\n')
  #cat('------------')
  #cat(sim)
  #cat('------------')
  #cat('\n')
  #cat('\n')
  seed <- as.numeric(rownames(grid))[sim]
  set.seed(seed)
  
  ##### simulate #####
  simdata <- Simulate_CIVA(Nr = grid$Nr[sim], V = 2000, Ti = grid$Ti[sim], 
                            Q = grid$Q[sim],R = grid$R[sim], D = 3, E = grid$Err[sim])
 
  X <- simdata$Xe
  X <- lapply(X, xscale) 
  
  ##### analyse #######
  
  ptm <- proc.time()
  civa <- CIVA(X = X, nc = grid$Q[sim], R = grid$R[sim], starts = 30, 
               parallel = FALSE,verbose = FALSE)
  time <- proc.time() - ptm
  #time 
  trueP <- simdata$P
  civatrue <- CIVA(X = X, nc = grid$Q[sim], R = grid$R[sim], starts = 1, 
                   parallel = FALSE,verbose = FALSE, initP = trueP)
  

  Pests <- sapply(seq_along(civa), function(x) civa[[x]]$P)
  
  ARI_starts <- apply(Pests, MARGIN = 2, adjustedRandIndex, y = simdata$P)
  #ARI_starts
  #cat(length( which(ARI_starts==1) ) /30 * 100,'%')

  losses <- sapply(seq_along(civa), function(x) civa[[x]]$LossStart)

  Optimal <- which.min(losses)

  OptStart <- civa[[Optimal]]

  ARIopt <- adjustedRandIndex(simdata$P, OptStart$P)

  #####Tucker S#####
  p_per <- FindOptimalClusPermut(OptStart$P, simdata$P)

  if(grid[sim,]$R == 2){
    # recov Q1
    SresTuck <- mean(
      mean(c(
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[1] ]][,,1]) , t(OptStart$IVA$Qr[[1]][,,1]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[1] ]][,,2]) , t(OptStart$IVA$Qr[[1]][,,2]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[1] ]][,,3]) , t(OptStart$IVA$Qr[[1]][,,3]) )$BestRecov)
      )
      ,
      # recov Q2
      mean(c(
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[2] ]][,,1]) , t(OptStart$IVA$Qr[[2]][,,1]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[2] ]][,,2]) , t(OptStart$IVA$Qr[[2]][,,2]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[2] ]][,,3]) , t(OptStart$IVA$Qr[[2]][,,3]) )$BestRecov)
      )
    )
    
    Aest1 <- inv_W(OptStart$IVA$Wr[[1]])
    Aest1 <- Asplit(Aest1, dimsplit = grid$Ti[sim])

    Aest2 <- inv_W(OptStart$IVA$Wr[[2]])
    Aest2 <- Asplit(Aest2, dimsplit = grid$Ti[sim])

    Ais1 <- numeric()
    for(i in 1:length(simdata$Ais[[ p_per$BestPerm[1]]]) ){
    Ais1[i] <- mean(c(
      FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,1] , Aest1[[i]][,,1] )$BestRecov,
      FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,2] , Aest1[[i]][,,2] )$BestRecov,
      FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,3] , Aest1[[i]][,,3] )$BestRecov)
    )
    }

    Ais2 <- numeric()
    for(i in 1:length(simdata$Ais[[ p_per$BestPerm[2]]]) ){
      Ais2[i] <- mean(c(
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,1] , Aest2[[i]][,,1] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,2] , Aest2[[i]][,,2] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,3] , Aest2[[i]][,,3] )$BestRecov)
      )
    }
    # Res A
    AresTuck <- mean(mean(Ais1), mean(Ais2))
    
  }else{
    SresTuck <- mean(
      mean(c(
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[1] ]][,,1]) , t(OptStart$IVA$Qr[[1]][,,1]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[1] ]][,,2]) , t(OptStart$IVA$Qr[[1]][,,2]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[1] ]][,,3]) , t(OptStart$IVA$Qr[[1]][,,3]) )$BestRecov)
      )
      ,
      # recov Q2
      mean(c(
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[2] ]][,,1]) , t(OptStart$IVA$Qr[[2]][,,1]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[2] ]][,,2]) , t(OptStart$IVA$Qr[[2]][,,2]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[2] ]][,,3]) , t(OptStart$IVA$Qr[[2]][,,3]) )$BestRecov)
      )
      ,
      # recov Q3
      mean(c(
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[3] ]][,,1]) , t(OptStart$IVA$Qr[[2]][,,1]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[3] ]][,,2]) , t(OptStart$IVA$Qr[[2]][,,2]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[3] ]][,,3]) , t(OptStart$IVA$Qr[[2]][,,3]) )$BestRecov)
      )
      ,
      # recov Q2
      mean(c(
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[4] ]][,,1]) , t(OptStart$IVA$Qr[[4]][,,1]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[4] ]][,,2]) , t(OptStart$IVA$Qr[[4]][,,2]) )$BestRecov,
        FindOptimalPermutSingle( t(simdata$Qr[[ p_per$BestPerm[4] ]][,,3]) , t(OptStart$IVA$Qr[[4]][,,3]) )$BestRecov)
      )
    )
    
    Aest1 <- inv_W(OptStart$IVA$Wr[[1]])
    Aest1 <- Asplit(Aest1, dimsplit = grid$Ti[sim])
    
    Aest2 <- inv_W(OptStart$IVA$Wr[[2]])
    Aest2 <- Asplit(Aest2, dimsplit = grid$Ti[sim])
    
    Aest3 <- inv_W(OptStart$IVA$Wr[[3]])
    Aest3 <- Asplit(Aest3, dimsplit = grid$Ti[sim])
    
    Aest4 <- inv_W(OptStart$IVA$Wr[[4]])
    Aest4 <- Asplit(Aest4, dimsplit = grid$Ti[sim])
    
    Ais1 <- numeric()
    for(i in 1:length(simdata$Ais[[ p_per$BestPerm[1]]]) ){
      Ais1[i] <- mean(c(
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,1] , Aest1[[i]][,,1] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,2] , Aest1[[i]][,,2] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,3] , Aest1[[i]][,,3] )$BestRecov)
      )
    }
    
    Ais2 <- numeric()
    for(i in 1:length(simdata$Ais[[ p_per$BestPerm[2]]]) ){
      Ais2[i] <- mean(c(
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,1] , Aest2[[i]][,,1] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,2] , Aest2[[i]][,,2] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,3] , Aest2[[i]][,,3] )$BestRecov)
      )
    }
    
    Ais3 <- numeric()
    for(i in 1:length(simdata$Ais[[ p_per$BestPerm[3]]]) ){
      Ais3[i] <- mean(c(
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[3] ]][[i]][,,1] , Aest3[[i]][,,1] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[3] ]][[i]][,,2] , Aest3[[i]][,,2] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[3] ]][[i]][,,3] , Aest3[[i]][,,3] )$BestRecov)
      )
    }
    
    Ais4 <- numeric()
    for(i in 1:length(simdata$Ais[[ p_per$BestPerm[4]]]) ){
      Ais4[i] <- mean(c(
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[4] ]][[i]][,,1] , Aest4[[i]][,,1] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[4] ]][[i]][,,2] , Aest4[[i]][,,2] )$BestRecov,
        FindOptimalPermutSingle( simdata$Ais[[ p_per$BestPerm[4] ]][[i]][,,3] , Aest4[[i]][,,3] )$BestRecov)
      )
    }
    # Res A
    AresTuck <- mean(mean(Ais1), mean(Ais2), mean(Ais3), mean(Ais4))
  }
  
  
  #### output list object #####
  output <- list()
  output$id <- grid[sim,]
  output$seed <- seed
  output$simdata <- simdata
  output$ari <- ARIopt
  output$ariall <- ARI_starts
  output$TuckS <- SresTuck
  output$TuckA <- AresTuck
  output$civaTrueP <- civatrue
  output$time <- time
  
  ext <- '/exports/fsw/durieuxj/P4/Sim1/'  
  ext <- paste(ext,'CIVA_sim1_',rows[sim], '.Rdata',sep = '')
  
  save(output,file = ext)
  
}