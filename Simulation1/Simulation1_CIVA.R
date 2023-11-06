# Mon Nov  6 08:43:06 2023 ------------------------------
# JD

# What: simulation script for CIVA
# simulation run on ALICE

#change path
ext <- '/exports/fsw/durieuxj/P4/Sim1/'  

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

#### Design ####

Q <- c(2,4,6)
R <- c(2, 4)
Nr <- 50 #c(50, 75)
Err <- c(.1, .3, .60)
rep <- 1:20

grid <- expand.grid(Q=Q, R=R, Nr=Nr, Err=Err, rep = rep)

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
    
    tempRecovBlock[1] = mean( abs( diag( Tucker(Strue ,
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
  Out$TuckerMatrix = Tucker(Strue , Sest[, BestPerm] )
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
  testdata <- Simulate_CIVA(Nr = 5, V = 2000, Ti = Ti, Q = 2,R = 2,D = 3, E = 0.1)
  simdata <- Simulate_CJICA(Nk = grid[sim,]$N, 
                            Vm = 2500,
                            K = grid[sim, ]$R,
                            Qm = grid[sim, ]$Q,
                            Err = grid[sim, ]$Err,
                            M = 2,
                            cor = grid[sim, ]$rho,
                            con = grid[sim, ]$lap
  )
  
  ##### analyse #######
  
  ptm <- proc.time()
  cjica <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                           nc = grid[sim, ]$Q, starts = 100, scale = T)
  time <- proc.time() - ptm
  
  cjicatrue <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                               nc = grid[sim, ]$Q, starts = 1, scale = T,
                               rational = simdata$P )
  
  cjicaperbs <- list()
  for(perbs in 1:10){
    pert <- perturbation(p = simdata$P, percentage = 0.1)
    cjicaPert <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                                 nc = grid[sim, ]$Q, starts = 1, scale = T,
                                 rational = pert)
    cjicaperbs[[perbs]] <- cjicaPert[[1]]
  }
  
  ##### cjica on mod one #####
  f1 <- sqrt(5000/sum(simdata$Xe[1:2500,]^2))
  X1 <- f1*simdata$Xe[1:2500,]
  cjica_m1 <- ClusterwiseJICA(X = X1, k = grid[sim,]$R,
                              nc = grid[sim, ]$Q, starts = 100, scale = F)
  
  ##### cjica on mod two #####
  f2 <- sqrt(5000/sum(simdata$Xe[2501:5000,]^2))
  X2 <- f2*simdata$Xe[2501:5000,]
  cjica_m2 <- ClusterwiseJICA(X = X2, k = grid[sim,]$R,
                              nc = grid[sim, ]$Q, starts = 100, scale = F)
  
  #### rational starts based on KM and HCL
  km <- kmeans(x = t(simdata$Xe), centers = grid[sim, ]$R, nstart = 100)
  
  cjica_km <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                              nc = grid[sim, ]$Q, starts = 1, scale = T,
                              rational = km$cluster)
  
  hcl <- hclust(dist( t(simdata$Xe)), method = 'ward.D2' )
  cut <- cutree(hcl, k = grid[sim,]$R)
  cjica_hcl <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                               nc = grid[sim, ]$Q, starts = 1, scale = T,
                               rational = cut )
  
  
  ##### evaluate ######
  loss100 <- sapply(seq_along(cjica), function(anom) tail(cjica[[anom]]$lossiter, n = 1))
  optimal <- cjica[[which.min(loss100)]]
  
  loss100m1 <- sapply(seq_along(cjica_m1), function(anom) tail(cjica_m1[[anom]]$lossiter, n = 1))
  optimalm1 <- cjica_m1[[which.min(loss100m1)]]
  
  loss100m2 <- sapply(seq_along(cjica_m2), function(anom) tail(cjica_m2[[anom]]$lossiter, n = 1))
  optimalm2 <- cjica_m2[[which.min(loss100m2)]]
  
  ### adjusted rand ###
  ari <- adjustedRandIndex(simdata$P, optimal$p)
  aritrueP <- adjustedRandIndex(simdata$P, cjicatrue[[1]]$p)
  arim1 <- adjustedRandIndex(simdata$P, optimalm1$p)
  arim2 <- adjustedRandIndex(simdata$P, optimalm2$p)
  
  
  ### Tucker S ###
  tucker_cor_lap <- unlist(TuckCheck(simdata$S))
  
  # add tucker congru 
  clusper <- FindOptimalClusPermut(optimal$p, simdata$P)
  
  if(grid[sim,]$R == 3){
    tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov))
    
  }else{
    tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[4]]], Strue = simdata$S[[4]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[5]]], Strue = simdata$S[[5]])$BestRecov))
    
  }
  
  #### output list object #####
  output <- list()
  output$id <- grid[sim,]
  output$seed <- seed
  output$simdata <- simdata
  output$ari <- ari
  output$S <- tucker_S
  output$Araw <- optimal$ica$Mr
  output$time <- time
  output$optimal <- optimal
  output$loss100 <- loss100
  output$tuckercheck <- tucker_cor_lap
  output$cjicatrue <- cjicatrue[[1]]
  output$cjicaperbs <- cjicaperbs
  output$m1 <- list(optimalm1, arim1)
  output$m2 <- list(optimalm2, arim2)
  output$km <- cjica_km
  output$km <- cjica_hcl
  
  ext <- '/exports/fsw/durieuxj/P4/Sim1/'  
  ext <- paste(ext,'CJICA_sim1_',rows[sim], '.Rdata',sep = '')
  save(output,file = ext)
  
}