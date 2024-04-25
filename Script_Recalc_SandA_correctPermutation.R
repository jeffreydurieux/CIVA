library(gtools)

setwd('/home/durieuxj/data1/CIVA_sim/')
files <- dir()
files <- mixedsort(files)

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


calc <- function(outputfile){
  R <- length(output$simdata$Qr)
  
  if(R == 4){
    SresTuck <- mean( c(
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[1]]][,,1]) , Strue = t(simdata$Qr[[1]][,,1]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[1]]][,,2]) , Strue = t(simdata$Qr[[1]][,,2]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[1]]][,,3]) , Strue = t(simdata$Qr[[1]][,,3]))$BestRecov,
      
      
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[2]]][,,1]) , Strue = t(simdata$Qr[[2]][,,1]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[2]]][,,2]) , Strue = t(simdata$Qr[[2]][,,2]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[2]]][,,3]) , Strue = t(simdata$Qr[[2]][,,3]))$BestRecov,
      
      
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[3]]][,,1]) , Strue = t(simdata$Qr[[3]][,,1]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[3]]][,,2]) , Strue = t(simdata$Qr[[3]][,,2]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[3]]][,,3]) , Strue = t(simdata$Qr[[3]][,,3]))$BestRecov,
      
      
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[4]]][,,1]) , Strue = t(simdata$Qr[[4]][,,1]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[4]]][,,2]) , Strue = t(simdata$Qr[[4]][,,2]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[4]]][,,3]) , Strue = t(simdata$Qr[[4]][,,3]))$BestRecov
    ))
    
    
    Aest1 <- inv_W(OptStart$IVA$Wr[[p_per$BestPerm[1]]])
    Aest1 <- Asplit(Aest1, dimsplit = grid$Ti[sim])
    
    Aest2 <- inv_W(OptStart$IVA$Wr[[p_per$BestPerm[2]]])
    Aest2 <- Asplit(Aest2, dimsplit = grid$Ti[sim])
    
    Aest3 <- inv_W(OptStart$IVA$Wr[[p_per$BestPerm[3]]])
    Aest3 <- Asplit(Aest3, dimsplit = grid$Ti[sim])
    
    Aest4 <- inv_W(OptStart$IVA$Wr[[p_per$BestPerm[4]]])
    Aest4 <- Asplit(Aest4, dimsplit = grid$Ti[sim])
    
    Ais1 <- numeric()
    for(i in 1:length(simdata$Ais[[ 1 ]]) ){
      Ais1[i] <- mean(c(
        FindOptimalPermutSingle( Aest1[[ i ]][,,1] , simdata$Ais[[ 1 ]][[i]][,,1])$BestRecov,
        FindOptimalPermutSingle( Aest1[[ i ]][,,2] , simdata$Ais[[ 1 ]][[i]][,,2])$BestRecov,
        FindOptimalPermutSingle( Aest1[[ i ]][,,3] , simdata$Ais[[ 1 ]][[i]][,,3])$BestRecov)
      )
    }
    
    Ais2 <- numeric()
    for(i in 1:length(simdata$Ais[[ 2 ]]) ){
      Ais2[i] <- mean(c(
        FindOptimalPermutSingle( Aest2[[ i ]][,,1] , simdata$Ais[[ 2 ]][[i]][,,1])$BestRecov,
        FindOptimalPermutSingle( Aest2[[ i ]][,,2] , simdata$Ais[[ 2 ]][[i]][,,2])$BestRecov,
        FindOptimalPermutSingle( Aest2[[ i ]][,,3] , simdata$Ais[[ 2 ]][[i]][,,3])$BestRecov)
      )
    }
    
    Ais3 <- numeric()
    for(i in 1:length(simdata$Ais[[ 3 ]]) ){
      Ais3[i] <- mean(c(
        FindOptimalPermutSingle( Aest3[[ i ]][,,1] , simdata$Ais[[ 3 ]][[i]][,,1])$BestRecov,
        FindOptimalPermutSingle( Aest3[[ i ]][,,2] , simdata$Ais[[ 3 ]][[i]][,,2])$BestRecov,
        FindOptimalPermutSingle( Aest3[[ i ]][,,3] , simdata$Ais[[ 3 ]][[i]][,,3])$BestRecov)
      )
    }
    
    Ais4 <- numeric()
    for(i in 1:length(simdata$Ais[[ 4 ]]) ){
      Ais4[i] <- mean(c(
        FindOptimalPermutSingle( Aest4[[ i ]][,,1] , simdata$Ais[[ 4 ]][[i]][,,1])$BestRecov,
        FindOptimalPermutSingle( Aest4[[ i ]][,,2] , simdata$Ais[[ 4 ]][[i]][,,2])$BestRecov,
        FindOptimalPermutSingle( Aest4[[ i ]][,,3] , simdata$Ais[[ 4 ]][[i]][,,3])$BestRecov)
      )
    }
    # Res A
    AresTuck <- mean(mean(Ais1), mean(Ais2), mean(Ais3), mean(Ais4))
  }else{
    SresTuck <- mean( c(
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[1]]][,,1]) , Strue = t(simdata$Qr[[1]][,,1]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[1]]][,,2]) , Strue = t(simdata$Qr[[1]][,,2]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[1]]][,,3]) , Strue = t(simdata$Qr[[1]][,,3]))$BestRecov,
      
      
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[2]]][,,1]) , Strue = t(simdata$Qr[[2]][,,1]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[2]]][,,2]) , Strue = t(simdata$Qr[[2]][,,2]))$BestRecov,
      FindOptimalPermutSingle(Sest = t(Sest[[p_per$BestPerm[2]]][,,3]) , Strue = t(simdata$Qr[[2]][,,3]))$BestRecov
      
    ))
    
    
    Aest1 <- inv_W(OptStart$IVA$Wr[[p_per$BestPerm[1]]])
    Aest1 <- Asplit(Aest1, dimsplit = grid$Ti[sim])
    
    Aest2 <- inv_W(OptStart$IVA$Wr[[p_per$BestPerm[2]]])
    Aest2 <- Asplit(Aest2, dimsplit = grid$Ti[sim])
    
    
    Ais1 <- numeric()
    for(i in 1:length(simdata$Ais[[ 1 ]]) ){
      Ais1[i] <- mean(c(
        FindOptimalPermutSingle( Aest1[[ i ]][,,1] , simdata$Ais[[ 1 ]][[i]][,,1])$BestRecov,
        FindOptimalPermutSingle( Aest1[[ i ]][,,2] , simdata$Ais[[ 1 ]][[i]][,,2])$BestRecov,
        FindOptimalPermutSingle( Aest1[[ i ]][,,3] , simdata$Ais[[ 1 ]][[i]][,,3])$BestRecov)
      )
    }
    
    Ais2 <- numeric()
    for(i in 1:length(simdata$Ais[[ 2 ]]) ){
      Ais2[i] <- mean(c(
        FindOptimalPermutSingle( Aest2[[ i ]][,,1] , simdata$Ais[[ 2 ]][[i]][,,1])$BestRecov,
        FindOptimalPermutSingle( Aest2[[ i ]][,,2] , simdata$Ais[[ 2 ]][[i]][,,2])$BestRecov,
        FindOptimalPermutSingle( Aest2[[ i ]][,,3] , simdata$Ais[[ 2 ]][[i]][,,3])$BestRecov)
      )
    }
    
    
    # Res A
    AresTuck <- mean(mean(Ais1), mean(Ais2))
  }
  
  
  output$TuckS <- SresTuck
  output$TuckA <- AresTuck
  
  return(output)  
}

for(i in 1:2){
  load(files[i])
  output <- calc(outputfile=output)
  save(output, file = files[i])
}


