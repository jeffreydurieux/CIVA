# Fri Nov  3 15:11:18 2023 ------------------------------
# JD
# What: first basic version of CIVA

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

CIVA <- function(X, nc, R, starts, parallel = FALSE, verbose = TRUE){
 
  X <- lapply(X, xscale) 
  
  Losstotal <- sum( sapply( seq_along(X),
                            function(i) sum( X[[i]]^2)) ) + 1
  allstarts <- list()
  for(i in 1:starts){
    
    iter <- 0
    repeat({
      iter <- iter + 1
      
      #1 init P
      if(iter == 1){
        Losstotalstart <- Losstotal
        P <- clusf(nBlocks = length(X), nClus = R)
        DataL <- ConcData(X,P)
      }else{
        DataL <- ConcData(X,UpdatePInfo$newclus)
      }
      
      #2 extract IVA
      IVAL <- ExtractIVA(DataList = DataL, nComp = nc, par = parallel, vars = ls(), seed = sample(1:100000,size = 1))
      
      #3 Reclus
      UpdatePInfo <- Reclus(DataList = X, QrL = IVAL$Qr)
      Losstotalstart<- c(Losstotalstart,UpdatePInfo$Loss)
      
      if(verbose == TRUE){
        cat('Loss ', tail(Losstotalstart, n=1), '\n')
      }
      
      #4 Convergence
      criterion <- tail(Losstotalstart, n = 2)
      if(criterion[2] == criterion[1]){
        cat('Convergence of start ', i, '\n')
        start <- list()
        start$P <- UpdatePInfo$newclus
        start$IVA <- IVAL
        start$LossStart <- UpdatePInfo$Loss
        break()
      }
      
      

    })# end repeat
  allstarts[[i]] <- start
  }# end for
  
  #out <- list()
  # extract optimal start
  
   return(allstarts)
}

set.seed(2407*3)
Ti <- 10
starts <- 10
Q <- 5
testdata <- Simulate_CIVA(Nr = 20, V = 1000, Ti =Ti, Q = Q, R = 2,D = 3, E = 0.01)

tmp <- proc.time()
testciva <- CIVA(X = testdata$Xe, nc = Q, R = 2, starts = starts)
time <- proc.time() - tmp
time

Pests <- sapply(seq_along(testciva), function(x) testciva[[x]]$P)
Pests
ARI_starts <- apply(Pests, MARGIN = 2, mclust::adjustedRandIndex, y = testdata$P)
ARI_starts
cat(length( which(ARI_starts==1) ) /starts * 100,'%')

losses <- sapply(seq_along(testciva), function(x) testciva[[x]]$LossStart)

Optimal <- which.min(losses)

OptStart <- testciva[[Optimal]]

ARI <- mclust::adjustedRandIndex(testdata$P, OptStart$P)

#source('https://raw.githubusercontent.com/jeffreydurieux/usefulcode/master/FindOptimalClusPermut.R')
#source('https://raw.githubusercontent.com/jeffreydurieux/usefulcode/master/FindOptimalPermutSources.R')
#library(gtools)
p_per <- FindOptimalClusPermut(OptStart$P, testdata$P)

# recov Q1
SresTuck <- mean(
  mean(c(
FindOptimalPermutSingle( t(testdata$Qr[[ p_per$BestPerm[1] ]][,,1]) , t(OptStart$IVA$Qr[[1]][,,1]) )$BestRecov,
FindOptimalPermutSingle( t(testdata$Qr[[ p_per$BestPerm[1] ]][,,2]) , t(OptStart$IVA$Qr[[1]][,,2]) )$BestRecov,
FindOptimalPermutSingle( t(testdata$Qr[[ p_per$BestPerm[1] ]][,,3]) , t(OptStart$IVA$Qr[[1]][,,3]) )$BestRecov)
)
,
# recov Q2
mean(c(
  FindOptimalPermutSingle( t(testdata$Qr[[ p_per$BestPerm[2] ]][,,1]) , t(OptStart$IVA$Qr[[2]][,,1]) )$BestRecov,
  FindOptimalPermutSingle( t(testdata$Qr[[ p_per$BestPerm[2] ]][,,2]) , t(OptStart$IVA$Qr[[2]][,,2]) )$BestRecov,
  FindOptimalPermutSingle( t(testdata$Qr[[ p_per$BestPerm[2] ]][,,3]) , t(OptStart$IVA$Qr[[2]][,,3]) )$BestRecov)
)
)
Aest1 <- inv_W(OptStart$IVA$Wr[[1]])
Aest1 <- Asplit(Aest1, dimsplit = Ti)

Aest2 <- inv_W(OptStart$IVA$Wr[[2]])
Aest2 <- Asplit(Aest2, dimsplit = Ti)

Ais1 <- numeric()
for(i in 1:length(testdata$Ais[[ p_per$BestPerm[1]]]) ){
Ais1[i] <- mean(c(
  FindOptimalPermutSingle( testdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,1] , Aest1[[i]][,,1] )$BestRecov,
  FindOptimalPermutSingle( testdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,2] , Aest1[[i]][,,2] )$BestRecov,
  FindOptimalPermutSingle( testdata$Ais[[ p_per$BestPerm[1] ]][[i]][,,3] , Aest1[[i]][,,3] )$BestRecov)
)  
}

Ais2 <- numeric()
for(i in 1:length(testdata$Ais[[ p_per$BestPerm[2]]]) ){
  Ais2[i] <- mean(c(
    FindOptimalPermutSingle( testdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,1] , Aest2[[i]][,,1] )$BestRecov,
    FindOptimalPermutSingle( testdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,2] , Aest2[[i]][,,2] )$BestRecov,
    FindOptimalPermutSingle( testdata$Ais[[ p_per$BestPerm[2] ]][[i]][,,3] , Aest2[[i]][,,3] )$BestRecov)
  )  
}
# Res A
AresTuck <- mean(mean(Ais1), mean(Ais2))


ARI; SresTuck; AresTuck
