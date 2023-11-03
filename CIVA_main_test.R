# Thu Nov  2 13:30:38 2023 ------------------------------
# JD

# what: CIVA mainscript first test

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

set.seed(42)
Ti <- 2
testdata <- Simulate_CIVA(Nr = 5, V = 2000, Ti = Ti, Q = 2,R = 2,D = 3, E = 0.1)
X <- testdata$Xe
sum(X[[1]]^2)
X <- lapply(X, xscale) 
sum(X[[1]]^2)


Losstotal <- sum( sapply( seq_along(X),
                          function(i) sum( X[[i]]^2)) ) + 1
Losstotal
# init P
P <- clusf(nBlocks = length(X), nClus = 2)
P
Pw <- c(1,2)
Pw <- c(1,1,1,1, 2,2,2,2)

DataL <- ConcData(X,P)

DataL <- ConcData(X,UpdatePInfo$newclus)

IVAL <- ExtractIVA(DataList = DataL, nComp = 2, par = TRUE, vars = ls(), seed = sample(1:100000,size = 1))


UpdatePInfo <- Reclus(DataList = X, QrL = IVAL$Qr)
UpdatePInfo$newclus
Losstotal <- c(Losstotal,UpdatePInfo$Loss)
Losstotal

Losstotalold <- tail(Losstotal,n = 1)

S1 <- t(IVAL$Qr[[1]][,,1])
Strue <- t(testdata$Qr[[2]][,,1])
congru(S1,Strue)

S2 <- t(IVAL$Qr[[2]][,,2])
Strue2 <- t(testdata$Qr[[2]][,,2])
congru(S2,Strue2)

Aest1 <- inv_W(IVAL$Wr[[2]])
Aest1 <- Asplit(Aest = Aest1, dimsplit = Ti)
Aest1 <- Aest1[[1]]

Aest2 <- inv_W(IVAL$Wr[[2]])
Aest2 <- Asplit(Aest = Aest2, dimsplit = Ti)

Atrue <- testdata$Ais
Atrue1 <- Atrue[[1]][[1]]
Atrue2 <- Atrue[[2]]

#A1 <- t(IVAL$Qr[[1]][,,1])

congru(Aest1[,,1],Atrue1[,,1])
congru(Aest1[,,2],Atrue1[,,2])
congru(Aest1[,,3],Atrue1[,,3])

congru(Aest1[[2]][,,1],Atrue1[[2]][,,1])
congru(Aest1[[1]][,,1],Atrue1[[1]][,,1])
congru(Aest1[[2]][,,1],Atrue1[[2]][,,1])

congru(Aest2[[1]][,,1],Atrue2[[1]][,,1])
congru(Aest2[[2]][,,1],Atrue2[[2]][,,1])
congru(Aest2[[3]][,,1],Atrue2[[3]][,,1])
congru(Aest2[[4]][,,1],Atrue2[[4]][,,1])
