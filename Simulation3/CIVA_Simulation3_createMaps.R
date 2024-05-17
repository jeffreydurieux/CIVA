# Thu Nov  9 10:17:25 2023 ------------------------------
# JD
# What: CIVA simulation 3. Realistic maps and dependency over sessions

# use two clusters

# both clusters have med visual and DMN intact in ses1

# for intact: assume same map across clusters
# one cluster have intact med visual and DMN intact in ses2
# one cluster have intact med visual and slightly affected DMN in ses2

# one cluster have intact med visual and DMN intact in ses3
# one cluster have intact med visual and severe affected DMN in ses2

# hypo: clustering will be driven completely by disrupted DMN


library(CICA)
n <- readNifti('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Data/23x_all_beckmann_8_and_csf_and_wm.nii.gz')
mask <- readNifti('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Data/brainmask_23_28_23.nii.gz')
mask <- matrix(mask)
maskid <- which(mask == 1)

dim(n)
Nif <- matrix(n, ncol = dim(n)[4])
#view(n[,,,6])
#view(n[,,,5])
#view(mask)



for(i in 1:10){
  set.seed(i)
  S1 <- Nif[maskid,c(1,5)] #### med vis stable / DMN deteriorate
  S2 <- Nif[maskid,c(1,5)] #### med vis stable / DMN stable
  
  
  #cluster 1 corr signals
  # session 2
  #set.seed(2407)
  
  #stable
  sig1 <- which(S1[,1]!=0)
  zeroid <- sample(sig1, size = ceiling(length(sig1)*.05) ) #~ 91
  S1_1ses2 <- S1[,1]
  S1_1ses2[zeroid] <- 0 
  cor(S1[,1], S1_1ses2)
  
  # deteriorate here
  sig1 <- which(S1[,2]!=0)
  zeroid <- sample(sig1, size = ceiling(length(sig1)*.4) ) #~ 91
  S1_2ses2 <- S1[,2]
  S1_2ses2[zeroid] <- 0 
  cor(S1[,2], S1_2ses2)
  
  
  
  # session 3
  #set.seed(2407)
  
  #stable
  sig1 <- which(S1_1ses2!=0)
  zeroid <- sample(sig1, size = ceiling(length(sig1)*.05) ) #~ 91
  S1_1ses3 <- S1_1ses2
  S1_1ses3[zeroid] <- 0 
  cor(S1[,1], S1_1ses2) #.95
  cor(S1_1ses2, S1_1ses3) #.97
  cor(S1[,1], S1_1ses3) #.93
  
  #set.seed(2407)
  #deteriorate
  sig2 <- which(S1_2ses2!=0)
  zeroid <- sample(sig2, size = ceiling(length(sig2)*.4) ) #~ 91
  S1_2ses3 <- S1_2ses2
  S1_2ses3[zeroid] <- 0 
  cor(S1[,2], S1_2ses2) #.74
  cor(S1_2ses2, S1_2ses3) #.86
  cor(S1[,2], S1_2ses3) #.64
  
  
  
  #cluster 2 corr signals ### make a relatively stable cluster
  # session 2
  #set.seed(2407)

  # #stable
  sig1 <- which(S2[,1]!=0)
  zeroid <- sample(sig1, size = ceiling(length(sig1)*.05) ) #~ 91
  S2_1ses2 <- S2[,1]
   S2_1ses2[zeroid] <- 0
   cor(S2[,1], S2_1ses2) #.96
  #
  # #stable
   sig1 <- which(S2[,2]!=0)
   zeroid <- sample(sig1, size = ceiling(length(sig1)*.05) ) #~ 91
   S2_2ses2 <- S2[,2]
   S2_2ses2[zeroid] <- 0
   cor(S2[,2], S2_2ses2) #.97
  #
  #
  #
  # # session 3
  # #set.seed(2407)
  #
  # #stable
   sig1 <- which(S2_1ses2!=0)
   zeroid <- sample(sig1, size = ceiling(length(sig1)*.05) ) #~ 91
   S2_1ses3 <- S2_1ses2
   S2_1ses3[zeroid] <- 0
   cor(S2[,1], S2_1ses2) #.90
   cor(S2_1ses2, S2_1ses3) #.97
   cor(S2[,1], S2_1ses3) #.93

  # #set.seed(2407)
  # #stable
   sig2 <- which(S2_2ses2!=0)
   zeroid <- sample(sig2, size = ceiling(length(sig2)*.05) ) #~ 91
   S2_2ses3 <- S2_2ses2
   S2_2ses3[zeroid] <- 0
   cor(S2[,2], S2_2ses2) #.97
   cor(S2_2ses2, S2_2ses3) #.97
   cor(S2[,2], S2_2ses3) #.94
  
  
  #### combine signals in array
  library(abind)
  #S1 
  QS11 <- S1
  QS12 <- cbind(S1_1ses2, S1_2ses2)
  QS13 <- cbind(S1_1ses3, S1_2ses3)
  
  #S2
  QS21 <- S1
  QS22 <- cbind(S2_1ses2, S2_2ses2)
  QS23 <- cbind(S2_1ses3, S2_2ses3)
  
  S1 <- abind(t(QS11),t(QS12),t(QS13), along = 3)
  S2 <- abind(t(QS21),t(QS22),t(QS23), along = 3)
  
   cor(t(S1[,,1]), t(S1[,,2]))
   cor(t(S1[,,2]), t(S1[,,3]))
   cor(t(S1[,,1]), t(S1[,,3]))
  # 
   
   cor(t(S2[,,1]), t(S2[,,2]))
   cor(t(S2[,,2]), t(S2[,,3]))
   cor(t(S2[,,1]), t(S2[,,3]))
  # 
   cor(t(S1[,,1]), t(S2[,,1])) #diag should be ones
   
   ext <- paste('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/Simulation3/Maps/',
               'Maps_S1_',i,'.Rdata', sep = '')
  save(S1, file = ext)
  
  ext <- paste('C:/Users/jeffr/OneDrive - Erasmus University Rotterdam/Documents/CIVA/CIVA/Simulation3/Maps/',
               'Maps_S2_',i,'.Rdata', sep = '')
  save(S2, file = ext)
}
