library(multiway)

# sca-p

Xl <- list()
for(i in 1:3){
  Xl[[i]] <- t(X[,,i])
}

scares <- sca(Xl, nfac = 2, type = 'sca-p', rotation = 'none')

arr <- array(data = 0, dim = c(4,1000,3))

for(i in 1:3){
  arr[,,i] <- t( scares$D[[i]] )
}

set.seed(2407)
tmp <- proc.time()
res <- fastIVA(arr, source_density = 'laplace_diag')
time <- proc.time()-tmp
time

library(multiway)

congru(t(test$Qr[[1]][,,1]), t(res$S[,,1])) 
congru(t(test$Qr[[1]][,,2]), t(res$S[,,2])) 
congru(t(test$Qr[[1]][,,3]), t(res$S[,,3])) 


#### svd version

XX <- do.call(rbind, Xl)

svdr <- svd(scale(XX))
svdr$d
XR <- svdr$u[,1:2]%*%diag(svdr$d[1:2])

cor(XR[1:1000,],scares$D[[1]])
