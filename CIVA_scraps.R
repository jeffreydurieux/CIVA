library(ivaBSS)
library(LaplacesDemon)

### check: https://www.mdpi.com/1424-8220/23/11/5333


# Generate sources from multivariate Laplace distribution
P <- 5; N <- 50000; D <- 3; # D is number of data sets # P dimension, # N sample
S <- array(NA, c(P, N, D))
  
  for (i in 1:P) {
    S[i, , ] <- rmvl(N, rep(0, D), diag(D))
  }
  
  # Generate mixing matrices from standard normal distribution
  A <- array(runif(100 * P * D), c(100, P, D))
  dim(A)
  # Generate mixtures
  X <- array(NaN, c(100, N, D))
  dim(X)
  for (d in 1:D) {
    X[, , d] <- A[, , d] %*% S[, , d]
  }
  
  # Estimate sources and unmixing matrices
  set.seed(2407)
  tmp <- proc.time()
  res <- fastIVA(X, source_density = 'laplace_diag')
  time <- proc.time()-tmp
  time
  
  library(multiway)

congru(t(S[,,1]), t(res$S[,,1])) 
congru(t(S[,,2]), t(res$S[,,2])) 
congru(t(S[,,3]), t(res$S[,,3])) 


A[,,1]
t=coef(res, which.dataset = 1)
congru(solve(t),A[,,1])
cor(A[,,1],res$W[,,1])
A[,,1]
res$W[,,1]

Xhat <- array(NaN, c(P, N, D))
for(d in 1:D){
  Ah <-   crossprod(t(X[,,d]), t(res$S[,,d])) %*% mpinv(res$S[,,d] %*% t(res$S[,,d]))
  Xhat[,,d] <- Ah %*% res$S[,,d]
}

sum((X[,,1]-Xhat[,,1])^2)+
sum((X[,,2]-Xhat[,,2])^2)+
sum((X[,,3]-Xhat[,,3])^2)

sum((X-Xhat)^2)


#### pca preproc test###