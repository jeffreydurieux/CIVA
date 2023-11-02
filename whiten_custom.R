whiten_custom <- function (X, nc) 
{
  P <- dim(X)[1]
  N <- dim(X)[2]
  D <- dim(X)[3]
  X_centered <- array(0, c(P, N, D))
  Z <- array(0, c(nc, N, D))
  V <- array(0, c(P, nc, D))
  X_means <- array(0, c(P, D))
  for (d in 1:D) {
    Xd <- t(X[,,d])
    nobs <- nrow(Xd)-1

    Xd <- scale(Xd, scale = FALSE)
    xeig <- eigen(crossprod(Xd)/nobs, symmetric = TRUE)
    
    Dmat <- diag(sqrt(xeig$values[1:nc]))
    Mprt <- tcrossprod(Dmat, xeig$vectors[, 1:nc, drop = FALSE])
    diag(Dmat) <- 1/diag(Dmat)
    
    if(nc < P){
      Pmat <- xeig$vectors[, 1:nc, drop = FALSE] %*% Dmat 
    }else{
      Pmat <- xeig$vectors[, 1:nc, drop = FALSE] %*% Dmat %*% t(xeig$vectors[, 1:nc, drop = FALSE])
      
    }
    
    
    Xw <- Xd %*% Pmat
    Z[,,d] <- t(Xw)
    V[,,d] <- Pmat
  }
  return(list(Z = Z, V = V))
}


# test <- whiten_custom(X,2)
# cov(t(test$Z[,,3]))
# plot(t(test$Z[,,3]), asp=T)
# 
# test2 <- ivaBSS:::whiten(X)
# cov(t(test2$Z[,,1]))
# plot(t(test2$Z[,,3]), asp=T)
