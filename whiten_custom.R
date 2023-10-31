whiten_custom <- function (X, nc) 
{
  P <- dim(X)[1]
  N <- dim(X)[2]
  D <- dim(X)[3]
  X_centered <- array(0, c(P, N, D))
  Z <- array(0, c(nc, N, D))
  V <- array(0, c(P, P, D))
  X_means <- array(0, c(P, D))
  for (d in 1:D) {
    Xd <- t(X[,,d])
    nobs <- nrow(Xd)

    Xd <- scale(Xd, scale = FALSE)
    xeig <- eigen(crossprod(Xd)/nobs, symmetric = TRUE)
    
    Dmat <- diag(sqrt(xeig$values[1:nc]))
    Mprt <- tcrossprod(Dmat, xeig$vectors[, 1:nc, drop = FALSE])
    diag(Dmat) <- 1/diag(Dmat)
    Pmat <- xeig$vectors[, 1:nc, drop = FALSE] %*% Dmat
    Xw <- Xd %*% Pmat
    Z[,,d] <- t(Xw)
  }
  return(list(Z = Z))
}

# test <- whiten_custom(X,4)
# cov(t(test$Z[,,2]))
# plot(t(test$Z[,,3]), asp=T)
# 
# test2 <- ivaBSS:::whiten(X)
# cov(t(test2$Z[,,1]))
# plot(t(test2$Z[,,3]), asp=T)
