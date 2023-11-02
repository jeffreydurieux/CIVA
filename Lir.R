# Thu Nov  2 12:41:15 2023 ------------------------------
# JD

#What partition criterion CIVA


Lir <- function(Xi, QrL)
{
  # compute Ahats
  Ah <- lapply(seq_along(QrL), function(i){
    EstAih(Xi = Xi, Sr = QrL[[i]])
  })
    
  
  ### for i in ah
  # Xih <- Ahs %*% Qr
  XhatL <- list()
  for(i in 1:length(Ah)){
    Xhat <- array(NA, dim = dim(Xi))  
    for(d in 1:dim(Xi)[3]){
      Xhat[,,d] <- Ah[[i]][,,d] %*% QrL[[i]][,,d]
    } # end d
    XhatL[[i]] <- Xhat
  }
  
  # loss per data matrix
  SS <- sapply(seq_along(XhatL), function(i) {
    sum((Xi - XhatL[[i]]) ^ 2)
  })
  
  return(SS)
}#function