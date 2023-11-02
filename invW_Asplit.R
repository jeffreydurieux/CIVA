# Thu Nov  2 11:07:27 2023 ------------------------------
# JD:
# What: helper functions CIVA: inverse of unmixing matrices and split A into list

inv_W <- function(W){
  # W = unmixing matrices from CIVA (res$W)
  dims <- dim(W)
  arrA <- array(data = NA, dim = c(dims[2],dims[1],dims[3]))
  for(i in 1:dim(W)[3]){
    arrA[,,i] <- mpinv(W[,,i])
  }
  return(arrA)
}

Asplit <- function(Aest, dimsplit){ 
  # Aest = result from inv_W
  #dimsplit = Ti
  len <- dim(Aest)[1]
  ids <- seq(1,(len-dimsplit+1), by = dimsplit)
  
  Al <- list()
  for(i in 1:length(ids)){
    Al[[i]] <- Aest[ids[i]:(ids[i]+(dimsplit-1)),,]
  }
  return(Al)
}