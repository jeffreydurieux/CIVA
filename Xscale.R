# Thu Nov  2 14:15:10 2023 ------------------------------
# JD

# What: scaling function for blockscaling

xscale <- function(Xi, value = 1000){
  dims <- dim(Xi)
  Xin <- array(data = NA, dim = dims)
  
  for(i in 1:dims[3]){
    f <- sqrt(value/sum(Xi[,,i]^2))
    Xin[,,i] <- f*Xi[,,i]
  }
  
  return(Xin)
}