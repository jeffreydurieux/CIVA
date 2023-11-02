# Mon Oct 30 16:50:02 2023 ------------------------------
# JD

# What: function for Aihat reconstruction


EstAih <- function(Xi, Sr){
  dims <- dim(Xi)
  dimss <- dim(Sr)
  Aih <- array(data = NA, dim = c(dims[1],dimss[1],dims[3]))
  for(i in 1:dims[3]){
    Aih[ , , i] <- Xi[,,i] %*% t(Sr[,,i]) %*% NMFN::mpinv(Sr[,,i] %*% t(Sr[,,i]))
  }
  
  return(Aih)
}

# test <- EstAih(X, Sr = res$S)            
# # fastIVA W are unmixing matrices, take inverse 
# apply(res$W, 3, solve, simplify = F)
