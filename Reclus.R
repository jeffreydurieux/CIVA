# Thu Nov  2 13:24:11 2023 ------------------------------
# JD:
# What: recluster based on Lir

Reclus <- function(DataList, QrL) {
  
  SSList <- lapply(DataList, FUN = Lir, QrL = QrL)
  SSminIndexVec <- sapply(SSList, FUN = which.min)
  SSminVec <- sapply(SSList, FUN = min)
  Loss <- sum(SSminVec)
  
  ResList <- list("newclus" = SSminIndexVec,'SSminVec' = SSminVec,
                  "Loss" = Loss, "SSList" = SSList)
  
  return(ResList)
}
