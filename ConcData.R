# Thu Nov  2 11:11:03 2023 ------------------------------
# JD
# What: function to concatenate arrays based on P

ConcData <- function(DataList, ClusVec){
  ClusL <- length( unique(ClusVec) )
  
  NewList <- list()
  
  for(i in 1:ClusL){
    NewList[[i]] <- do.call(abind, list(DataList[ClusVec == i], along = 1))
  }
  return(NewList)
} 
