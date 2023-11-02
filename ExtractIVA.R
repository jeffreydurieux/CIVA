# Thu Nov  2 11:34:47 2023 ------------------------------
# JD

# what: perform IVA per cluster

ExtractIVA <- function(DataList, nComp, par = F, vars=NULL,seed=1){
  
    if(par == TRUE){
      R <- length(DataList)
      if(R >= 15){
        R=15
      }
      cl <- makeCluster(R)
      clusterExport(cl = cl, varlist = vars)
      ivaListCluster <- parLapply(cl = cl, 1:length(DataList), fun = function(x){
        set.seed(seed+x)
        fastIVA_custom(X = DataList[[x]], nc = nComp)
      })
      stopCluster(cl)
    }else{
      ivaListCluster <- lapply(DataList, fastIVA_custom, nc = nComp)
    }
    
    IVA_S <- lapply(ivaListCluster, function(anom) anom$S)
    IVA_W <- lapply(ivaListCluster, function(anom) anom$W)
    ListRes <- list("Qr" = IVA_S, "Wr" = IVA_W)

  return(ListRes)
}
