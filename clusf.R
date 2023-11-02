clusf <- function(nBlocks, nClus) {
  #simplyfied cluster generation function using an equal probability
  clus <- GenerateRandomClustering(nBlocks, nClus, rep(c(1 / nClus), nClus))
  return(clus)
}


AdjustProb <- function(v , MaxElem)
{
  # INPUT
  #   v (1 x nElem): vector
  #   MaxElem: number of elements
  #
  # OUTPUT
  #   v (1 x nElem): vector with 'sum(out)=MaxElem' and no 0-cells
  
  nElem = length(v)
  
  if (any(v < 1))
    #add 1 to 0-cells
  {
    tempv = rep(0,nElem)
    tempv[v < 1] = rep(1 , length(which(v < 1)))
    v = v + tempv
  }#### replace(v,which(v<1),1)
  
  it <- 1
  while (!(sum(v) == MaxElem) | it == 500000)
  {
    if(it == 499999){
      stop('No suitable start clustering found, select fewer clusters or select clusters as max number of objects (i.e., length(DataList)')
    }
    
    it <- it + 1
    diff = sum(v) - MaxElem
    if (diff < 0)
      #add elements
    {
      if (abs(diff) < (nElem - sum(v == 1)))
      {
        for (tel in 1:abs(diff))
        {
          tempcl = ceiling(runif(1) * nElem)
          while (v[tempcl] == 1)
          {
            tempcl = ceiling(runif(1) * nElem)
          }
          v[tempcl] = v[tempcl] + 1
        }
      }
      else
      {
        v[which(!(v == 1))] = v[which(!(v == 1))] + rep(1 , length(which(!(v == 1))))
      }
    }
    else
      #delete elements
    {
      if (abs(diff) < (nElem - sum(v == 1)))
      {
        for (tel in 1:abs(diff))
        {
          tempcl = ceiling(runif(1) * nElem)
          while (v[tempcl] == 1)
          {
            tempcl = ceiling(runif(1) * nElem)
          }
          v[tempcl] = v[tempcl] - 1
        }
      }
      else
      {
        v[which(!(v == 1))] = v[which(!(v == 1))] - rep(1 , length(which(!(v == 1))))
      }
    }
  }
  
  return(v)
}


GenerateRandomClustering <- function(nElement , nClust , Prob = NULL)
{
  ####GenerateRandomClustering = for Random Starts
  
  # Author: Tom F. Wilderjans
  # nElement: number of elements to be clustered
  # nClust: number of clusters
  # Prob (1 x nClust): proportion of elements in each cluster
  
  # Added by Jeffrey Durieux: default Prob = equal cluster prob
  # This done to adjust code later on for potential cluster perbutation
  
  if(is.null(Prob))
  {
    Prob <- rep(x = (1/nClust) , nClust)
  }
  
  
  BestClust = NULL
  ErrorEncountered = F
  
  if (!(length(Prob) == nClust))
  {
    ErrorEncountered = T
  }
  
  if ((abs(sum(Prob) - 1) > .000000001) | (any(Prob < 0)))
  {
    
    ErrorEncountered = T
  }
  
  if (!(any(nClust == 1:nElement)))
  {
    ErrorEncountered = T
  }
  
  if (!(ErrorEncountered))
  {
    if (nElement > nClust)
    {
      if (nClust == 1)
      {
        BestClust = rep(1 , times = nElement)
      }
      else
      {
        ProbVV = round(Prob * nElement)
        if (!(sum(ProbVV) == nElement) |
            (any(ProbVV < 1)))
          #not enough elements, or empty clusters
        {
          ProbVV = AdjustProb(ProbVV , nElement)
        }
        
        tempclus = rep(1:length(ProbVV) , ProbVV)
        BestClust = tempclus[sample(1:nElement,size = nElement,replace =
                                      FALSE)]
      }
    }
    else
    {
      BestClust = 1:nClust
    }
  }
  
  if (!(length(unique(BestClust)) == nClust))
  {
    BestClust = NULL
  }
  
  return(BestClust)
}