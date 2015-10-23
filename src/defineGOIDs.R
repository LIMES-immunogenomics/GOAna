################## Function: defineGOIDs ####################
#
# This function defines the list of GOIDs the analysis will be run on
# Input: - chip: Array type
#        - numprobesets_low: the minimal number of probesets per GOID
#        - numprobesets_high: the maximum number of probesets per GOID
#        - category: GO category of IDs
#        - GO.node: point in the GO hierarchy to start
#
#############################################################

defineGOIDs <- function(chip="hgu133a",numprobesets_low=5,numprobesets_high=10000,is.list=F,list=NULL,category="BP",GO.node="GO:0008150")
{  

  # load required libraries and functions
  eval(parse(text=paste("require(" , chip , ".db)",sep="")))
  require(GO.db)
  require(annotate)

  # get GO-Offspring for specified category and root node
  eval(parse(text=paste("GOIDs.category.node = GO",category,"OFFSPRING$","'",GO.node,"'", sep="")))

  # get GOIDs for specified array
  eval(parse(text=paste("GOIDs.array = ls(", paste(chip,"GO2ALLPROBES",sep=''), ")", sep="")))

  # get list of probesets for GOIDs in intersection of both
  eval(parse(text=paste("listGO <- mget(intersect(GOIDs.array,GOIDs.category.node),env=", paste(chip,"GO2ALLPROBES",sep=''),")",sep="")))

  GOIDs <- vector()
  for (i in 1:length(listGO))
  {
    if (length(listGO[[i]]) > numprobesets_low && length(listGO[[i]]) < numprobesets_high)
    {
      if (is.list)
      {
        if(any(is.element(listGO[[i]],list)))
        {
          GOIDs <- c(GOIDs,names(listGO[i]))
        }
      }
      else
      {
        GOIDs <- c(GOIDs,names(listGO[i]))
      }
    }
  }

  return(GOIDs)
}


