################## Function: getdistances ####################
#
# This function calculates distances 
# Input: - inputdata: Expressionmatrix with probesets as rownames
#	     - calls: matrix of calls, sorted exactly like inputdata, probesets as rownames
#        - GOIDs: GOIDs to use as returned by defineGOIDs
#        - classes: sample affiliations
#        - outputfile: store output file
#        - chip: Array type to use, e.g. "hgu133a", "illuminaHumanv3"
#        - permutations: number of permutations to calculate significance
#
#############################################################

getdistances_calls <- function(inputdata, calls=T, calls.file, numpresent=5, GOIDs, classes, outputfile, chip="hgu133a",permutations=1000) 
{
  # load required libraries and functions
  eval(parse(text=paste("require(", chip , ".db)",sep="")))
  require(affy)
    
  newifelse <- function(row1,row2) 
  {
    return(ifelse(row1>row2,1,0))
  }

  checkpresent <- function(vec,nrgroups)
  {
    absent <- length(which(vec=="A"))
    if(absent > length(vec)/nrgroups)
    {
      return(FALSE)
    }
    else
    {
      return(TRUE)
    }
  }  

  getGOterms <- function(GOIDs)
  {
    # load required libraries and functions
    require(GO.db)
    require(annotate)

    description <- mget(GOIDs,env=GOTERM)
    terms <- vector()

    for (i in 1:length(description))
    {
      terms <- c(terms,Term(description[[i]] ))
    }
    return(terms)
  }

  ###########################################################################
  
  # start with function

  # get group and sample information #
  groups <- as.integer(as.factor(classes))
  classes_r <-  classes[order(groups)]
  nrgroups <- length(levels(as.factor(classes)))
  nrsamples <- length(as.factor(classes))

  # define new list of GOIDs
  #checked <- sapply(GOIDs,checkGOID,calls,chip)
  #GOIDs <- GOIDs[checked]

  errors <- vector()

  for(k in 1:length(GOIDs))
  {
    if (calls)
    {
	  eval(parse(text=paste("APcalls <- calls.file[unique(get(GOIDs[k],env=",chip,"GO2ALLPROBES)),]",sep="")))
      APcalls.checked <- apply(APcalls,1,checkpresent,nrgroups)
      if( length(which(APcalls.checked)) > numpresent )
      {
        ids <- names(which(APcalls.checked))
        data <- inputdata[ids,]
      }
    }  
    else
    {
      eval(parse(text=paste("data <- inputdata[unique(get(GOIDs[k],env=",chip,"GO2ALLPROBES)),]",sep="")))	  
	    #data <- inputdata
    }
            
    print(date())
    print(paste("Work on GOID #",k,":",GOIDs[k],"..."))
    

    mean_GO_all<-matrix(nrow=nrow(data),ncol=nrgroups)

    for (l in 1:nrgroups)
    {
      mean_GO_all[,l] <- apply(as.matrix(data[,groups==l]),1,mean)
    }

    mean_GO.dist<-dist(t(mean_GO_all))

    perm.mean.dist_GO = matrix(nrow=permutations, ncol=choose(nrgroups,2))


    for (i in 1:permutations)
    {
      index.rand = sample(1:nrsamples,nrsamples)
      data = data[,index.rand]

      means = matrix(nrow=nrow(data), ncol=nrgroups)
      for (j in 1:nrgroups)
      {
        means[,j] = apply(as.matrix(data[,groups==j]),1, mean)
      }

      perm.mean.dist_GO[i,] = as.vector(dist(t(means)))
    }  

    # compare distances of "real" subgroups to distances of permutations #
    dist.comb.GO<-t(apply(perm.mean.dist_GO,1,newifelse,mean_GO.dist))
    # configuring the output #
    error <- GOIDs[k]
    error <- c(error,getGOterms(GOIDs[k]))
    error1 <- GOIDs[k]
    error1 <- cbind(error1,getGOterms(GOIDs[k]))
    if (nrgroups > 2)
    {
      for (m in 1:choose(nrgroups,2))
      {
        error1 <- cbind(error1,(sum(dist.comb.GO[,m]))/permutations)
        error <- c(error,(sum(dist.comb.GO[,m]))/permutations)
      }
      errors <- rbind(errors,error)
    }
    else
    {
      error1 <- cbind(error1,(sum(dist.comb.GO))/permutations)
    }

### GET COLUMN NAMES FOR P-VALUES
    p1 <- NULL
    p2 <- NULL
    for(i in 1:nrgroups){
       if(i < nrgroups){
          p1 <- c(p1, unique(classes_r)[rep(i,(nrgroups - i))]) 
          p2 <- c(p2, unique(classes_r)[(i+1):nrgroups])           
       }
    }
     if(k == 1) {
        colnames(error1) <- c("GO ID","GO Term",paste("p-value ", paste(p1,p2,sep=" vs "), sep = ""))
        write.table(error1,file=outputfile,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
     } else {
        write.table(error1,file=outputfile,append = TRUE, sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)  
     }
  }
  return(errors)
}

