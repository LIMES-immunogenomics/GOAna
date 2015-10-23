################## Function: CreateCytoscape ####################
#
# This function extends the significant GOIDs and creates a contrib-network 
# Input: - GOIDs: a vector of GOIDs to analyze
#	       - calls: T or F (TRUE or FALSE), wether you have a call (A,P,M) for each ProbeID on the chip
#        - datacalls: if calls == T, then you have to give the file which includes the calls
#        - intersection: determines the number of shared GOIDs between two genes
#        - chip: Array type to use, e.g. "hgu133a", "illuminaHumanv3"
#        - comparison: in the case of more than two groups, specify the pairwise comparison of interest 	
#						as character, it will be used for the output file names	(e.g. "A_vs_B")
# 
#
#############################################################



createCytoscape<-function(GOIDs,calls=T,datacalls,intersection=2,chip="hgu133a", comparison = "1")
{
  # load required libraries and functions
  eval(parse(text=paste("require(", chip , ".db)",sep="")))
  require(GO.db)
  require(annotate)
  
  
  getSymbol <- function(vec,chip)
  {
    eval(parse(text=paste("require(", chip , ".db)",sep="")))
    eval(parse(text=paste("SYMBOLs <- apply(as.matrix(vec),1,retrieveSYMBOL<-function(ID){SYMBOL <- get(ID,env=",chip,"SYMBOL)})",sep="")))
    return(SYMBOLs)
  }

  checkpresent <- function(vec)
  {
    absent <- length(which(vec=="A"))
    if(absent > length(vec)/3)
    {
      return(FALSE)
    } 
    else
    {
      return(TRUE)
    }
  }

  is.greater <- function(vec,number)
  {
    if (length(vec)>number)
    {
      return(TRUE)
    }
    else
    {
      return(FALSE)
    }
  }
  

#  if (!exists("newlist"))
#  {
    eval(parse(text=paste("WorkList <- mget(GOIDs,env=",chip,"GO2ALLPROBES)",sep="")))
    print("Creating worklist ...") 
#    print(WorkList)
    if(calls)
    {
      for ( i in 1:length(WorkList) )
      {
        APcalls <- datacalls[WorkList[[i]],]
        APcalls.checked <- apply(APcalls,1,checkpresent)
        WorkList[[i]] <- WorkList[[i]][APcalls.checked]
      }
    }
  
  mat <- vector()
  print("Creating opposite view ...") 
  
  # work on List of GODIDs to create the opposite view 
  # -> a list of genes with resp. GOIDs (newlist)
  for ( i in 1:length(WorkList))  
  {  
    mat <- rbind(mat,cbind(unique(as.vector(unlist(WorkList[i]))),rep(as.character(GOIDs[i]),length(unique(as.vector(unlist(WorkList[i])))))))
  }

  base <- as.integer(as.factor(mat[,1]))
  newlist <- list()
  eval(parse(text=paste("newlist <- list(\'" ,unique(mat[base==1,1]), "\'=mat[base==1,2])",sep="")))
  for ( i in 2:max(base))
  {
    eval(parse(text=paste("newlist <- c(newlist,list(\'" ,unique(mat[base==i,1]), "\'=mat[base==i,2]))",sep="")))
  }
  assign("newlist",newlist,env=.GlobalEnv)
  
#  }
#  else
#  {
#    newlist <- get("newlist",env=.GlobalEnv)
#  }

  # Create output files for Cytoscape:
  # Genes will be nodes; two genes are connected by an edge if the list of GOIDs
  # they're assigned to overlap to a specified extend (intersection)
  
  genes <- vector()
  counter <- 0
  print("Checking for intersections ...")

  newlist <- newlist[which(sapply(newlist,is.greater,1))]

  for ( i in 1:(length(newlist)-1) )
  {
    for ( j in (i+1):length(newlist) )
    {
      if( length(intersect(newlist[[i]],newlist[[j]]))>= intersection )
      {
        #print(paste("i",i,"j",j))
        counter <- counter+1
        weight <- length(intersect(newlist[[i]],newlist[[j]]))
        labeli <- getSymbol(names(newlist[i]),chip)
        labelj <- getSymbol(names(newlist[j]),chip)
        genes <- c(genes,labeli,labelj)
        if (is.na(labeli))
        {
          labeli <- names(newlist[i])
        }
        if (is.na(labelj))
        {
          labelj <- names(newlist[j])
        }
        if (labeli != labelj)
        {  
          cat(labeli,"(pp)",labelj, "=", weight, file=paste("output_",comparison,"_",intersection,"intersections.weight",sep=""),"\n", append = TRUE)
          cat(labeli,"pp",labelj, file=paste("output_",comparison,"_",intersection,"intersections.sif",sep=""),"\n", append = TRUE)
        }  
      }  
    } 
  }

  genes <- unique(genes)
  write.table(genes,file=paste(intersection,"intersections_genes_",comparison,".txt",sep=""),quote=F,col.names=F,row.names=F)  
  if (counter == 0)
  {
    print("No Gene satisfies intersection criteria.")
  }
  else
  {
    # run Cytoscape
#    system("java -jar C:/Programme/Cytoscape-v3.2.0/cytoscape.jar -i output.sif")

    #After cytoscape is closed manually by the user, the output files are deleted.
#    system("rm output.sif")
#    system("rm output.weight")
  }
}
