"pumaClustii" <-function(e=NULL, se=NULL, efile=NULL, sefile=NULL, 
                       subset=NULL, gsnorm=FALSE, mincls, maxcls, 
                       conds, reps,  
                       verbose=FALSE,
                       eps=1.0e-5, del0=0.01)
{
  if (length(e)==0)
  {
    e <- read.csv(efile,check.names=FALSE,row.names=1)
  }
  
  dim_e <- dim(e)
  genes <- dim_e[1]
  chips <- dim_e[2]

  if (length(se)==0 & length(sefile)>0)
  {
    se <- read.csv(sefile,check.names=FALSE,row.names=1)
    se <- se^2
    m <- 1; 
  } 
  else
  {
    if (length(se)==0 & length(sefile)==0)
    {
      se <- matrix(0, genes, chips)
      m <- 0;  
    }
    else
    {
      se <- se^2
      m <- 1;
    }
  }

  if(gsnorm==TRUE)
    e<-normalisation.gs(e)
 
  if (length(subset)>0)
  {
    e <- e[subset,]
    se<- se[subset,]
    dim_e <- dim(e)
    genes <- dim_e[1]
    chips <- dim_e[2]
  }

  if (m==1)
  {
    se <- t(apply(cbind(e,se), 1, clusterNormVar))
  }
  e <- t(apply(e, 1, clusterNormE))
  ee<-matrix(0,genes,conds)
  for (i in 1:conds)
  {
  	ind<-which(reps==i)
  	if (length(ind)==1)
  	  ee[,i]<-e[,i]
  	else
      ee[,i]<-t(apply(e[,ind],1,mean))
  }

  
  cl<-Mclust(ee,G=maxcls) 
  centers <- matrix(0,maxcls,conds)
  clsig <- matrix(0,maxcls,conds)
  for (i in 1:maxcls)
  {
    if (length(which(cl$classification==i))<1)
       stop("please specify a smaller maxcls")
       
    centers[i,] <- apply(ee[which(cl$classification==i),],2,mean)
    clsig[i,]<-diag(cov(ee[which(cl$classification==i),]))

  }

  res<-.Call("pumaclustii_c", as.matrix(e), as.matrix(se),conds,reps, mincls, maxcls, as.matrix(centers), as.matrix(clsig), verbose, eps, del0, PACKAGE="puma")
  
  names(res[[1]]) <- rownames(e)
  colnames(res[[2]]) <- paste(1:conds)
  rownames(res[[2]]) <- paste(1:res[[5]])
  colnames(res[[3]]) <- paste(1:conds)
  rownames(res[[3]]) <- paste(1:res[[5]])
  rownames(res[[4]]) <- rownames(e)
  colnames(res[[4]]) <- paste(1:res[[5]])

  return(list(cluster=res[[1]], centers=res[[2]], centersigs=res[[3]], likelipergene=res[[4]], optK=res[[5]], optF=res[[6]]))
}
