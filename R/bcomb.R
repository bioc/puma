bcomb <-function(e,se,replicates,method=c("map","em"), 
                 gsnorm=FALSE, nsample=1000, eps=1.0e-6)
{

  #e<-read.csv(efile,check.names=FALSE,row.names=1)
  #se<-read.csv(sefile,check.names=FALSE,row.names=1)
  
  dim_e<-dim(e)
  genes<-dim_e[1]
  chips<-dim_e[2]
  conds<-max(replicates)
  
  if(conds<2){
    cat("The number of conditions is too small. It should be at least two.
      Check the input of replicates.\n")
    return();
  }

  if(method=="em")
    m<-1
  else if(method=="map")
    m<-4
  else{
    cat("Unrecognised method:", method, "\n")
    return();
  }
  
  if(gsnorm==TRUE)
    e<-normalisation.gs(e)
  
  se<-se^2
  
  set.seed(123456)
  res<-.Call(
    "bcomb_c"
  , as.matrix(e)
  , as.matrix(se)
  , replicates
  , m
  , conds
  , nsample
  , eps
  , PACKAGE="puma"
  )
  
  out<-matrix(res, genes, 2*conds)
  
  rownames(out)<-rownames(e)
  
  a<-c("M1")
  b<-c("Std1")
  for(i in 2:conds){
    a<-c(a,paste("M",i,sep=""))
    b<-c(b,paste("Std",i,sep=""))
  }
  colnames(out)<-c(a,b)
  
  return(as.data.frame(out))
}

orig_pplr <-function(e, control, experiment){

  dim_e<-dim(e)
  genes<-dim_e[1]
  conds<-dim_e[2]/2
  
  res<-matrix(0,genes,9)
  res[,2]<-e[,control]
  res[,3]<-e[,experiment]
  res[,4]<-e[,conds+control]
  res[,5]<-e[,conds+experiment]
  res[,6]<-e[,experiment]-e[,control]
  res[,7]<-sqrt(e[,conds+control]^2+e[,conds+experiment]^2)
  res[,8]<--res[,6]/(res[,7]*sqrt(2))
  res[,9]<-erfc(res[,8])/2
  
  sortout<-sort(res[,8],method="quick",index.return=TRUE)
  
  res[,c(2:9)]<-res[sortout$ix,c(2:9)]
  res[,1]<-sortout$ix
  
  colnames(res)<-c("index","cM","sM","cStd","sStd","LRM","LRStd","stat","PPLR")
  rownames(res)<-rownames(e)[sortout$ix]
   
  return(as.data.frame(res)  )
}

erfc <- function(x) 2 * pnorm(x * sqrt(2), lower=FALSE)

normalisation.gs <-function(x)
{
  x<-as.data.frame(2^x)
  d<-dim(x)
  c<-d[2]
  m<-mean(x)
  m<-m/m[1]
  
  for(i in 2:c)
    x[[i]]<-x[[i]]/m[i]
  x<-log2(x)  
  
  return(x)
}
