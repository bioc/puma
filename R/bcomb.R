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
