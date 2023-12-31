hcomb <-function(e,se,replicates, max_num=c(200,500,1000),gsnorm=FALSE,  eps=1.0e-6)
{

  dim_e<-dim(e)
  genes<-dim_e[1]
  chips<-dim_e[2]
  conds<-max(replicates)
  
  if(conds<2){
    cat("The number of conditions is too small. It should be at least two. Check the input of replicates.\n")
    return();
  }

  if(gsnorm==TRUE)
    e<-normalisation.gs(e)
  
  se<-se^2
         
  set.seed(123456)
  res<-.Call("hcomb_c", as.matrix(e), as.matrix(se), replicates,conds, max_num, eps, PACKAGE="puma")
  
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