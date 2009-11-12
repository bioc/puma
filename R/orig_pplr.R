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
