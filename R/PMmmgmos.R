PMmmgmos<-function(
	object
,	background=TRUE
,	replaceZeroIntensities=TRUE
,	gsnorm=c("median", "none", "mean", "meanlog")
,	savepar=FALSE
,	eps=1.0e-6
,        addConstant = 0
)
{   
    
       
 	conds <- length(sampleNames(object));
	genes <-length(unique(oligo:::probeNames(object)));
    
   
          probe<-(oligo:::probeNames(object));        
          pm_g<-oligo:::pm(object);
      
          index<-order(probe);
          probe_sort<-probe[index];
          probe<-probe_sort;
          pm_g<-pm_g[index,];
   

  if(conds==1){
    pm_g<-as.matrix(pm_g);
  }
        
	

  if(background==TRUE)
  {
     for (i in c(1:conds))
    {
      m<-min(pm_g[,i])
      pm_g[,i]<-pm_g[,i]-m+1
     
    }
  }

  if (replaceZeroIntensities)
  {
    pm_g[which(pm_g==0)] <- 1
   
  }
  prctiles <- 0.01*c(5, 25, 50, 75, 95);
 
  res <-
 	.Call(
  	  "pmmmgmos_c"
  	 , pm_g
  	 , genes
         , probe
         , prctiles
  	 , length(prctiles)
   	 , savepar
  	 , eps
  	 , PACKAGE="puma"
  	 )
     
     expr <- matrix(res[c(1:genes),],genes,conds);
     se <- matrix(res[c((genes+1):(2*genes)),],genes,conds);
     prc5 <- matrix(res[c((2*genes+1):(3*genes)),],genes,conds)
     prc25 <- matrix(res[c((3*genes+1):(4*genes)),],genes,conds)
     prc50 <- matrix(res[c((4*genes+1):(5*genes)),],genes,conds)
     prc75 <- matrix(res[c((5*genes+1):(6*genes)),],genes,conds)
     prc95 <- matrix(res[c((6*genes+1):(7*genes)),],genes,conds)
  
    rm(res)
   if (gsnorm[1]=="mean")
   {
     expr <- as.data.frame(2^expr)
    
     chipm <- apply(expr,2,mean)
     chipm <- chipm/chipm[1]

     expr <- as.matrix(log2(expr))
   
    for (i in 1:conds)
    {
      expr[,i] <- expr[,i]-log2(chipm[i])
      prc5[,i] <- prc5[,i]-log2(chipm[i])
      prc25[,i] <- prc25[,i]-log2(chipm[i])
      prc50[,i] <- prc50[,i]-log2(chipm[i])
      prc75[,i] <- prc75[,i]-log2(chipm[i])
      prc95[,i] <- prc95[,i]-log2(chipm[i])
    }
  }
  else if (gsnorm[1]=="median")
  {
     expr <- as.data.frame(2^expr)
    
     chipm <- apply(expr,2,median)
     chipm <- chipm/chipm[1]

     expr <- as.matrix(log2(expr))
    for (i in 1:conds)
    {
      expr[,i] <- expr[,i]-log2(chipm[i])
      prc5[,i] <- prc5[,i]-log2(chipm[i])
      prc25[,i] <- prc25[,i]-log2(chipm[i])
      prc50[,i] <- prc50[,i]-log2(chipm[i])
      prc75[,i] <- prc75[,i]-log2(chipm[i])
      prc95[,i] <- prc95[,i]-log2(chipm[i])
    }
  }
  else if (gsnorm[1]=="meanlog")
  {
    chipm <- apply(expr,2,mean)
    chipm <- chipm-chipm[1]

    for (i in 1:conds)
    {
      expr[,i] <- expr[,i]-chipm[i]
      prc5[,i] <- prc5[,i]-chipm[i]
      prc25[,i] <- prc25[,i]-chipm[i]
      prc50[,i] <- prc50[,i]-chipm[i]
      prc75[,i] <- prc75[,i]-chipm[i]
      prc95[,i] <- prc95[,i]-chipm[i]
    }
  }
   
   probe_names<-unique(probe);
   rownames(expr) <- probe_names
   colnames(expr) <- sampleNames(object)
   rownames(se) <- probe_names
   colnames(se) <- sampleNames(object)
   rownames(prc5) <- probe_names
   colnames(prc5) <- sampleNames(object)
   rownames(prc25) <- probe_names
   colnames(prc25) <- sampleNames(object)
   rownames(prc50) <- probe_names
   colnames(prc50) <- sampleNames(object)
   rownames(prc75) <- probe_names
   colnames(prc75) <- sampleNames(object)
   rownames(prc95) <- probe_names
   colnames(prc95) <- sampleNames(object)


    phenodata <- phenoData(object)
    annotation <- annotation(object)
    escription <- description(object)
  #notes <- notes(object)
	
    return_exprReslt <- new(
		"exprReslt"
	,	exprs=log2((2^expr)+addConstant)
	,	se.exprs=se
 	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData(object)
		,	varMetadata=data.frame(labelDescription=varLabels(phenoData(object)))
		)
	# , notes=notes
	)
       
          prcfive(return_exprReslt) <- prc5
	prctwfive(return_exprReslt) <- prc25
          
          prcfifty(return_exprReslt) <- prc50
	prcsevfive(return_exprReslt) <- prc75
	prcninfive(return_exprReslt) <- prc95
	#phenoData(return_exprReslt) <- phenoData(object)
	#pData(return_exprReslt) <- pData(object)
    
     annotation(return_exprReslt) <- annotation(object)
     description(return_exprReslt) <- description(object)
     notes(return_exprReslt) <- notes(object)
     rm(object);
     return_exprReslt
     
     ##write.table(se,'ddd.txt',quote=FALSE,sep='\t',row.names=FALSE,col.name=FALSE);

}
