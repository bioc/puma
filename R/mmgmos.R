mmgmos <- function(
	object
,	background=FALSE
,	replaceZeroIntensities=TRUE
,	gsnorm=c("median", "none", "mean", "meanlog")
,	savepar=FALSE
,	eps=1.0e-6	
,	addConstant = 0
)
{

if (class(object)[1]=='ExpressionFeatureSet'){           # if the user use the oligo package loading the cel files
   
   
          conds <- length(sampleNames(object));
	genes <-length(unique(oligo:::probeNames(object)));
	phis <- c(0,0,0);

          pm_g<-oligo:::pm(object);
          mm_g<-oligo:::mm(object);
          probe<-(oligo:::probeNames(object));   

          index<-order(probe);
          probe_sort<-probe[index];
          probe<-probe_sort;
          pm_g<-pm_g[index,];
          mm_g<-mm_g[index,];


       if(conds==1){
          pm_g<-as.matrix(pm_g);
          mm_g<-as.matrix(mm_g);
       }

           
      if (background == TRUE)
      {
         for (i in c(1:conds)){
             m<-min(c(min(pm_g[,i]),min(mm_g[,i])))
             pm_g[,i]<-pm_g[,i]-m+1
             mm_g[,i]<-mm_g[,i]-m+1
         }
      }

      if (replaceZeroIntensities)
     {
         pm_g[which(pm_g==0)] <- 1
         mm_g[which(mm_g==0)] <- 1
     }

      prctiles <- 0.01*c(5, 25, 50, 75, 95);

  res <-
  	.Call(
  	  "mmgmos_c"
  	 , pm_g
  	 , mm_g
  	 , genes
  	 , probe
  	 , phis
  	 , prctiles
  	 , length(prctiles)
  	 , savepar
  	 , eps
  	 , PACKAGE="puma"
  	 )

}else if (class(object)[1]=='AffyBatch')             ##if the user use the affy package loading the cel files 
{
        probes <- length(affy:::probeNames(object))
	conds <- length(object)
	genes <- length(featureNames(object))
        orig.phis = FALSE
	cdf <- cleancdfname(cdfName(object))
	phiname <- paste(substr(cdf,1,nchar(cdf)-3), "phis", sep="")
	if(phiname %in% do.call("data", list(package="puma"))$results[, 3])
	{
	    do.call("data", list(phiname))
	    phis <- eval(parse(text=phiname))
	}
	else
	{
		if(orig.phis)
		{
			do.call("data", list("hgu95aphis"))
			phis <- eval(parse(text="hgu95aphis"))
		}
		else
			phis <- c(0,0,0)
        }
        prctiles <- 0.01*c(5, 25, 50, 75, 95);
  
       if (background == TRUE)
       {
          for (i in c(1:conds)){
             m<-min(c(min(affy:::pm(object)[,i]),min(affy:::mm(object)[,i])))
             affy:::pm(object)[,i]<-affy:::pm(object)[,i]-m+1
             affy:::mm(object)[,i]<-affy:::mm(object)[,i]-m+1
           }
       }

     if (replaceZeroIntensities)
     {
        affy:::pm(object)[which(affy:::pm(object)==0)] <- 1
        affy:::mm(object)[which(affy:::mm(object)==0)] <- 1
     }


    res <-
  	.Call(
  	  "mmgmos_c"
  	 , affy:::pm(object)
  	 , affy:::mm(object)
  	 , genes
  	 , affy:::probeNames(object)
  	 , phis
  	 , prctiles
  	 , length(prctiles)
  	 , savepar
  	 , eps
  	 , PACKAGE="puma"
  	 )

}        


  expr <- matrix(res[c(1:genes),],genes,conds)
  se <- matrix(res[c((genes+1):(2*genes)),],genes,conds)
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

if (class(object)[1]=='ExpressionFeatureSet'){


  probe_names<-unique(probe);
  rownames(expr) <-probe_names
  colnames(expr) <- sampleNames(object)
  rownames(se) <-probe_names
  colnames(se) <- sampleNames(object)
  rownames(prc5) <-probe_names
  colnames(prc5) <- sampleNames(object)
  rownames(prc25) <-probe_names
  colnames(prc25) <- sampleNames(object)
  rownames(prc50) <-probe_names
  colnames(prc50) <- sampleNames(object)
  rownames(prc75) <-probe_names
  colnames(prc75) <- sampleNames(object)
  rownames(prc95) <-probe_names
  colnames(prc95) <- sampleNames(object)
}else if (class(object)[1]=='AffyBatch')
{
    rownames(expr) <- featureNames(object)
  colnames(expr) <- sampleNames(object)
  rownames(se) <- featureNames(object)
  colnames(se) <- sampleNames(object)
  rownames(prc5) <- featureNames(object)
  colnames(prc5) <- sampleNames(object)
  rownames(prc25) <- featureNames(object)
  colnames(prc25) <- sampleNames(object)
  rownames(prc50) <- featureNames(object)
  colnames(prc50) <- sampleNames(object)
  rownames(prc75) <- featureNames(object)
  colnames(prc75) <- sampleNames(object)
  rownames(prc95) <- featureNames(object)
  colnames(prc95) <- sampleNames(object)
}
  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)
  # notes <- notes(object)

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
	# phenoData(return_exprReslt) <- phenoData(object)
	# pData(return_exprReslt) <- pData(object)
	annotation(return_exprReslt) <- annotation(object)
	description(return_exprReslt) <- description(object)
	notes(return_exprReslt) <- notes(object)
	return_exprReslt;
}
