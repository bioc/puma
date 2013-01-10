mmgmos <- function(
	object
,	background=FALSE
,	replaceZeroIntensities=TRUE
,	gsnorm=c("median", "none", "mean", "meanlog")
,	savepar=FALSE
,	eps=1.0e-6
,	orig.phis = FALSE
,	addConstant = 0
)
{

	probes <- length(probeNames(object))
	conds <- length(object)
	genes <- length(featureNames(object))

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
      m<-min(c(min(pm(object)[,i]),min(mm(object)[,i])))
      pm(object)[,i]<-pm(object)[,i]-m+1
      mm(object)[,i]<-mm(object)[,i]-m+1
    }
  }

  if (replaceZeroIntensities)
  {
    pm(object)[which(pm(object)==0)] <- 1
    mm(object)[which(mm(object)==0)] <- 1
  }

  res <-
  	.Call(
  	  "mmgmos_c"
  	 , pm(object)
  	 , mm(object)
  	 , genes
  	 , probeNames(object)
  	 , phis
  	 , prctiles
  	 , length(prctiles)
  	 , savepar
  	 , eps
  	 , PACKAGE="puma"
  	 )

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
	return_exprReslt
}
