
pumaCombImproved <- function (
	eset
,	design.matrix=NULL
,	numOfChunks=1000
,       maxOfIterations=200 
,	save_r=FALSE
,	cl=NULL
,	parallelCompute=FALSE
# ,	parallelCompute=if(
# 		"Rmpi" %in% installed.packages() & "snow" %in% installed.packages()
# 	)
# 	as.logical(length(grep("origin",system("lamnodes",TRUE,TRUE))))
# 	else FALSE
)
{
	if(!(is(eset, "ExpressionSet") || is(eset, "ExpressionSetIllumina")))
		stop("eset is not a valid ExpressionSet object!")
	if(is.null(design.matrix) && is.null(pData(eset)))
		stop("ExpressionSet contains no valid phenoData and no design.matrix
			given!")
	if(parallelCompute && is.null(cl))
	{
		# library(Rmpi)
		library(snow)
		cl <- getMPIcluster()
		if(is.null(cl))
			cl <- makeCluster(length(system("lamnodes",TRUE,TRUE))-1)
		# clusterEvalQ(cl, library(puma))
	}
	if(!is.null(cl))
		clusterEvalQ(cl, library(puma))	
	if(is.null(design.matrix))
		design.matrix <- createDesignMatrix(eset)
	# if(length(names(which(colSums(dm2)==1))) > 1)
	# 	warning("The following conditions have no replicates", names(which(colSums(dm2)==1)),". While it is possible to use PPLR without any biological replicates, this is not advised.")
	e <- exprs(eset)
	se <- assayDataElement(eset,"se.exprs")
	# se <- se.exprs(eset)
	numOfGenes <- dim(e)[1]
	if(numOfChunks > (numOfGenes/10))
		numOfChunks <- max(floor(numOfGenes/10),1)
	numPerChunk <- c(0, rep(floor(numOfGenes / numOfChunks),numOfChunks))
	numOfExtras <- numOfGenes - sum(numPerChunk)
	if(numOfExtras > 0)
		numPerChunk[2:(numOfExtras+1)] <- numPerChunk[2:(numOfExtras+1)] + 1
	eList <- list()
	seList <- list()
	rl <- list()
	rFirst <- list()
	paramsList <- list()
	for (i in 1:numOfChunks)
	{
		eList[[i]] <- e[(1 + sum(numPerChunk[1:i])):sum(numPerChunk[1:(i+1)]),]
		seList[[i]] <- se[(1 + sum(numPerChunk[1:i])):sum(numPerChunk[1:(i+1)]),]
#		seList[[i]] <- se[(1 + (numPerChunk[i] * (i-1))):min((numPerChunk[i+1] * i),numOfGenes),]
		paramsList[[i]] <- list(
			eList[[i]]
		,	seList[[i]]
		,	replicates=apply(design.matrix,1,function(x) which(x==1))
                ,       max_num=maxOfIterations
		)
	}

	if(is.null(cl))
	{
		expectedCompletionTime <-	system.time(
										rFirst <- list(
											hcomb(
												paramsList[[1]][[1]]
											,	paramsList[[1]][[2]]
											,	paramsList[[1]][[3]]
											,	paramsList[[1]][[4]]
											)
										)
									)[3] * numOfChunks
	}else
	{
		expectedCompletionTime <- system.time(
									rFirst <- clusterApplyLB(
										cl
									,	paramsList[1:length(cl)]
									,	function(paramsList) 
										{
											b <- hcomb(
												paramsList[[1]]
											,	paramsList[[2]]
											,	paramsList[[3]]
											,	paramsList[[4]]
											)
										}
									)
								)[3] * numOfChunks /length(cl)
	}
	if (expectedCompletionTime < 120)
		expectedTimeString <- paste(round(expectedCompletionTime, 0), "seconds")
	else if (expectedCompletionTime < 7200)
		expectedTimeString <- paste(round(expectedCompletionTime/60, 0), "minutes")
	else if (expectedCompletionTime < 172800)
		expectedTimeString <- paste(round(expectedCompletionTime/3600, 0), "hours")
	else
		expectedTimeString <- paste(round(expectedCompletionTime/172800, 0), "days")
	cat(paste("pumaComb expected completion time is", expectedTimeString,"\n"))
	cat(".......20%.......40%.......60%.......80%......100%\n")
	i <- 0
	chunksPerDot <- ceiling(numOfChunks / 50)
	if(is.null(cl) & numOfChunks > 1)
    	rl <- lapply(
    			paramsList[2:length(paramsList)]
    		,	function(paramsList) 
    			{
    				b <- hcomb(
						paramsList[[1]]
					,	paramsList[[2]]
					,	paramsList[[3]]
                                        ,	paramsList[[4]]
					)
					i <<- i + 1
					dotsPerChunk <- floor( (50/(numOfChunks-1))*i ) - floor( (50/(numOfChunks-1))*(i-1) )
					#if(i %% chunksPerDot == 0) for(d in 1:dotsPerChunk) cat(".")
					if(dotsPerChunk > 0) for(d in 1:dotsPerChunk) cat(".")
					return(b)
				}
			)
	if(is.null(cl) & numOfChunks == 1)
		cat("..................................................")
    if(!is.null(cl))
    {
    	if(numOfChunks > length(cl))
			rl <-	clusterApplyLBDots(
						cl
					,	paramsList[(length(cl)+1):length(paramsList)]
					,	function(paramsList) 
						{
							b <- hcomb(
								paramsList[[1]]
							,	paramsList[[2]]
							,	paramsList[[3]]
                                                        ,	paramsList[[4]]
							)
						}
					)
    	if(numOfChunks <= length(cl))
			cat("..................................................")
	}
	cat("\n")
	rl <- c(rFirst, rl)
	r <- rl[[1]]
	if(numOfChunks > 1)
		for (i in 2:numOfChunks)
			r <- rbind(r,rl[[i]])
	colnames(r) <- c(paste("M", colnames(design.matrix)), paste("Std", colnames(design.matrix)))
	if(save_r)
		save(r, file="r.rda")
	eset_r <- create_eset_r(eset, r, design.matrix)
	# colnames(exprs(eset_r)) <- paste("M", colnames(design.matrix), sep="_")
	# colnames(assayDataElement(eset_r,"se.exprs")) <- paste("Std", colnames(design.matrix), sep="_")
	# design.matrix.ordering <- replace(
	# 	which(unique(design.matrix)==1)%%(dim(design.matrix)[2])
	# ,	which(unique(design.matrix)==1)%%(dim(design.matrix)[2])==0
	# ,	dim(design.matrix)[2]
	# )
	# design.matrix.ordering <- which(unique(design.matrix)==1, arr.ind=TRUE)[,1][which(unique(design.matrix)==1, arr.ind=TRUE)[,1]]
	# rownames(pData(eset_r)) <- colnames(design.matrix)[design.matrix.ordering]
	# varLabels(eset_r) <- colnames(pData(eset)[1:numOfFactorsToUse(removeUninformativeFactors(eset))])
	# colnames(pData(eset_r)) <- colnames(pData(eset))
	return(eset_r)
             cat('\nDone.')
}
