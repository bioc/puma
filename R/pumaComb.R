clusterApplyLBDots <- function (cl, x, fun, ...) 
{
    chunksCount <- 0
    numOfChunks <- length(x)
    numOfNodes <- length(cl)
    #fun <- fun
    if (numOfChunks > 0 && numOfNodes > 0) {
        wrap <- function(x, i, ...) list(value = try(fun(x, ...)), 
            index = i)
        submit <- function(node, job) {
            args <- c(list(x[[job]]), list(job), list(...))
            sendCall(cl[[node]], wrap, args)
        }
        for (i in 1:min(numOfChunks, numOfNodes)) submit(i, i)
        val <- vector("list", length(x))
        for (i in seq(along = x)) {
			chunksCount <- chunksCount + 1
			dotsPerChunk <- floor( (50/numOfChunks)*chunksCount ) - floor( (50/numOfChunks)*(chunksCount-1) )
			if(dotsPerChunk > 0) for(dotsCount in 1:dotsPerChunk) cat(".")
            d <- recvOneResult(cl)
            j <- i + min(numOfChunks, numOfNodes)
            if (j <= numOfChunks) 
                submit(d$node, j)
            val[d$value$index] <- list(d$value$value)
        }
        val
    }
}

pumaComb <- function (
	eset
,	design.matrix=NULL
,	method="em"
,	numOfChunks=1000
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
		,	method=method
		)
	}
	cat("Calculating expected completion time\n")
	if(is.null(cl))
	{
		expectedCompletionTime <-	system.time(
										rFirst <- list(
											bcomb(
												paramsList[[1]][[1]]
											,	paramsList[[1]][[2]]
											,	paramsList[[1]][[3]]
											,	paramsList[[1]][[4]]
											)
										)
									)[3] * numOfChunks
	}
	else
	{
		expectedCompletionTime <- system.time(
									rFirst <- clusterApplyLB(
										cl
									,	paramsList[1:length(cl)]
									,	function(paramsList) 
										{
											b <- bcomb(
												paramsList[[1]]
											,	paramsList[[2]]
											,	paramsList[[3]]
											,	paramsList[[4]]
											)
										}
									)
								)[3] * numOfChunks / length(cl)
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
    				b <- bcomb(
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
							b <- bcomb(
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
}

removeUninformativeFactors <- function(eset)
{
	informativeFactors <- which(
		apply(
			as.matrix(pData(eset))
		,	2
		,	function(x) length(levels(as.factor(x)))
		) > 1
	)
	phenoData(eset) <- phenoData(eset)[,informativeFactors]
	return(eset)
}

numOfFactorsToUse <- function (eset)
{
	numOfFactors <- dim(pData(eset))[2]
	if(numOfFactors == 0)
		stop("eset has no informative factors! Does eset have an appropriate phenoData slot?")
	while (
		length(
			levels(
				as.factor(
					apply(
						as.matrix(pData(eset)[,1:numOfFactors])
						,	1
						,	function(items) paste(items, collapse=".")
					)
				)
			)
		)
		<
		prod(
			apply(
				as.matrix(pData(eset)[,1:numOfFactors])
			,	2
			,	function(x) length(levels(as.factor(x)))
			)
		)
	) numOfFactors <- numOfFactors - 1
	return(numOfFactors)
}

createDesignMatrix <- function (eset)
{
	if(!(is(eset, "ExpressionSet") || is(eset, "ExpressionSetIllumina")))
		stop("eset is not a valid ExpressionSet object!")
	numOfFactors <- numOfFactorsToUse(removeUninformativeFactors(eset))
	if(numOfFactors == 0)
		stop("eset has no informative factors! Does eset have an appropriate phenoData slot?")
	pd <- as.data.frame(pData(removeUninformativeFactors(eset))[,1:numOfFactors])
	for (i in 1:numOfFactors)
	## coerce all columns of pheno data to factors - needed because
	## read.phenoData doesn't treat numerics as factors
	{
		pd[,i] <- factor(pd[,i])
	}
	# levelsOfFactors <- rev(lapply(pd,levels))
	levelsOfFactors <- lapply(pd,levels)
	lengthsOfLevels <- sapply(levelsOfFactors,length)

#	levelsOfFactors <- rep(0,numOfFactors)
#	for(i in 1:numOfFactors)
#		levelsOfFactors[i] <- length(levels(pd[,i]))
#	for(firstFactor in 1:levelsOfFactors[1])
#	levelsOfFactors <- rev(levelsOfFactors)
	##create an n dimensional array where n is number of factors. Each element
	##this array is a name combining the levels of each factor

	arrayOfNames <- levelsOfFactors[[1]]
	if(numOfFactors > 1)
	{
		for(i in 2:numOfFactors)
		{
			arrayOfNames <- sapply(
				arrayOfNames
			,	function(x) paste(x,levelsOfFactors[[i]],sep=".")
			)
		}
	}
	arrayOfNames <- aperm(array(arrayOfNames, rev(lengthsOfLevels)))
	# arrayOfNames <- aperm(array(colnames(design), rev(lengthsOfLevels)))
	# vectorOfNames <- levelsOfFactors[[1]]
	# if(numOfFactors > 1)
	# {
	# 	for(i in 2:numOfFactors)
	# 	{
	# 		vectorOfNames <- sapply(
	# 			vectorOfNames
	# 		,	function(x) paste(x,levelsOfFactors[[i]],sep=".")
	# 		)
	# 	}
	# }
	vectorOfNames <- as.character(as.vector(arrayOfNames))
	conditionLevels <-	factor(
							apply(
								as.matrix(pd[,1:numOfFactors])
							,	1
							,	function(items) paste(items, collapse=".")
							)
						,	levels = vectorOfNames
						# ,	levels = unique(
						# 		apply(
						# 			as.matrix(pData(eset)[,1:numOfFactors])
						# 		,	1
						# 		,	function(items) paste(items, collapse=".")
						# 		)
						# 	)
						)
	design.matrix <- model.matrix(~0+conditionLevels)
	colnames(design.matrix) <- vectorOfNames
	# colnames(design.matrix) <- unique(conditionLevels)
	rownames(design.matrix) <- sampleNames(eset)
	return(design.matrix)
}

createContrastMatrix <- function (eset, design=NULL)
{
	if(!(is(eset, "ExpressionSet") || is(eset, "exprSet") || is(eset, "ExpressionSetIllumina")))
		stop("eset is not a valid ExpressionSet object!")
	if(is.null(design))
		design <- createDesignMatrix(eset)
	numOfFactors <- numOfFactorsToUse(removeUninformativeFactors(eset))
	pd <- as.data.frame(pData(removeUninformativeFactors(eset))[,1:numOfFactors])
	for (i in 1:numOfFactors)
	## coerce all columns of pheno data to factors - needed because
	## read.phenoData doesn't treat numerics as factors
	{
		pd[,i] <- as.factor(pd[,i])
	}
	# levelsOfFactors <- rev(lapply(pd,levels))
	levelsOfFactors <- lapply(pd,levels)
	lengthsOfLevels <- sapply(levelsOfFactors,length)
	
#	levelsOfFactors <- rep(0,numOfFactors)
#	for(i in 1:numOfFactors)
#		levelsOfFactors[i] <- length(levels(pd[,i]))
#	for(firstFactor in 1:levelsOfFactors[1])
#	levelsOfFactors <- rev(levelsOfFactors)
	##create an n dimensional array where n is number of factors. Each element
	##this array is a name combining the levels of each factor
	
	# arrayOfNames <- aperm(array(colnames(design), rev(lengthsOfLevels)))
	arrayOfNames <- levelsOfFactors[[1]]
	if(numOfFactors > 1)
	{
		for(i in 2:numOfFactors)
		{
			arrayOfNames <- sapply(
				arrayOfNames
				,	function(x) paste(x,levelsOfFactors[[i]],sep=".")
			)
		}
	}
	arrayOfNames <- aperm(array(arrayOfNames, rev(lengthsOfLevels)))
	
	############################################################################
	## 1-1 contrasts where only one factor is changing for fixed levels of all other
	############################################################################
	
	## create a 2n dimensional array each element of which corrresponds to a
	## combination of two sets of levels. The idea is we then loop through this
	## array looking for combinations which are of interest, i.e. ones for which
	## only one dimension has changed between the two sets, and which have not
	## already been included in the "reverse" direction
	comp_array <- array(1:prod(lengthsOfLevels)^2,dim=c(lengthsOfLevels, lengthsOfLevels))
	
	## find which elements (contrasts) of comp_array are of interest
	cont_vector <- vector()
	for(i in 1:length(comp_array))
	{
		if(
			length(
				which(
					which(comp_array==i, arr.ind=TRUE)[1:numOfFactors]
						==	which(comp_array==i, arr.ind=TRUE)[(numOfFactors+1):(numOfFactors*2)]
				)
			) == (numOfFactors - 1)
			&&
			length(
				which(
					which(comp_array==i, arr.ind=TRUE)[1:numOfFactors]
						<	which(comp_array==i, arr.ind=TRUE)[(numOfFactors+1):(numOfFactors*2)]
				)
			) == 0
			)
			cont_vector <- c(cont_vector,i)
	}
	
	## get dimensions for each contrast of interest, where first n dimensions
	## specify the "experiment" condition, and second n dimensions specify the
	## "control" condition
	cl <- lapply(cont_vector, function(x) which(comp_array == x, arr.ind=TRUE))
	## create a list, each item of which is a contrast, specified as an n-
	## dimensional array
	ct <- list()
	for(i in 1:length(cl))
	{
		ct[[i]] <- array(0,dim=lengthsOfLevels)
		ct[[i]][t(as.matrix(cl[[i]][1:numOfFactors]))] <- 1
		ct[[i]][t(as.matrix(cl[[i]][(numOfFactors+1):(numOfFactors*2)]))] <- -1
	}
	
	## make a matrix out of the list
	cm <- sapply(ct, as.vector)
	colnames(cm) <- sapply(
		cl
		,	function(x) paste(
				arrayOfNames[t(as.matrix(x[,1:numOfFactors]))]
				,	arrayOfNames[t(as.matrix(x[,(numOfFactors+1):(numOfFactors*2)]))]
				,	sep="_vs_"
			)
	)
	rownames(cm) <- arrayOfNames
	
	############################################################################
	## 1-other contrasts where only one factor is changing for fixed levels of all other
	############################################################################
	
	if(numOfFactors == 1)
	{
		ct <- list()
		count <- 0
		contrastNames <- vector()
		if(lengthsOfLevels[1] > 2)
		{
			cm8 <- matrix(0, lengthsOfLevels[1], lengthsOfLevels[1])
			for(k in 1:lengthsOfLevels[1])
			{
				cm8[k,k] <- 1
				cm8[(1:lengthsOfLevels[1])[-k],k] <- -1
			}
			colnames(cm8) <- paste(levelsOfFactors[[1]], "_vs_others", sep="")
			cm <- cbind(cm,cm8)
		}
	}
	
	if(numOfFactors > 1)
	{
		for(index_Factor in 1:length(lengthsOfLevels))
		{
			if(lengthsOfLevels[index_Factor] > 2)
			{
				temp <- array(0,dim=lengthsOfLevels[-index_Factor])
				perms <- which(temp==0, arr.ind=TRUE)
				ct <- list()
				count <- 0
				contrastNames <- vector()
				for(index_unchangingLevels in 1:(dim(perms)[1]))
				{
					for(index_changingLevels in 1:lengthsOfLevels[index_Factor])
					{
						count <- count+1
						ct[[count]] <- array(0,dim=lengthsOfLevels)
						dims<-vector(length=length(lengthsOfLevels))
						dims[-index_Factor] <- perms[index_unchangingLevels,]
						dims[index_Factor] <- index_changingLevels
						dims_for_names <- dims
						dimsmat <- t(matrix(dims))
						ct[[count]][dimsmat] <- 1
						dimsmat <- matrix(0,lengthsOfLevels[index_Factor],length(lengthsOfLevels))
						count2 <- 0
						for(index_changingFactor in
							(1:lengthsOfLevels[index_Factor])[-index_changingLevels])
						{
							dims[index_Factor] <- index_changingFactor
							dimsmat[index_changingFactor,] <- t(matrix(dims))
						}
						ct[[count]][dimsmat] <- -1
						if(index_Factor == 1)
							contrastNames[count] <- levelsOfFactors[[1]][dims_for_names[1]]
						else
							contrastNames[count] <- levelsOfFactors[[1]][dims_for_names[1]]
						for(m in 2:numOfFactors)
						{
							if(index_Factor == m)
								contrastNames[count] <- paste(contrastNames[count], "."
									,  levelsOfFactors[[m]][dims_for_names[m]], sep="")
							else
								contrastNames[count] <- paste(contrastNames[count], "."
									,  levelsOfFactors[[m]][dims_for_names[m]], sep="")
						}
						contrastNames[count] <- paste(contrastNames[count], "_vs_", sep="")
						if(index_Factor == 1)
							contrastNames[count] <- paste(contrastNames[count]
								, "others", sep="")
						else
							contrastNames[count] <- paste(contrastNames[count]
								, levelsOfFactors[[1]][dims_for_names[1]], sep="")
						for(m in 2:numOfFactors)
						{
							if(index_Factor == m)
								contrastNames[count] <- paste(contrastNames[count]
									, ".",  "others", sep="")
							else
								contrastNames[count] <- paste(contrastNames[count]
									, ".",  levelsOfFactors[[m]][dims_for_names[m]], sep="")
						}
					}
				}
				## make a matrix out of the list
				cm7 <- sapply(ct, as.vector)
				colnames(cm7) <- contrastNames
				cm <- cbind(cm,cm7)
			}
		}
	}
	
	
	############################################################################
	## contrasts where only one factor is changing for all levels of all other
	############################################################################
	
	if(numOfFactors > 1)
	{
		i <- 0
		ct <- list()
		contrastNames <- vector()
		for(factorIndex in 1:numOfFactors)
		{
			for(firstLevel in 1:(lengthsOfLevels[factorIndex]-1))
			{
				for(secondLevel in (firstLevel+1):lengthsOfLevels[factorIndex])
				{
					i <- i+1
					ct[[i]] <- array(0,dim=lengthsOfLevels)
					ct[[i]][which(slice.index(ct[[i]],factorIndex)==firstLevel)] <- 1
					ct[[i]][which(slice.index(ct[[i]],factorIndex)==secondLevel)] <- -1
					contrastNames[i] <- paste(
						names(levelsOfFactors)[factorIndex]
						,	"_"
						,	levelsOfFactors[[factorIndex]][firstLevel]
						,	"_vs_"
						,	levelsOfFactors[[factorIndex]][secondLevel]
						,	sep=""
					)
				}
			}
		}
		cm2 <- sapply(ct, as.vector)
		colnames(cm2) <- contrastNames
		cm <- cbind(cm,cm2)
	}
	
	############################################################################
	## 2-way interaction terms when only two factors
	############################################################################
	
	if(numOfFactors == 2)
	{
		i <- 0
		ct <- list()
		contrastNames <- vector()
		for(factor1 in 1:(numOfFactors-1))
		{
			for(factor2 in (factor1+1):numOfFactors)
			{
				for(factor1level1 in 1:(lengthsOfLevels[factor1]-1))
				{
					for(factor1level2 in (factor1level1+1):lengthsOfLevels[factor1])
					{
						for(factor2level1 in 1:(lengthsOfLevels[factor2]-1))
						{
							for(factor2level2 in (factor2level1+1):lengthsOfLevels[factor2])
							{
								i <- i+1
								ct[[i]] <- array(0,dim=lengthsOfLevels)
								ct[[i]][factor1level1,factor2level1] <- 1
								ct[[i]][factor1level1,factor2level2] <- -1
								ct[[i]][factor1level2,factor2level1] <- -1
								ct[[i]][factor1level2,factor2level2] <- 1
								contrastNames[i] <- paste(
									"Int__"
									,	names(levelsOfFactors)[factor1]
									,	"_"
									,	levelsOfFactors[[factor1]][factor1level1]
									,	"."
									,	levelsOfFactors[[factor1]][factor1level2]
									,	"_vs_"
									,	names(levelsOfFactors)[factor2]
									,	"_"
									,	levelsOfFactors[[factor2]][factor2level1]
									,	"."
									,	levelsOfFactors[[factor2]][factor2level2]
									,	sep=""
								)
							}
						}
					}
				}
			}
		}
		cm3 <- sapply(ct, as.vector)
		colnames(cm3) <- contrastNames
		cm <- cbind(cm, cm3)
	}
	
	############################################################################
	## 2-way interaction terms for each level of third factor
	############################################################################
	
	if(numOfFactors == 3)
	{
		i <- 0
		ct <- list()
		contrastNames <- vector()
		for(factor1 in 1:(numOfFactors-1))
		{
			for(factor2 in (factor1+1):numOfFactors)
			{
				for(factor1level1 in 1:(lengthsOfLevels[factor1]-1))
				{
					for(factor1level2 in (factor1level1+1):lengthsOfLevels[factor1])
					{
						for(factor2level1 in 1:(lengthsOfLevels[factor2]-1))
						{
							for(factor2level2 in (factor2level1+1):lengthsOfLevels[factor2])
							{
								factor3 <- (1:numOfFactors)[-c(factor1, factor2)]
								for(factor3level in 1:lengthsOfLevels[factor3])
								{
									i <- i+1
									ct[[i]] <- array(0,dim=lengthsOfLevels)
									arrayDims <- array(0,dim=c(1,3))
									arrayDims[1,factor1] <- factor1level1
									arrayDims[1,factor2] <- factor2level1
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- 1
									arrayDims[1,factor1] <- factor1level1
									arrayDims[1,factor2] <- factor2level2
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- -1
									arrayDims[1,factor1] <- factor1level2
									arrayDims[1,factor2] <- factor2level1
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- -1
									arrayDims[1,factor1] <- factor1level2
									arrayDims[1,factor2] <- factor2level2
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- 1
									contrastNames[i] <- paste(
										"Int__"
										,	names(levelsOfFactors)[factor3]
										,	"_"
										,	levelsOfFactors[[factor3]][factor3level]
										,	"__"
										,	names(levelsOfFactors)[factor1]
										,	"_"
										,	levelsOfFactors[[factor1]][factor1level1]
										,	"."
										,	levelsOfFactors[[factor1]][factor1level2]
										,	"_vs_"
										,	names(levelsOfFactors)[factor2]
										,	"_"
										,	levelsOfFactors[[factor2]][factor2level1]
										,	"."
										,	levelsOfFactors[[factor2]][factor2level2]
										,	sep=""
									)
								}
							}
						}
					}
				}
			}
		}
		cm4 <- sapply(ct, as.vector)
		colnames(cm4) <- contrastNames
		cm <- cbind(cm, cm4)
	}
	
	############################################################################
	## 2-way interaction terms for all levels of third factor
	############################################################################
	
	if(numOfFactors == 3)
	{
		i <- 0
		ct <- list()
		contrastNames <- vector()
		for(factor1 in 1:(numOfFactors-1))
		{
			for(factor2 in (factor1+1):numOfFactors)
			{
				for(factor1level1 in 1:(lengthsOfLevels[factor1]-1))
				{
					for(factor1level2 in (factor1level1+1):lengthsOfLevels[factor1])
					{
						for(factor2level1 in 1:(lengthsOfLevels[factor2]-1))
						{
							for(factor2level2 in (factor2level1+1):lengthsOfLevels[factor2])
							{
								factor3 <- (1:numOfFactors)[-c(factor1, factor2)]
								i <- i+1
								ct[[i]] <- array(0,dim=lengthsOfLevels)
								arrayDims <- array(0,dim=c(1,3))
								for(factor3level in 1:lengthsOfLevels[factor3])
								{
									arrayDims[1,factor1] <- factor1level1
									arrayDims[1,factor2] <- factor2level1
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- 1
									arrayDims[1,factor1] <- factor1level1
									arrayDims[1,factor2] <- factor2level2
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- -1
									arrayDims[1,factor1] <- factor1level2
									arrayDims[1,factor2] <- factor2level1
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- -1
									arrayDims[1,factor1] <- factor1level2
									arrayDims[1,factor2] <- factor2level2
									arrayDims[1,factor3] <- factor3level
									ct[[i]][arrayDims] <- 1
								}
								contrastNames[i] <- paste(
									"Int__"
									,	names(levelsOfFactors)[factor3]
									,	":all-levels"
									,	"__"
									,	names(levelsOfFactors)[factor1]
									,	"_"
									,	levelsOfFactors[[factor1]][factor1level1]
									,	"."
									,	levelsOfFactors[[factor1]][factor1level2]
									,	"_vs_"
									,	names(levelsOfFactors)[factor2]
									,	"_"
									,	levelsOfFactors[[factor2]][factor2level1]
									,	"."
									,	levelsOfFactors[[factor2]][factor2level2]
									,	sep=""
								)
							}
						}
					}
				}
			}
		}
		cm5 <- sapply(ct, as.vector)
		colnames(cm5) <- contrastNames
		cm <- cbind(cm, cm5)
	}
	
	############################################################################
	## 3-way interaction terms
	############################################################################
	
	if(numOfFactors == 3)
	{
		i <- 0
		ct <- list()
		contrastNames <- vector()
		for(factor1 in 1:(numOfFactors-2))
		{
			for(factor2 in (factor1+1):(numOfFactors-1))
			{
				for(factor3 in (factor2+1):numOfFactors)
				{
					for(factor1level1 in 1:(lengthsOfLevels[factor1]-1))
					{
						for(factor1level2 in (factor1level1+1):lengthsOfLevels[factor1])
						{
							for(factor2level1 in 1:(lengthsOfLevels[factor2]-1))
							{
								for(factor2level2 in (factor2level1+1):lengthsOfLevels[factor2])
								{
									for(factor3level1 in 1:(lengthsOfLevels[factor3]-1))
									{
										for(factor3level2 in (factor3level1+1):lengthsOfLevels[factor3])
										{
											i <- i+1
											ct[[i]] <- array(0,dim=lengthsOfLevels)
											ct[[i]][factor1level1,factor2level1,factor3level1] <- 1
											ct[[i]][factor1level1,factor2level1,factor3level2] <- -1
											ct[[i]][factor1level1,factor2level2,factor3level1] <- -1
											ct[[i]][factor1level1,factor2level2,factor3level2] <- 1
											ct[[i]][factor1level2,factor2level1,factor3level1] <- -1
											ct[[i]][factor1level2,factor2level1,factor3level2] <- 1
											ct[[i]][factor1level2,factor2level2,factor3level1] <- 1
											ct[[i]][factor1level2,factor2level2,factor3level2] <- -1
											contrastNames[i] <- paste(
												"Int__"
												,	names(levelsOfFactors)[factor1]
												,	"_"
												,	levelsOfFactors[[factor1]][factor1level1]
												,	"."
												,	levelsOfFactors[[factor1]][factor1level2]
												,	"_vs_"
												,	names(levelsOfFactors)[factor2]
												,	"_"
												,	levelsOfFactors[[factor2]][factor2level1]
												,	"."
												,	levelsOfFactors[[factor2]][factor2level2]
												,	"_vs_"
												,	names(levelsOfFactors)[factor3]
												,	"_"
												,	levelsOfFactors[[factor3]][factor3level1]
												,	"."
												,	levelsOfFactors[[factor3]][factor3level2]
												,	sep=""
											)
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cm6 <- sapply(ct, as.vector)
		colnames(cm6) <- contrastNames
		cm <- cbind(cm, cm6)
	}
	
	return(cm)
}

