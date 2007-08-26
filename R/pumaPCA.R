pumaPCA <- function
(
    eset
,   latentDim           = 	if(dim(exprs(eset))[2] <= 3)
 								dim(exprs(eset))[[2]]-1
							else
								3
,   sampleSize          =	if(dim(exprs(eset))[1] <= 1000)
								 dim(exprs(eset))[[1]]
							else
								1000	## Set to integer or FALSE for all
,   initPCA             =	TRUE	## Initialise parameters with PCA
,   randomOrder         =	FALSE	## Update parameters in random order
,   optimMethod         =	"BFGS"  ## ?optim for details of methods
,   stoppingCriterion   =	"deltaW"## can also be "deltaL"
,   tol                 =	1e-3	## Stop when delta update < this
,   stepChecks          =	FALSE	## Check likelihood after each update?
,   iterationNumbers    =	TRUE	## Show iteration numbers?
,   showUpdates         =	FALSE	## Show values after each update?
,   showTimings         =	FALSE	## Show timings after each update?
,   showPlot            =	FALSE	## Show projection plot after each update?
,   maxIters            =	500		## Number of EM iterations.
,   transposeData       =	FALSE	## Transpose eset matrices?
,   returnExpectations  =	FALSE
,   returnData          =	FALSE
,   returnFeedback      =	FALSE
,	pumaNormalize		=	TRUE
)
{
    if(pumaNormalize)
	{
		eset <- pumaNormalize(
			eset
		,	arrayScale="median"
		,	probesetScale="median"
		)
	}
	if(sampleSize)
    {
		sampleY <- sample(dim(exprs(eset))[1], sampleSize)
		Y <- exprs(eset)[sampleY, ]
		varY <- assayDataElement(eset,"se.exprs")[sampleY, ]
    }
    else
    {
		Y <- exprs(eset)
		varY <- assayDataElement(eset,"se.exprs")
    }
    
    if(transposeData)
    {
        Y <- t(Y)
        varY <- t(varY)
    }

    model <- new("pumaPCAModel")
	# model <- list()
    expectations <- new("pumaPCAExpectations")
    # expectations <- list()
    
    numData <- dim(Y)[1]
    dataDim <- dim(Y)[2]
    
    ##  Initialise the model
    if (initPCA)
    {
        avgVarY     <-  colMeans(varY)
        ## scale the rows to transform FA to PCA
        scaledY     <-  t(solve(diag(sqrt(avgVarY)),t(Y)))
        pca         <-  prcomp(scaledY)
        model@sigma <-  sqrt(mean(pca$sdev[(latentDim + 1):dataDim] ^ 2))
        model@m     <-  matrix(0, 1, latentDim)
        model@Cinv  <-  diag(latentDim)
        model@W     <-  diag(sqrt(avgVarY)) %*%
                        (   pca$rotation[, 1:latentDim] %*%
                            diag(
                                sqrt(
                                    pca$sdev[1:latentDim] ^ 2 - model@sigma ^ 2
                                )
                            )
                        )
        ## initialise mean to weighted mean
        weightedY <- matrix(0, numData, dataDim)
        for (i in 1:numData)
        {
            weightedY[i, ] <- varY[i, ] * Y[i, ]
        }
        model@mu <- matrix(
                        solve(
                            diag(colSums(varY))
                        ,   colSums(weightedY)
                        )
                    ,   1
                    ,   dataDim
                    )
    }
    else
    {
    precYPlusSigma <- 1/varY
    model@mu <- matrix(
                    colSums(
                        (Y * (1/varY))
                        /
                        rep(colSums(1/varY),each=numData)
                    )
                ,   1
                ,   dataDim
                )
    model@sigma <-  mean
                    (
                        colSums
                        (
                            (
                                (Y - rep(model@mu, each=numData)) ^ 2
                                *
                                (1 / varY)
                            )
                            /
                            colSums(1 / varY)
                        )
                        /
                        colSums(1 / varY)
                    )
    model@W <- matrix(rnorm(dataDim * latentDim) * 0.1, dataDim, latentDim)
    ##model@W <- matrix(0.1, dataDim, latentDim)
    model@m <- matrix(0, 1, latentDim)
    model@Cinv <- diag(latentDim)
    }
    
    ##  Initialise the expectations
    expectations@x <- matrix(0, numData, latentDim)
    expectations@xxT <- array(0, c(latentDim, latentDim, numData))
    expectations <- pumaPCAEstep(model, expectations, varY, Y)
    
    ##  Compute the starting likelihood.
    maxDeltaL <- 1
    deltaW <- 1
    counter <- 0
    params <- c('mu', 'W', 'sigma', 'estep')
    oldL <- pumaPCALikelihoodBound(model, expectations, varY, Y)
    likelihoodHistory <- list()
    timingHistory <- list()
    modelHistory <- list()
    exitReason <- "unknown exit reason"
    if(stoppingCriterion=="deltaL")
    {
        exitReason <- sprintf(
            'Update of Likelihood less than tolerance %f'
        ,   tol
        )
    }
    else if(stoppingCriterion=="deltaW")
    {
        exitReason <- sprintf(
            'Update of W less than tolerance %f'
        ,   tol
        )
    }
    
    timeToCompute <- system.time(
    {
    tryCatch(
    {
        while (
            (   ( stoppingCriterion == "deltaL" & maxDeltaL > tol ) |
                ( stoppingCriterion == "deltaW" & deltaW > tol )
            )
            &   counter < maxIters
        )
        {
            maxDeltaL <- 0
            counter <- counter + 1
            if(randomOrder) { params <- sample(params) }
            for (param in params)
            {
                if (param == 'mu')
                {
                    model@m <- pumaPCAUpdateM(model, expectations, varY, Y)
                    model@mu <- pumaPCAUpdateMu(model, expectations, varY, Y)
                }
                else if (param == 'W')
                {
                    model@Cinv <-   pumaPCAUpdateCinv(
                                        model
                                    ,   expectations
                                    ,   varY
                                    ,   Y
                                    )
                    model@W <- pumaPCAUpdateW(model, expectations, varY, Y)
                    modelHistory[[counter]] <- model@W
                    if (counter > 1)
                    {
                        deltaW <-   matrixDistance(
                                        model@W[,1:2]
                                    ,   modelHistory[[counter-1]][,1:2]
                                    )
                        if(showUpdates)
							cat(sprintf('Movement in W %s\n', deltaW))
                    }
                
                }
                else if (param == 'sigma')
                {
                    if(
                        optimMethod == "optimise"
                        |
                        optimMethod == "optimize"
                    )
                    {
                         model@sigma <- optimise(
                                            function(sigma)
                                            pumaPCASigmaObjective
                                            (
                                                sigma
                                            ,   model
                                            ,   expectations
                                            ,   varY
                                            ,   Y
                                            )
                                        ,   c(0, model@sigma)
                                        ,   tol=1e-12
                                        )$minimum
                     }
                     else if (optimMethod == "newton")
                     {
                        model@sigma <-  pumaPCANewtonUpdateLogSigma(
                                            model
                                        ,   expectations
                                        ,   varY
                                        ,   Y
                                        )
                     }
                     else if (optimMethod == "0")
                     {
                        model@sigma <- 0
                     }
                     else
                     {
                        optimOutput <- optim(
                                            par=model@sigma
                                        ,   fn=pumaPCASigmaObjective
                                        ,   gr=pumaPCASigmaGradient
                                        ,   method=optimMethod
                                        ,   model=model
                                        ,   expectations=expectations
                                        ,   varY=varY
                                        ,   Y=Y
                                        )
                        model@sigma <- optimOutput$par
                        if(showUpdates)
							cat(optimOutput$counts, "\n")
                    }
                }
                else if (param == 'estep')
                {
                    model <- pumaPCARemoveRedundancy(model)
                    expectations <- pumaPCAEstep(model, expectations, varY, Y)
                }
                if (stepChecks)
                {
    
                    deltaLoldL <-   pumaPCALikelihoodCheck(
                                        model
                                    ,   expectations
                                    ,   varY
                                    ,   Y
                                    ,   oldL
                                    ,   i
									,	stepChecks
                                    )
                    deltaL <- deltaLoldL[[1]]
                    oldL <- deltaLoldL[[2]]
                    maxDeltaL <- max(maxDeltaL, deltaL)
                    likelihoodHistory[[param]][counter] <- deltaL
                }
                timingHistory[[param]][counter] <- date()
                if (showTimings)
                {
                    cat(sprintf('%s after update of %s\n', date(), param))
                }
                if (showPlot)
                {
                    par(mfcol=c(1,2))
                    plot(
                        x   =   model@W[,1]
                    ,   y   =   model@W[,2]
                    ,   pch =   unclass(
							as.factor(
			 					apply(
			 						as.matrix(pData(eset))
			 					,	1
			 					,	function(items) paste(items, collapse=".")
			 					)
			 				)
						)
                    )
                    plot(
                        x   =   model@W[,1]
                    ,   y   =   model@W[,3]
                    ,   pch =   unclass(
							as.factor(
			 					apply(
			 						as.matrix(pData(eset))
			 					,	1
			 					,	function(items) paste(items, collapse=".")
			 					)
			 				)
						)
                    )
                }
            }
            if (!stepChecks)
            {
                deltaLoldL <-   pumaPCALikelihoodCheck(
                                    model
                                ,   expectations
                                ,   varY
                                ,   Y
                                ,   oldL
                                ,   param
								,	stepChecks
                                )
                maxDeltaL <- deltaLoldL[[1]]
                oldL <- deltaLoldL[[2]]
                likelihoodHistory[["estep"]][counter] <- maxDeltaL
            }
    
            if (iterationNumbers)
            {
                cat(sprintf('Iteration number: %d\n', counter))
            }
        }
    }
    ,   interrupt = function(interrupt)
        {
            cat("User interrupt!\n")
            exitReason <<- "User interrupt"
        }
    )
    }
    )[3]
    if (counter >= maxIters)
    {
        cat(sprintf('Warning maximum iterations exceeded.\n'))
        exitReason <- "Iterations exceeded"
    }
    model <- pumaPCARemoveRedundancy(model)
    expectations <- pumaPCAEstep(model, expectations, varY, Y)
	new(
		"pumaPCARes"
    ,   model               = model
    ,   expectations        = if(returnExpectations) expectations
								else new("pumaPCAExpectations")
    ,   varY                = if(returnData) varY else matrix()
    ,   Y                   = if(returnData) Y else matrix()
	,	phenoData			= phenoData(eset)
    ,   timeToCompute       = if(returnFeedback) timeToCompute else 0
    ,   numberOfIterations  = if(returnFeedback) counter else 0
    ,   likelihoodHistory   = if(returnFeedback) likelihoodHistory else list()
    ,   timingHistory       = if(returnFeedback) timingHistory else list()
    ,   modelHistory        = if(returnFeedback) modelHistory else list()
    ,   exitReason          = if(returnFeedback) exitReason else ""
    )
}

