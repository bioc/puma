plotErrorBars <- function(
	eset
,	probesets = if(dim(exprs(eset))[1] <= 12) 1:dim(exprs(eset))[1] else 1
	# just plot the first probe set if none supplied and more than 12 available
,	arrays = 1:dim(pData(eset))[1] # default is to use all
,	xlab = paste(colnames(pData(eset))[1:numOfFactorsToUse(eset)], collapse=":")
,	ylab = "Expression Estimate"
,	xLabels = apply(
                  as.matrix(pData(eset)[arrays,1:numOfFactorsToUse(eset)])
                , 1
                , function(mat){paste(mat, collapse=":")}
                )
,	ylim = NA
,	numOfSEs = qnorm(0.975)
,	globalYlim = FALSE # Not yet implemented!
,	plot_cols = NA
,	plot_rows = NA
,	featureNames = NA
,	showGeneNames = FALSE
# ,	showGeneNames = if(length(featureNames)==1 & is.na(featureNames[1]))
# 		FALSE else TRUE
,	showErrorBars = if(
					length(assayDataElement(eset,"se.exprs"))==0 ||
					length(assayDataElement(eset,"se.exprs")) == sum(is.na(assayDataElement(eset,"se.exprs")))
					) FALSE else TRUE
,	plotColours = FALSE
,	log.it = if(max(exprs(eset)) > 32) TRUE else FALSE
,	eset_comb = NULL
,	jitterWidth = NA
,	qtpcrData = NULL
, ...
)
{
  if(!is.na(featureNames[1]))
  {
    geneRowNumber <- function(geneName) {which(featureNames(eset)==geneName)}
    probesets <- sapply(featureNames,geneRowNumber)
  }
  cols = rep(1,length(sampleNames(eset)))
  multiplier <- 1
  for(i in 1:numOfFactorsToUse(eset))
  { levelValues <- unclass(as.factor(pData(eset)[,i])) - 1
    cols <- cols + levelValues * multiplier
    levels <- length(levels(pData(eset)[,i]))
    multiplier <- multiplier * levels
  }
  # palette(rainbow(multiplier))
  # palette(rainbow(5))
  if(is.na(plot_cols)&is.na(plot_rows))
  { plot_cols <- ceiling(sqrt(length(probesets))) }
  if(is.na(plot_cols)&!is.na(plot_rows))
  { plot_cols <- ceiling(length(probesets)/plot_rows) }
  if(is.na(plot_rows))
  { plot_rows <- ceiling(length(probesets)/plot_cols) }
  par(mfrow=c(plot_rows, plot_cols))
  for(i in probesets)
  {
    if(log.it)
	    expvals <- log2(exprs(eset)[i, arrays])
    else
	    expvals <- exprs(eset)[i, arrays]
	if(length(assayDataElement(eset,"se.exprs"))==0)
	{
		se_expvals <- expvals
		se_expvals[] <- 0
    }
	else
	{	
		se_expvals <- assayDataElement(eset,"se.exprs")[i, arrays]
    	se_expvals[which(is.na(se_expvals))] <- 0
	}
	if(numOfFactorsToUse(eset) == 1)
	{
		xValues <- unclass(pData(eset)[arrays,1])
		# xLabels <- xValues
		# xValuesForAxis <- xValues
	}
	else
		xValues <- cols
	xJitters <- rep(0,length(xValues))
	jitterCount <- 1
	for(levelCount in unique(xValues))
	{
		numInLevel <- sum(xValues==levelCount)
		if(is.null(eset_comb))
		{
			tempOffsets <- seq(
				-((numInLevel-1)/2)
			,	(numInLevel-1)/2
			)
		}
		else
		{
			if(even(numInLevel))
			{
				tempOffsets <- seq(
					-((numInLevel)/2)
				,	(numInLevel)/2
				)[-((numInLevel/2)+1)]
			}
			if(odd(numInLevel))
			{
				tempOffsets <- seq(
					-((numInLevel-1)/2)
				,	(numInLevel+1)/2
				)[-((numInLevel+1)/2)]
			}
		}
		xJitters[jitterCount:(jitterCount+(numInLevel-1))] <- tempOffsets
		jitterCount <- jitterCount + numInLevel
	}
	if(is.na(jitterWidth))
		jitterWidth <- (max(xValues) - min(xValues))/100
	xValuesJittered <- xValues + (xJitters * jitterWidth)
	# xValues <- rep(1:5,each=3)
	# xValues <- 1:length(arrays)
    plot(
      xValuesJittered
    , expvals
    , axes = FALSE
    , pch = cols[arrays]
	, type = "p"
#    , cex = 1.5
    , col = if(plotColours[1]) plotColours else cols[arrays]
    , xlab = xlab
    , ylab = ylab
    , ylim =  if(is.na(ylim[1]))
			{
				c(
					min(expvals-numOfSEs*se_expvals)
				,	max(expvals+numOfSEs*se_expvals)
				)
			}
              else
              ylim
    , ...
    )
    if(showErrorBars)
	{
		arrows(
		  xValuesJittered
		, expvals[1:length(arrays)]
		, xValuesJittered
		, expvals[1:length(arrays)] + numOfSEs*se_expvals[1:length(arrays)]
		, angle = 90
		, length = 0.05
		, col = if(plotColours[1]) plotColours else cols[arrays]
		, ...
		)
		arrows(
		  xValuesJittered
		, expvals[1:length(arrays)]
		, xValuesJittered
		, expvals[1:length(arrays)] - numOfSEs*se_expvals[1:length(arrays)]
		, angle = 90
		, length = 0.05
		, col = if(plotColours[1]) plotColours else cols[arrays]
		, ...
		)
	}
	if(!is.null(eset_comb))
	{
		expvals_comb <- exprs(eset_comb)[i,]
	    se_expvals_comb <- assayDataElement(eset_comb,"se.exprs")[i,]
	    # se_expvals_comb[which(is.na(se_expvals_comb))] <- 0
	    # xValues_comb <- unclass(pData(eset_comb)[,1])
		xValues_comb <- sort(unique(xValues))
		points(
	      xValues_comb
	      # xValues
	    , expvals_comb
	    , pch = cols[arrays]
		, type = "p"
	    # , col = if(plotColours[1]) unique(plotColours) else unique(cols[arrays])
	    , col = if(plotColours[1])
	 				unique(plotColours)
				else
					sort(unique(cols[arrays]))
	    , ...
	    )
		if(numOfFactorsToUse(eset) == 1 & is.numeric(pData(eset)[,1]))
			points(
		     	xValues_comb
		    ,	expvals_comb
		    ,	pch = cols[arrays]
			,	type = "l"
		    ,	col = "black"
			,	lwd=3
		    , ...
		    )
	    if(showErrorBars)
		{
			arrows(
			 	xValues_comb
			,	expvals_comb[1:length(arrays)]
			,	xValues_comb
			,	expvals_comb[1:length(arrays)]
					+ numOfSEs*se_expvals_comb[1:length(arrays)]
			,	angle = 90
			,	length = 0.05
			,	col = if(plotColours[1])
					unique(plotColours)
				else
					sort(unique(cols[arrays]))
			,	lwd = 3
			, ...
			)
			arrows(
			 	xValues_comb
			,	expvals_comb[1:length(arrays)]
			,	xValues_comb
			,	expvals_comb[1:length(arrays)]
					- numOfSEs*se_expvals_comb[1:length(arrays)]
			,	angle = 90
			,	length = 0.05
			,	col = if(plotColours[1])
					unique(plotColours)
				else
					sort(unique(cols[arrays]))
			,	lwd = 3
			, ...
			)
		}
	}
	if(!is.null(qtpcrData))
	{
		points(
	     	qtpcrData[,1]
	    ,	qtpcrData[,2]
	    ,	pch = 1
		,	type = "b"
	    ,	col = "black"
		,	lty=2
	    , ...
	    )
	}

    axis(
      1
    , at = xValues
#    , las = 3
    , labels = xLabels
    , ...
    )
	# axis(1, ...)
	axis(2, ...)
#    FC <- round(
#            mean(expvals[which(pheno_simple$uv[arrays]=="Yes")])
#            -
#            mean(expvals[which(pheno_simple$uv[arrays]=="No")])
#          , 3
#          )
#    title(paste("Probeset ", i, ", log(fold change)=", FC, sep = ""))
    if(showGeneNames)
    {
      title(featureNames(eset[i]))
    }
    else
    {
      title(paste("Probe set number", i, sep = ""))
    }
  }
}

