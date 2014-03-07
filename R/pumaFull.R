pumaFull <- function(
	ExpressionFeatureSet = NULL
,	data_dir = getwd()
,	load_ExpressionFeatureSet = FALSE
,	calculate_eset = TRUE
,	calculate_pumaPCAs = TRUE
,	calculate_bcomb = TRUE
,	mmgmosComparisons = FALSE
)
{
	if(is.null(ExpressionFeatureSet) && !load_ExpressionFeatureSet && calculate_eset)
	{
                    setwd(data_dir);
	 	ExpressionFeatureSet<-read.celfiles(list.celfiles());
		save(ExpressionFeatureSet, file=paste(data_dir,"ExpressionFeatureSet.rda",sep="/"))
	}

	if(load_ExpressionFeatureSet && calculate_eset)
	{
	  load(paste(data_dir,"ExpressionFeatureSet.rda",sep="/"))
	}

	if(calculate_eset)
	{
		cat("calculating eset_mmgmos\n")
		eset_mmgmos <- mmgmos(ExpressionFeatureSet)
		save(eset_mmgmos, file=paste(data_dir,"/eset_mmgmos.rda",sep=""))
		if(mmgmosComparisons)
		{
			cat("calculating eset_mmgmos_bg\n")
			eset_mmgmos_bg <- mmgmos(
				ExpressionFeatureSet
			,	background=TRUE
			,	replaceZeroIntensities=FALSE
			)
			save(
				eset_mmgmos_bg
			,	file=paste(data_dir,"/eset_mmgmos_bg.rda",sep="")
			)
		}
		cat("calculating eset_rma\n")
		eset_rma <- oligo:::rma(ExpressionFeatureSet)
		save(eset_rma, file=paste(data_dir,"/eset_rma.rda",sep=""))
	}

	if(!calculate_eset)
	{
		load(paste(data_dir,"/eset_mmgmos.rda",sep=""))
		load(paste(data_dir,"/eset_rma.rda",sep=""))
		if(mmgmosComparisons)
		{
			load(paste(data_dir,"/eset_mmgmos_bg.rda",sep=""))
		}
	}

	##############################
	## pumaPCA
	##############################
	if(calculate_pumaPCAs)
	{
		pumaPCA_results <- pumaPCA(eset_mmgmos)
		save(pumaPCA_results, file=paste(data_dir,"/pumaPCA_results.rda",sep=""))
	}
	if(!calculate_pumaPCAs)
	{
		load(paste(data_dir,"/pumaPCA_results.rda",sep=""))
	}

	pdf("pca.pdf",width=7, height=10)
	par(mfrow=c(2,1))
	plot(pumaPCA_results,legend1pos="topleft",legend2pos="top", main="pumaPCA")
	pca_results <- prcomp(t(exprs(eset_rma)))
	if(dim(pData(eset_rma))[2] > 1)
	{	
		plot(
			pca_results$x
		,	xlab="Component 1"
		,	ylab="Component 2"
		,	pch=unclass(as.factor(pData(eset_rma)[,1]))
		,	col=unclass(as.factor(pData(eset_rma)[,2]))
		,	main="PCA"
		)
	}
	else
	{
		plot(
			pca_results$x
		,	xlab="Component 1"
		,	ylab="Component 2"
		,	pch=unclass(as.factor(pData(eset_rma)[,1]))
		,	col=unclass(as.factor(pData(eset_rma)[,1]))
		,	main="PCA"
		)
	}
	dev.off()

	##############################
	## Genes with highest interaction
	##############################

	pdf("boxplots.pdf",width=10, height=7)
	par(mfcol=c(2,2))
	boxplot(data.frame(exprs(eset_mmgmos)), main="mmgMOS - No norm")
	boxplot(data.frame(exprs(eset_rma)), main="Standard RMA")
	eset_mmgmos_normd <- pumaNormalize(eset_mmgmos, "median")
	boxplot(data.frame(exprs(eset_mmgmos_normd)), main="mmgMOS - median scaling")
	dev.off()

	if(calculate_bcomb)
	{
		eset_comb <- pumaComb(eset_mmgmos)
		save(eset_comb, file=paste(data_dir,"/eset_comb.rda",sep=""))
		eset_normd_comb <- pumaComb(eset_mmgmos_normd)
		save(
			eset_normd_comb
		,	file=paste(data_dir,"/eset_normd_comb.rda",sep="")
		)
	}
	if(!calculate_bcomb)
	{
		load(paste(data_dir,"/eset_comb.rda",sep=""))
		load(paste(data_dir,"/eset_normd_comb.rda",sep=""))
	}

	compareLimmapumaDE(
		eset_mmgmos
	,	eset_comb
	,	eset_rma
	,	contrastsFilename="contrasts_mmgmos"
	,	plotOther=TRUE
	,	otherFilename="contrasts_rma"
	)
	compareLimmapumaDE(
		eset_mmgmos_normd
	,	eset_normd_comb
	,	eset_rma
	,	contrastsFilename="contrasts_mmgmos_normd"
	)
}
