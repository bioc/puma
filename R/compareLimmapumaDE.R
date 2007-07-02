compareLimmapumaDE <- function (
	eset_mmgmos
,	eset_comb = NULL
,	eset_other = eset_mmgmos
,	limmaRes = calculateLimma(eset_other)
,	pumaDERes = pumaDE(eset_comb)
,	contrastMatrix = createContrastMatrix(eset_mmgmos)
,	numberToCompareForContrasts = 3
,	numberToCompareForVenn = 100
,	plotContrasts = TRUE
,	contrastsFilename = NULL
,	plotOther = FALSE
,	otherFilename = "other"
,	plotBcombContrasts = FALSE
,	bcombContrastsFilename = "bcomb_contrasts"
,	plotVenn = FALSE
,	vennFilename = "venn.pdf"
,	showTopMatches = FALSE
,	returnResults = FALSE
)
{
	numberOfContrasts <- dim(contrastMatrix)[2]
	numberOfGenes <- numberOfGenes(limmaRes)
	if(numberOfGenes != numberOfGenes(pumaDERes))
	{ stop("Different numbers of genes") }
	for(i in 1:numberOfContrasts)
	{
		genes_to_plot <- as.integer(
			c(
				topGenes(pumaDERes,numberOfGenes=numberToCompareForContrasts,contrast=i)
			,	topGenes(limmaRes,numberOfGenes=numberToCompareForContrasts,contrast=i)
			)
		)
		if(plotContrasts)
		{
			if(!is.null(contrastsFilename))
			{
				pdf(
					paste(
						contrastsFilename
					,	"_"
					,	colnames(contrastMatrix)[i]
					,	".pdf"
					,	sep=""
					)
				,	width=10
				,	height=6
				,	pointsize=6
				)
			}
			# else
			# 	do.call(getOption("device"),list())
			if(is.null(eset_comb))
				plotErrorBars(eset_mmgmos, genes_to_plot, showGeneNames=TRUE)
			else
				plotErrorBars(eset_mmgmos, genes_to_plot, showGeneNames=TRUE, eset_comb=eset_comb)
			if(!is.null(contrastsFilename)) dev.off()
		}
		if(plotOther)
		{
			if(!is.null(otherFilename))
			{
				pdf(
					paste(
						otherFilename
					,	"_"
					,	colnames(contrastMatrix)[i]
					,	".pdf"
					,	sep=""
					)
				,	width=10
				,	height=6
				,	pointsize=6
				)
			}
			# else
			# 	do.call(getOption("device"),list())
			plotErrorBars(eset_other, genes_to_plot, showGeneNames=TRUE)
			if(!is.null(contrastsFilename)) dev.off()
		}
		if(plotBcombContrasts && !(is.null(eset_comb)))
		{
			if(!is.null(bcombContrastsFilename))
			{
				pdf(
					paste(
						bcombContrastsFilename
					,	"_"
					,	colnames(contrastMatrix)[i]
					,	".pdf"
					,	sep=""
					)
				,	width=10
				,	height=6
				,	pointsize=6
				)
			}
			else
				do.call(getOption("device"),list())
			plotErrorBars(eset_comb, genes_to_plot, showGeneNames=TRUE)
			if(!is.null(bcombContrastsFilename)) dev.off()
		}
		if(showTopMatches)
		{
			cat(paste("Contrast",colnames(contrastMatrix)[i]),"\n")
			print("Position of limma top matches in pumaDE list")
			print(
				match(
					topGenes(limmaRes,numberOfGenes=numberToCompareForContrasts,contrast=i)
				,	topGenes(pumaDERes,numberOfGenes=numberOfGenes,contrast=i)
				)
			)
			print("PPLR of limma top matches")
			print(
				statistic(pumaDERes)[match(
					topGenes(limmaRes,numberOfGenes=numberToCompareForContrasts,contrast=i)
				,	topGenes(pumaDERes,numberOfGenes=numberOfGenes,contrast=i)
				)
				,i]
			)
			print("")
		}
	}
	results <- list()
	for(i in 1:numberOfContrasts)
	{
		topLimma <- rep(0, numberOfGenes)
		topLimma[topGenes(limmaRes,numberOfGenes=numberToCompareForVenn,contrast=i)] <- 1
		toppumaDE <- rep(0, numberOfGenes)
		toppumaDE[topGenes(pumaDERes,numberOfGenes=numberToCompareForVenn,contrast=i)] <- 1
		results[[i]] <- cbind (topLimma, toppumaDE)
	}
	if(plotVenn)
	{
		plotColumns <- ceiling(sqrt(numberOfContrasts))
		plotRows <- ceiling(numberOfContrasts/plotColumns)
		pdf(vennFilename, paper="a4", pointsize=6)
		par(mfrow=c(plotRows, plotColumns))
		for(i in 1:numberOfContrasts)
		{
			vennDiagram(results[[i]], main=colnames(contrastMatrix)[i])
		}
		dev.off()
	}
	if(returnResults) return(results)
}

