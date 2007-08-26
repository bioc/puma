numFP <- function (
	scores
,	truthValues
,	TPRate = 0.5
)
{
	numberOfTPToCount <- length(which(truthValues)) * TPRate
	scoresIndex <- sort(scores,decreasing=TRUE,index.return=TRUE)$ix
	numberOfFP <- 0
	numberOfTP <- 0
	for(i in 1:length(scoresIndex))
	{
		if(truthValues[scoresIndex[i]]) numberOfTP <- numberOfTP+1 else numberOfFP <- numberOfFP+1
		if(numberOfTP >= numberOfTPToCount) return(numberOfFP)
	}
}

numTP <- function (
	scores
,	truthValues
,	FPRate = 0.5
)
{
	numberOfFPToCount <- length(which(!truthValues)) * FPRate
	scoresIndex <- sort(scores,decreasing=TRUE,index.return=TRUE)$ix
	numberOfFP <- 0
	numberOfTP <- 0
	for(i in 1:length(scoresIndex))
	{
		if(truthValues[scoresIndex[i]]) numberOfTP <- numberOfTP+1 else numberOfFP <- numberOfFP+1
		if(numberOfFP > numberOfFPToCount) return(numberOfTP)
	}
}

