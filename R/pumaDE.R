pumaDE <- function (
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
)
{
  p <- matrix(0, dim(exprs(eset))[1], dim(contrast.matrix)[2])
  genes <- matrix(0, dim(exprs(eset))[1], dim(contrast.matrix)[2])
  for(i in 1:dim(contrast.matrix)[2])
  {
    p[,i] <- pplr(
                cbind(exprs(eset), assayDataElement(eset,"se.exprs"))
              , which(contrast.matrix[,i]==-1)
              , which(contrast.matrix[,i]==1)
            )[,9]
    genes[,i] <-  pplr(
                    cbind(exprs(eset), assayDataElement(eset,"se.exprs"))
                  , which(contrast.matrix[,i]==-1)
                  , which(contrast.matrix[,i]==1)
                  )[,1]
  }
  list(p=p, genes=genes)
}

pumaDEUnsorted <- function (
	pp
)
{
	p <- matrix(0,dim(pp[[1]])[1], dim(pp[[1]])[2])
	for(i in 1:dim(pp[[1]])[2])
	{
		sortpumaDE <- sort(pp[[2]][,i], method="quick", index.return=TRUE)[[2]]
		p[,i] <- pp[[1]][sortpumaDE,i]
	}
	return(p)
}

pplrUnsorted <- function (
	p
)
{
	sortPplr <- sort(p[,1], method="quick", index.return=TRUE)[[2]]
	pSorted <- p[,9][sortPplr]
	names(pSorted) <- rownames(p[sortPplr,])
	return(pSorted)
}
