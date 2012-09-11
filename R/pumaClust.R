"pumaClust" <-
function (e = NULL, se = NULL, efile = NULL, sefile = NULL, subset = NULL, 
    gsnorm = FALSE, clusters = 10, iter.max = 100, nstart = 10, eps = 1e-06, 
    del0 = 0.01) 
{
    if(class(e)=="exprReslt" || class(e)=="ExpressionSet")
	{
		se <- assayDataElement(e,"se.exprs")
		e <- exprs(e)
	}

	if (length(e) == 0) {
        e <- read.csv(efile, check.names = FALSE, row.names = 1)
    }
    dim_e <- dim(e)
    genes <- dim_e[1]
    chips <- dim_e[2]
    if (length(se) == 0 & length(sefile) > 0) {
        se <- read.csv(sefile, check.names = FALSE, row.names = 1)
        se <- se^2
        m <- 1
    }
    else {
        if (length(se) == 0 & length(sefile) == 0) {
            se <- matrix(0, genes, chips)
            m <- 0
        }
        else {
            se <- se^2
            m <- 1
        }
    }
    if (gsnorm == TRUE) 
        e <- normalisation.gs(e)
    if (length(subset) > 0) {
        e <- e[subset, ]
        se <- se[subset, ]
        dim_e <- dim(e)
        genes <- dim_e[1]
        chips <- dim_e[2]
    }
    if (m == 1) {
        se <- t(apply(cbind(e, se), 1, clusterNormVar))
    }
    e <- t(apply(e, 1, clusterNormE))
 
    cl <- kmeans(e, clusters, iter.max = iter.max, nstart = nstart)
   
    clsig <- 0 * c(1:clusters)
    for (i in 1:clusters) {
        clsig[i] <- mean(diag(cov(e[which(cl$cluster == i), ])))
    }
    res <- .Call("pumaclust_c", as.matrix(e), as.matrix(se), 
        clusters, cl$centers, clsig, eps, del0, PACKAGE = "puma")
    names(res[[1]]) <- rownames(e)
    colnames(res[[2]]) <- colnames(e)
    rownames(res[[2]]) <- paste(1:clusters)
    names(res[[3]]) <- paste(1:clusters)
    rownames(res[[4]]) <- rownames(e)
    colnames(res[[4]]) <- paste(1:clusters)
    return(list(cluster = res[[1]], centers = res[[2]], centersigs = res[[3]], 
        likelipergene = res[[4]], bic = res[[5]]))
}

clusterNormE <- function(x)
{
    meanExp <- mean(x)
    stdExp <- sqrt(var(x))
    x <- (x-meanExp)/stdExp
    return(x)
}

clusterNormVar <- function(x)
{
    n <- length(x)
    stdExp <- var(x[c(1:(n/2))])
    x <- x[c((n/2+1):n)]/stdExp
    return(x)
}

