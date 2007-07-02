"pplr" <-
function (e, control, experiment, sorted=TRUE) 
{
    dim_e <- dim(e)
    genes <- dim_e[1]
    conds <- dim_e[2]/2
    res <- matrix(0, genes, 9)
    if(length(control) > 1)
    { res[, 2] <- rowMeans(e[, control]) }
    else
    { res[, 2] <- e[, control] }
    if(length(experiment) > 1)
    { res[, 3] <- rowMeans(e[, experiment]) }
    else
    { res[, 3] <- e[, experiment] }
    if(length(control) > 1)
    { res[, 4] <- sqrt(rowSums(e[, conds + control]^2)) }
    else
    { res[, 4] <- e[, conds + control] }
    if(length(experiment) > 1)
    { res[, 5] <- sqrt(rowSums(e[, conds + experiment]^2)) }
    else
    { res[, 5] <- e[, conds + experiment] }
    res[, 6] <- res[, 3] - res[, 2]
#    res[, 7] <- sqrt(e[, conds + control]^2 + e[, conds + experiment]^2)
    res[, 7] <- sqrt((res[, 4]^2) + (res[, 5]^2))
    res[, 8] <- -res[, 6]/(res[, 7] * sqrt(2))
    res[, 9] <- erfc(res[, 8])/2
	if(sorted)
    {
		sortout <- sort(abs(res[, 8]), method = "quick", index.return = TRUE, decreasing=TRUE)
    	res[, c(2:9)] <- res[sortout$ix, c(2:9)]
    	res[, 1] <- sortout$ix
	}
    colnames(res) <- c("index", "cM", "sM", "cStd", "sStd", "LRM", 
        "LRStd", "stat", "PPLR")
	if(sorted)
	   	rownames(res) <- rownames(e)[sortout$ix]
	if(!sorted)
		rownames(res) <- rownames(e)
    return(as.data.frame(res))
}

