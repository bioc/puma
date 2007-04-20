pumaPCAUpdateMu <- function
(
    model
,   expectations
,   varY
,   Y
)
{
    dataDim <- dim(Y)[2]
    numData <- dim(Y)[1]

    G <- matrix(nrow=dataDim, ncol=numData)

    for (i in 1:numData)
    {
        G[, i] <- ( 1 / (model@sigma^2 + varY[i,] ) ) * ( Y[i,] - model@W %*% expectations@x[i, ] )
    }
    Q <- colSums( 1 / ( model@sigma^2 + varY ) )
    mu <- t(matrix(rowSums(G) / Q))
    mu
}

