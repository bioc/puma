pumaPCAUpdateW <- function
(
    model
,   expectations
,   varY
,   Y
)
{
    dataDim <- dim(Y)[2]
    numData <- dim(Y)[1]
    latentDim <- dim(expectations@x)[2]

    H <- matrix(0,dataDim, latentDim)
    L <- array(0, c(latentDim, latentDim, dataDim))
    W <- matrix(0,dataDim, latentDim)

##    print(sprintf("Time before ifor : %s\n", date()))
    for (i in 1:numData)
    {
        H <- H +    matrix(
                        1 / 
                        ( model@sigma ^ 2 + varY[i,] ) *
                        (Y[i,] - model@mu)
                    ) %*%
                    expectations@x[i,]
        for (j in 1:dataDim)
        {
            L[,,j] <- L[,,j] + (1/(model@sigma^2 + varY[i, j])) *
                                 expectations@xxT[,,i]
        }
    }
##    print(sprintf('Time before sum  : %s\n', date()))
    for (j in 1:dataDim)
    {
        W[j,] <- solve(L[,,j],H[j,])
    }
    W
}
