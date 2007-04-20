pumaPCASigmaGradient <- function
(
    sigma
,   model
,   expectations
,   varY
,   Y
)
{
    dataDim <- dim(Y)[2]
    numData <- dim(Y)[1]

    Sigma <- 0
    curvature <- 0

    for (i in 1:numData)
    {
        Fun <- diag(1/(sigma^2+varY[i,]))
        F2 <- Fun*Fun
        ## computes the gradient with respect to sigma
        Sigma <- Sigma -
            2 * ( (Y[i,] - model@mu) %*% F2 %*% t(Y[i,] - model@mu) ) * sigma +
            2 * sum(diag(Fun)) * sigma +
            ( 4 * expectations@x[i,] %*% t(model@W) %*% F2 %*% t(Y[i,] - model@mu) ) * sigma -
            2 * sum(diag(t(model@W) %*% F2 %*% model@W %*% expectations@xxT[,,i])) * sigma
    }
    Grad <- 0.5 * Sigma
    as.vector(Grad)
}
