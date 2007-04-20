pumaPCANewtonUpdateLogSigma <- function
(
    model
,   expectations
,   varY
,   Y
)
{
    sigma <- model@sigma
    deltaL1 <- 1
    maxIters <- 10
    counter2 <- 0
    Step <- 1
    while (deltaL1 > 1e-5 & counter2 <= 10 & Step > 1e-5)
    {
        deltaL <- -1
        counter2 <- counter2 + 1
        oldL <- pumaPCASigmaObjective(sigma, model, expectations, varY, Y)
        Step <- newtonStep(sigma, model, expectations, varY, Y)
        counter <- 0
        while (deltaL < 0)
        {
            counter <- counter + 1
            L <- pumaPCASigmaObjective(sigma-Step, model, expectations, varY, Y)
            deltaL <- oldL - L
            Step <- Step/2
        }
        deltaL1 <- deltaL
        sigma <- sigma - 2*Step
    }
#   Step;
    as.vector(sigma)
}

newtonStep <- function
(
    sigma
,   model
,   expectations
,   varY
,   Y
)
{
    eps <- 0 ## ???
    if (sigma < eps)
    {
        sigma <- eps
    }
    eta <- 0.1

    dataDim <- dim(Y)[2]
    numData <- dim(Y)[1]

    Sigma <- 0
    curvature <- 0

    for (i in 1:numData)
    {
        Fun <- diag(1/(sigma^2+varY[i,]))
        F2 <- Fun*Fun
        F3 <- F2*Fun
        ## computes the gradient with respect to sigma
        Sigma <- Sigma -
            2 * ( (Y[i,] - model@mu) %*% F2 %*% t(Y[i,] - model@mu) ) * sigma +
            2 * sum(diag(Fun)) * sigma +
            ( 4 * expectations@x[i,] %*% t(model@W) %*% F2 %*% t(Y[i,] - model@mu) ) * sigma -
            2 * sum(diag(t(model@W) %*% F2 %*% model@W %*% expectations@xxT[,,i])) * sigma
        ## computes the curvature with respect to sigma
        curvature <- curvature +
            8 * ( (Y[i,] - model@mu) %*% F3 %*% t(Y[i,] - model@mu) ) * sigma -
            4 * sum(diag(F2)) * sigma +
            ( -16 * expectations@x[i,] %*% t(model@W) %*% F3 %*% t(Y[i,] - model@mu) ) * sigma +
            8 * sum(diag(t(model@W) %*% F3 %*% model@W %*% expectations@xxT[,,i])) * sigma
    }
    l <- log(sigma)
    Grad <- 0.5 * Sigma
    Curv <- ( 0.5 * curvature * sigma + Grad/sigma) * sigma * sigma / 4
    Step.l <- Grad*exp(l)/(Curv*exp(2*l)+Grad*exp(l))

    Step <- exp(l-Step.l)
}
