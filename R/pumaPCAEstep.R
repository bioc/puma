pumaPCAEstep <- function
(
    model
,   expectations
,   varY
,   Y
)
{
    dataDim <- dim(Y)[2]
##  latentDim <- size(expectations.x, 2);
    numData <- dim(Y)[1]

    for (i in 1:numData)   #computes the inverses of the matrices M_n
    {
        Binv <- diag(1/(model@sigma^2 + varY[i,]))
  ## This inverse should really be done across a q dimensional matrix.
        SigmaInv <- t(model@W) %*% Binv %*% model@W + model@Cinv
#        print(SigmaInv)
#        Sigma_x <- solve(SigmaInv)
#        print(Sigma_x)
        expectations@logDetCov[i] <- -log(det(SigmaInv))
#        print(model@W)
#        print(Binv)
#        print(model@mu)
#        expectations@x[i,] <- (Sigma_x %*% (t(model@W) %*% Binv %*% (Y[i,] - model@mu) + model@Cinv %*% model@m))
#        expectations@xxT[,,i] <- Sigma_x + expectations@x[i,] %*% t(expectations@x[i,])
        expectations@x[i,] <- (solve(SigmaInv, (t(model@W) %*% Binv %*% t(Y[i,] - model@mu) + model@Cinv %*% t(model@m))))
        expectations@xxT[,,i] <- solve(SigmaInv) + expectations@x[i,] %*% t(expectations@x[i,])
    ##  expectations.xTx(i) = trace(expectations.xxT(:,:,i));
    }
    expectations
}
