pumaPCALikelihoodBound <- function
(
    model
,   expectations
,   varY
,   Y
)
{
    dataDim <- dim(Y)[2]
    numData <- dim(Y)[1]

#    z <- rep(0,numData)
#    s <- rep(0,numData)
#    s1 <- rep(0,numData)
#    s2 <- rep(0,numData)
#    s3 <- matrix(0, numData, dataDim)
#    trxHatCxHatT <- rep(0, numData)
    sigma2 <- model@sigma * model@sigma

#    for (i in 1:numData)
#    {
#        yHat <- Y[i,] - model@mu
#        Binv <- 1/t(sigma2 + varY[i,])
#        z[i] <- -sum(log(Binv))
#        s[i] <- yHat %*% t(Binv * yHat)
#        s1[i] <- expectations@x[i,] %*% t(model@W) %*% t(Binv*yHat)
#		s2[i] <- sum(Binv[] * diag(model@W[,] %*% expectations@xxT[,,i] %*% t(model@W[,])))
#        for (j in 1:dataDim)
#        {
#            s2[i] <- s2[i] + Binv[j] * model@W[j,] %*% expectations@xxT[,,i] %*% model@W[j,]
#        }

#        trxHatCxHatT[i] <- sum(expectations@xxT[,,i] * t(model@Cinv)) - 2 * expectations@x[i,] %*% model@Cinv %*% t(model@m) + model@m %*% model@Cinv %*% t(model@m)
#    }

        yHat <- t(Y) - as.vector(model@mu)
        Binv <- 1/t(sigma2 + varY)
        z <- -sum(log(1/t(sigma2 + varY[,])))
        s <- sum(yHat * ( Binv * yHat))
##        s <- sum(diag(t(yHat) %*% ( Binv * yHat)))
##        s <- sum(diag(t(t(Y[,]) - as.vector(model@mu)) %*% ( (1/t(sigma2 + varY[,])) * (t(Y[,]) - as.vector(model@mu)))))
        s1 <- sum( t( expectations@x[,] %*% t(model@W) ) * (Binv*yHat) )
#        s1 <- sum(diag(expectations@x[,] %*% t(model@W) %*% (Binv*yHat) ) )
        s2 <- 0
        trxHatCxHatT <- as.vector(0)
	    for (i in 1:numData)
    	{
			s2 <- s2 + sum(Binv[,i] * diag(model@W[,] %*% expectations@xxT[,,i] %*% t(model@W[,])))
	        trxHatCxHatT <- trxHatCxHatT + as.vector(sum(expectations@xxT[,,i] * t(model@Cinv)) - 2 * expectations@x[i,] %*% model@Cinv %*% t(model@m) + model@m %*% model@Cinv %*% t(model@m))
		}


#    s2 <- rowSums(s3)
    f <- 0.5 * (z + trxHatCxHatT + s - 2*s1 + s2 - sum(expectations@logDetCov)) - 0.5 * numData * log(det(model@Cinv))
#    f <- f - 0.5 * sum(expectations@logDetCov)
}
