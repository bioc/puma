pumaPCARemoveRedundancy <- function
(
    model
)
{
    latentDim <- dim(model@W)[2]
    model@mu <- model@mu + model@m %*% t(model@W)
    Cmatrix <- model@W %*% solve(model@Cinv) %*% t(model@W)
    Cmatrix <- Cmatrix/2 + t(Cmatrix)/2
    eigenC <- eigen(Cmatrix)
    U <- eigenC$vectors
    V <- abs(eigenC$values)
##  V = diag(V)
##  [V, index] = sort(V);
##  index = index(end:-1:1);
##  V = V(end:-1:1);
##  U = U(:, index);

    oldW <- model@W
    model@W <- U[, 1:latentDim] %*% diag(sqrt(V[1:latentDim]))
    for(i in 1:dim(model@W)[2])
    {
        if(dist(rbind(model@W[,i],oldW[,i]))
            >
            dist(rbind(model@W[,i],-oldW[,i]))
        )
        {
            model@W[,i] <- -model@W[,i]
        }
    }
    model@Cinv <- diag(dim(model@Cinv)[1])
    model@m <- matrix(0,1,dim(model@m)[2])

    model
}
