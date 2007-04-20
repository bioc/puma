pumaPCAUpdateCinv <- function
(
    model
,   expectations
,   varY
,   Y
)
{
    C <- apply(expectations@xxT, 1:2, mean) -
            ( 2 * matrix(colMeans(expectations@x)) %*% model@m ) +
            ( t(model@m) %*% model@m )
    Cinv <- solve(C)
#    C
}
