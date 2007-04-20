matrixDistance <- function(
    matrixA
,   matrixB
)
{
     distance <- 0
     for(i in 1:dim(matrixA)[1])
     {  distance <- distance +
                    dist(
                        rbind(
                            matrixA[i,]
                        ,   matrixB[i,]
                        )
                    )[1]
    }
    distance <- distance / dim(matrixA)[1]
}
