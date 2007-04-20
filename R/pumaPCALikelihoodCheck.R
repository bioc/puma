pumaPCALikelihoodCheck <- function
(
    model
,   expectations
,   varY
,   Y
,   oldL
,   param
,	stepChecks
)
{
    L <- pumaPCALikelihoodBound(model, expectations, varY, Y)
    deltaL <- oldL - L
    oldL <- L
    if( stepChecks )
		cat(
	        sprintf(
	            'Likelihood change with update of %s: %2.4f\n'
	        ,   param
	        ,   deltaL
	        )
	    )
    if (deltaL < 0)
    {
        print(
            paste(
                'Likelihood drop of '
            ,   deltaL
            ,   ' after update of '
            ,   param
            ,   '.'
            )
        )
    }
    list(deltaL, L)
}
