pumaPCASigmaObjective <- function
(
    sigma
,   model
,   expectations
,   varY
,   Y
)
{
    model@sigma <- sigma
    f <- pumaPCALikelihoodBound(model, expectations, varY, Y)
}
