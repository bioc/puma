#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "multimgmos.h"
#include "pplr_c.h"
#include "pumaclust_c.h"

/*SEXP mmgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps);

SEXP mgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps);

SEXP bcomb_c(SEXP Mmat, SEXP Stdmat, SEXP rep, SEXP method, SEXP conds, SEXP nsample, SEXP eps);

SEXP pumaclust_c(SEXP Mmat, SEXP Stdmat, SEXP clusters, SEXP centers, SEXP clsig, SEXP eps, SEXP del0);

*/

static const R_CallMethodDef callMethods[] = {
	{"mmgmos_c", (DL_FUNC) &mmgmos_c, 9},
	{"mgmos_c", (DL_FUNC) &mgmos_c, 9},
	{"bcomb_c", (DL_FUNC) &bcomb_c, 7},
	{"pumaclust_c", (DL_FUNC) &pumaclust_c, 7},
	{NULL, NULL, 0}
};

void R_init_puma(DllInfo *info)
{
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}
