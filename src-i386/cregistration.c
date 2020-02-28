#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "multimgmos.h"
//#include "multimgmosii.h"
#include "pplr_c.h"
#include "pumaclust_c.h"
#include "ipplr_c.h"
#include "pumaclustii_c.h"
#include "PMmultimgmoshead.h"
#include  "gmehead.h"

/*SEXP mmgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps);

SEXP mgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps);

SEXP bcomb_c(SEXP Mmat, SEXP Stdmat, SEXP rep, SEXP method, SEXP conds, SEXP nsample, SEXP eps);

SEXP pumaclust_c(SEXP Mmat, SEXP Stdmat, SEXP clusters, SEXP centers, SEXP clsig, SEXP eps, SEXP del0);

SEXP hcomb_c(SEXP Mmat, SEXP Stdmat, SEXP rep,  SEXP conds,  SEXP max_num,  SEXP eps );

SEXP pumaclustii_c(SEXP Mmat, SEXP Stdmat, SEXP conds, SEXP reps, SEXP mincls, SEXP maxcls, SEXP centers, SEXP clsig, SEXP verbose, SEXP eps, SEXP del0);

SEXP gme_c(SEXP PMmat, SEXP GTmat, SEXP PNmat, SEXP GNmat, SEXP ANmat,SEXP totalgene, SEXP saveflag, SEXP eps, SEXP pm_index, SEXP gt_index);

SEXP pmmmgmos_c(SEXP PMmat, SEXP ngenes, SEXP probeNames,SEXP prctiles, SEXP nprc,  SEXP saveflag, SEXP eps);
*/

static const R_CallMethodDef callMethods[] = {
	{"mmgmos_c", (DL_FUNC) &mmgmos_c, 9},
 //         {"mmgmosii_c", (DL_FUNC) &mmgmosii_c, 9},
	{"mgmos_c", (DL_FUNC) &mgmos_c, 9},
	{"bcomb_c", (DL_FUNC) &bcomb_c, 7},
	{"pumaclust_c", (DL_FUNC) &pumaclust_c, 7},
	{"hcomb_c",(DL_FUNC)&hcomb_c,6},
	{"pumaclustii_c",(DL_FUNC)&pumaclustii_c,11},
          {"pmmmgmos_c",(DL_FUNC)&pmmmgmos_c,7},
          {"gme_c", (DL_FUNC) &gme_c,10},
	{NULL, NULL, 0}
};

void R_init_puma(DllInfo *info)
{
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  gme_expparam_init();
}

void R_unload_puma(DllInfo *info)
{
    gme_expparam_free();
}
