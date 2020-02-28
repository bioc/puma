/* PMmulti-mgMOS definition  */

#include <R.h> 
#include <Rdefines.h>
#include "global_constants.h"

typedef struct{
        /* input parameters */
	int conds;    /* number of conditions */
	int genes;   /* number of genes/probe-set */
	int chips;         /* the number of chips */
	long probes;            /* the number of probes */
	
	
	double *data_pm;   /* PM data for the whole data set */
	
	double pm[MAX_NUM_PROBE_PM][MAX_NUM_COND];         /* PM data for one gene */
	
	int *probesets;   /* the number of probe-pairs for each probe set */
	
	long totalprobe;    /* the total processed probe */
	int cur_gene;      /* current processed gene */
	int cur_cond;      /* current processed condition -- for mgMOS only */
	int num_probe;       /* the number of probe-pairs in the currently processed gene */
           double *prctiles;     /* the percentiles of expression */
	int num_prctile;     /* number of percentiles */
	
	double *outp;    /* intervals of expression,  default 50% */
	double **parameters;    /* estimated parameters */
	                               /* conds of alpha, chips of a, c, d, rc, fopt for each gene */
	double par_mgmos[5];    /* parameters of mgMOS */
	
	/* optimisation parameters */
	
	double eps;     /* optimisation stop criteria */
	
	/* expression levels calculation parameters */
	double step0;
	
	int flag;   /* 0: gene optimisation; 1: phi optimisation; 2: mgmos optimisation */
	int saveflag;   /* FALSE: not save parameters; TRUE: save parameters */
  
} pmexpparam;

double pmdierfc(double y);  /* inverse error function */

double pmerfc(double x);    /* complementary error function */

void mmgmos_initialparamspm();

void pmgetgenedata(int g);
	
/* model optimisation */
void pmcalparameters();
	
/* calculate percentiles of expression levels */
void pmcalexpression();
	
void mmgmos_allocatemempm();
void mmgmos_freemempm();



SEXP pmmmgmos_c(SEXP PMmat, SEXP ngenes, SEXP probeNames,SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps);




