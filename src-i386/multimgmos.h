/* multi-mgMOS definition  */

#include <R.h> 
#include <Rdefines.h>
#include "global_constants.h"

typedef struct{
        /* input parameters */
	int conds;    /* number of conditions */
	int genes;   /* number of genes/probe-set */
	int chips;         /* the number of chips */
	long probes;            /* the number of probes */
	double phi;
	double mu;     /* mean of fitted lognormal distribution of phi */
	double sigma;  /* variance of fitted lognormal distribution of phi */
	
	
	double *data_pm;   /* PM data for the whole data set */
	double *data_mm;   /* MM data for the whole data set */
	double pm[MAX_NUM_PROBE][MAX_NUM_COND];         /* PM data for one gene */
	double mm[MAX_NUM_PROBE][MAX_NUM_COND];         /* MM data for one gene */
	int *probesets;   /* the number of probe-pairs for each probe set */
	int *replicates;   /* replicates for each condition */
	long totalprobe;    /* the total processed probe */
	int cur_gene;      /* current processed gene */
	int cur_cond;      /* current processed condition -- for mgMOS only */
	int num_probe;       /* the number of probe-pairs in the currently processed gene */
	
	double *paramphi;   /* iteration number, phi, rc, the optimised likelihood function with respect to phi */
			    
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
  
} expparam;

double dierfc(double y);  /* inverse error function */

double erfc(double x);    /* complementary error function */

void mmgmos_initialparams();

void getgenedata(int g);
	
/* model optimisation */
void calparameters();
	
/* calculate percentiles of expression levels */
void calexpression();
	
void mmgmos_allocatemem();
void mmgmos_freemem();


void workout_mgmos();

SEXP mmgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps);

SEXP mgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps);

