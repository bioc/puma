/* PPLR definition  */

#include <R.h> 
#include <Rdefines.h>
#include "global_constants.h"


typedef struct{
        /* input parameters */
	int genes;   /* number of genes/probe-set */
	int chips;         /* the number of chips */
	int conds;         /* the number of conditions */
	int *replicates;    /* indicating the replicates for each conditions */
	
	double *data_m;   /* mean expression level */
	double *data_var;  /* variance of expression level */
	double pm[MAX_NUM_PROBE][MAX_NUM_COND];         /* PM data for one gene */
	double mm[MAX_NUM_PROBE][MAX_NUM_COND];         /* MM data for one gene */
	
	int cur_gene;    /* the current processing gene */
	double *expr;   /* expression levels for the current condition */
	double *var;  /* variance for the current condition */
	
	int method;   /* Bayesion method, 0-hier_shalam, 1-hier_conj_shalam, 2-hier_conj, 3-MCMC */
	
	
	double *outp;    /* output */

	/* parameters of variational algorithm */
	int sample_num;	
	double alpha0;
	double beta0;
	double mu0;
	double eta0;
	double *mu1sq;
	double *mu1;
	double lamda_m;
	double loglamda;
	
	/* optimisation parameters */
	double eps;     /* optimisation stop criteria */

} pplrparam;


void initialparams_pplr();

/* model optimisation */
void workout_pplr0();
void workout_pplr1();

void allocatemem_pplr();
void freemem_pplr();
void findeforc(double *x, int c, double *e_c, int *num_c);

double mean(double *x, int len);
double var(double *x, int len);
double fmaxn(double *x, int n);

SEXP bcomb_c(SEXP Mmat, SEXP Stdmat, SEXP rep, SEXP method, SEXP conds, SEXP nsample, SEXP eps);
