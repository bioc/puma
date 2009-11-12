/* IPPLR definition  */

#include <R.h> 
#include <Rdefines.h>
#include "global_constants.h"

typedef struct{
        /* input parameters */
	int genes;   /* number of genes/probe-set */
	int chips;         /* the number of chips */
	int conds;         /* the number of conditions */
	int *replicates;    /* indicating the replicates for each conditions */
	int max_num;              /*the number of iterations */
       
	double *data_m;   /* mean expression level */
	double *data_var;  /* variance of expression level */

	double pm[MAX_NUM_PROBE][MAX_NUM_COND];         /* PM data for one gene */
	double mm[MAX_NUM_PROBE][MAX_NUM_COND];         /* MM data for one gene */
	
	int cur_gene;    /* the current processing gene */
	double *expr;   /* expression levels for the current condition */
	double *var;  /* variance for the current condition */
        double *outp;    /* output */

    /* ----------------------parameters of variational algorithm----------------------------*/
	
    double mu0;      
	double eta0; 
	double alpha0;   
	double beta0;     
	
	double *mu1;
	double *mu1sq;
    double lamda_m;
    double loglamda;
      

    /* ------------------optimisation parameters ----------------------------*/
	
    double eps;     /* optimisation stop criteria */

} ipplrparam;


void initialparams_ipplr();

/* model optimisation */

void workout0();

void allocatemem_ipplr();
void freemem_ipplr();

/*extern void findeforc(double *x, int c, double *e_c, int *num_c);
extern double mean(double *x, int len);
extern double var(double *x, int len);
extern double fmaxn(double *x, int n);*/

void findeforc_ipplr(double *x, int c, double *e_c, int *num_c);
double mean_ipplr(double *x, int len);
double var_ipplr(double *x, int len);
double fmaxn_ipplr(double *x, int n);

SEXP hcomb_c(SEXP Mmat, SEXP Stdmat, SEXP rep,  SEXP conds,  SEXP max_num,  SEXP eps );
