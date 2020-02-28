/* pumaclustii  definition  */

#include <R.h> 
#include <Rdefines.h>
#include "global_constants.h"


typedef struct{
    /* input parameters */
	int genes;   /* number of genes/probe-set */
	int chips;         /* the number of chips */
	int conds;    /* the number of conditions */
	int *reps;     /* the vector indicating the condition each chip belongs to */
	int mincls;       /* minimum number of clusters */
	int maxcls;      /* maximum number of clusters */
	double **cmu;    /* centers of clusters */ 
	double **csig;      /* sig^2 of clusters */
	
	double *data_m;   /* mean expression level */
	double *data_var;  /* variance of expression level */
    
    int verbose;
	
	//int method;   /* 0-standardm, 1-augamented */
	//double *pj;
	
	/* variables during the optimisation */
	int *cond_n;
	int i_opt;
	int var_flag;
	double *v;
	double *calpha;
	double *cbeta;
	double **q_kn;
	double **un;
	double **logun;
	double **etan;
	double **logetan;
	
	/* outputs */
	int K_best;
	double F_value;
	double **q_kn_best;
	double *pi_k_best;
	double **cmu_best;
	double **csig_best;
	

	/* optimisation parameters */
	double eps;     /* optimisation stop criteria */
	double del0;

} clustiiparam;


void initialparams_pumaclustii();

/* model optimisation */
void workout_pumaclustii();

void allocatemem_pumaclustii();
void freemem_pumaclustii();

void calMU(double **mujd, double *clsig, int nstart);
double mean(double *x, int len);
double var(double *x, int len);
double fmaxn(double *x, int n);

double *CallocD(long s);
double **CallocDD(long s);
double ***CallocDDD(long s);
int *CallocI(long s);
void FreeD(double *p);
void FreeDD(double **p);
void FreeDDD(double ***p);
void FreeI(int *p);

SEXP pumaclustii_c(SEXP Mmat, SEXP Stdmat, SEXP conds, SEXP reps, SEXP mincls, SEXP maxcls, SEXP centers, SEXP clsig, SEXP verbose, SEXP eps, SEXP del0);

