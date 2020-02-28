/* pumaclust definition  */

#include <R.h> 
#include <Rdefines.h>
#include "global_constants.h"


typedef struct{
        /* input parameters */
	int genes;   /* number of genes/probe-set */
	int chips;         /* the number of chips */
	int clusters;      /* number of clusters */
	double *centers;    /* centers of clusters */ 
	double *clsig;      /* sig^2 of clusters */
	
	double *data_m;   /* mean expression level */
	double *data_var;  /* variance of expression level */
	
	int method;   /* 0-standardm, 1-augamented */
	double **pjPerData;
	double *pj;
	
	int *CIndex;    /* output */
	double *ClusterCenters;
	double *ClusterSig;
	double *LikeliPerGene;
	double *BIC;

	/* optimisation parameters */
	double eps;     /* optimisation stop criteria */
	double del0;

} clusterparam;


void initialparams_pumaclust();

/* model optimisation */
void workout();

void allocatemem_pumaclust();
void freemem_pumaclust();

void calMU(double **mujd, double *clsig, int nstart);
 double mean(double *x, int len);
 double var(double *x, int len);
 double fmaxn(double *x, int n);

SEXP pumaclust_c(SEXP Mmat, SEXP Stdmat, SEXP clusters, SEXP centers, SEXP clsig, SEXP eps, SEXP del0);
