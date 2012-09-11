/* GME definition  */

#include <R.h> 
#include <Rdefines.h>
#include "global_constants.h"

typedef struct{
        /* input parameters */
	int conds;         /* number of conditions */
	int genes;         /* number of total genes */
	int chips;         /* the number of chips */
        long probes; 
         int gtdim;/////
       /* number of total probes */
	//int num_probe;      /* the number of probes in current processed gene */
	//int num_gene;       /* num of total probes in current processed gene */
	//int num_alpha;      /* num of alpha in current processed gene */

        int cur_num_probe;   /* the number of probes in current processed gene */
        int cur_num_gene;   /* num of total probes in current processed gene */
        int cur_num_alpha;   /* num of alpha in current processed gene */
        int cur_gene;      
                          
        int numofgenes;
	
	double *data_GT;   /* corresponding bewteen gene and transcript for the whole data set */
	int *data_GN;   /* num of gene for the whole data set*/
	int *data_PN;   /* num of probe for the whole data set*/
        int  *data_AN;    /* num of alpha for everygene */
        double *data_pm;   /* PM data for the whole data set */
	double *data_gt;    /* corresponding bewteen everygene and transcript */


       
	double pm[MAX_NUM_PROBE_GME][MAX_NUM_COND_GME];         /* PM data for one gene */
	double gt[10000][3];         /* corresponding for one gene */
        int num_gene[MAX_NUM_GENE_GME];
        int num_probe[MAX_NUM_GENE_GME];
        int num_alpha[MAX_NUM_GENE_GME];
	double MB[MAX_NUM_PROBE_GME][MAX_NUM_COND_GME];
	double G[MAX_NUM_PROBE_GME][MAX_NUM_COND_GME]; 
    double G1[MAX_NUM_COND_GME]; 
	double alpha[MAX_NUM_COND_GME];
	double alphaii[MAX_NUM_COND_GME];
	double alpha2[MAX_NUM_COND_GME];
	double  alphaii2[MAX_NUM_COND_GME];
    double Dif[MAX_NUM_PROBE_GME][MAX_NUM_COND_GME];
    



	long totalprobe;    /* the total processed probe */
	long totalgene;    /* the total processed gene */
		
	double **parameters;    /* estimated parameters */
	                               /* conds of alpha, c, d, fopt for each gene */
	//double **cd;
	/* optimisation parameters */
	
	double eps;     /* optimisation stop criteria */
	              	/* expression levels calculation parameters */
	double step0;
	int saveflag;   /* FALSE: not save parameters; TRUE: save parameters */
        double *outp; 
  
} gme_expparam;

void gme_initialparams();

void gme_getgenedata(int g);
	
/* model optimisation */
void gme_calparameters();
	
/* calculate percentiles of expression levels */
//void calexpression();
	
void gme_allocatemem();
void gme_freemem();
double  ** allocate_matrix(int ,int ); 
//double  ** allocate_matrix_PM(int ,int );
//double ** dotproduct( double *,double *,int ,int );
SEXP gme_c(SEXP PMmat, SEXP GTmat, SEXP PNmat, SEXP GNmat, SEXP ANmat,SEXP totalgene, SEXP saveflag, SEXP eps);

