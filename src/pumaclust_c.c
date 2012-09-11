/* **************************************************************************** */
/*        pumaclust C implementation                                          */
/* **************************************************************************** */
#include <math.h>
#include "o8para.h"
#include "pumaclust_c.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

static clusterparam in_param;
void donlp2(void);

typedef void (*func_void_void_t)(void);
typedef void (*func_void_type_liste_donlp_x_err_t)(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
typedef void (*func_void_liste_shift_donlp_x_grad_t)(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
typedef void (*func_void_donlp_x_fx_t)(DDOUBLE donlp2_x[],DDOUBLE *fx);
typedef void (*func_void_donlp_x_gradf_t)(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
typedef void (*func_void_mode_t)(IINTEGER mode);
typedef void (*func_void_t)();

extern func_void_type_liste_donlp_x_err_t econ;
extern func_void_liste_shift_donlp_x_grad_t econgrad;
extern func_void_donlp_x_fx_t ef;
extern func_void_donlp_x_gradf_t egradf;
extern func_void_mode_t eval_extern;
extern func_void_t freemem;
extern func_void_t initialparams;
extern func_void_void_t setup;
extern func_void_void_t solchk;
extern func_void_void_t user_init;
extern func_void_void_t user_init_size;
extern func_void_t allocatemem;

void econ_pumaclust(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_pumaclust(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_pumaclust(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_pumaclust(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_pumaclust(IINTEGER mode);
void freemem_pumaclust();
void initialparams_pumaclust();
void setup_pumaclust();
void solchk_pumaclust();
void user_init_pumaclust(void);
void user_init_size_pumaclust(void);
void allocatemem_pumaclust();

void initialparams_pumaclust()
{
	in_param.centers = NULL;
	in_param.clsig = NULL;
	in_param.data_m = NULL;
	in_param.data_var = NULL;
	in_param.pjPerData = NULL;
	in_param.pj = NULL;
	in_param.CIndex = NULL;
	in_param.ClusterCenters = NULL;
	in_param.ClusterSig = NULL;
	in_param.LikeliPerGene = NULL;

	econ = econ_pumaclust;
	econgrad = econgrad_pumaclust;
	ef = ef_pumaclust;
	egradf = egradf_pumaclust;
	eval_extern = eval_extern_pumaclust;
	freemem = freemem_pumaclust;
	initialparams = initialparams_pumaclust;
	setup = setup_pumaclust;
	solchk = solchk_pumaclust;
	user_init = user_init_pumaclust;
	user_init_size = user_init_size_pumaclust;
	allocatemem = allocatemem_pumaclust;
}

void allocatemem_pumaclust()
{
	int i;
	in_param.pj = (double*)R_alloc(in_param.clusters, sizeof(double));
	in_param.pjPerData = (double**)R_alloc(in_param.genes, sizeof(double*));
	for(i=0; i<in_param.genes; i++)
		in_param.pjPerData[i] = (double*)R_alloc(in_param.clusters, sizeof(double));
}
  
void freemem_pumaclust()
{
	int i;
	for (i=0; i<in_param.genes; i++)
		if (in_param.pjPerData[i] != NULL) Free(in_param.pjPerData[i]);
	if (in_param.pjPerData != NULL) Free(in_param.pjPerData);
	if (in_param.pj != NULL) Free(in_param.pj);
}

/*double mean(double *x, int len)
{
	int i;
	double s=0.0;
	
	for (i=0; i<len; i++)
		s += x[i];
	return s/len;
}

double var(double *x, int len)
{
	int i;
	double m,s=0.0;
	
	m = mean(x,len);
	for (i=0; i<len; i++)
		s += (x[i]-m)*(x[i]-m);
		
	return s/(len-1);
}

double fmaxn(double *x, int n)
{
	int i;
	double m;
	
	m = x[0];
	for (i=1; i<n; i++)
		if (x[i]>m)
			m = x[i];
			
	return m;
}
*/

void calMU(double **mujd, double *clsig, int nstart)
{
	int i, j, k;
	double t1, t2, expr_temp, var_temp;
	
	for (k=0; k<in_param.clusters; k++)
	{
		for (j=0; j<in_param.chips; j++)
		{
			t1 = 0.0;
			t2 = 0.0;
			for (i=0; i<in_param.genes; i++)
			{
				expr_temp = in_param.data_m[j*in_param.genes+i];
				var_temp = in_param.data_var[j*in_param.genes+i];
				t1 += in_param.pjPerData[i][k]*expr_temp/(clsig[k+nstart]+var_temp);
				t2 += in_param.pjPerData[i][k]/(clsig[k+nstart]+var_temp);
			}
			mujd[k][j] = t1/t2;
		}
	}
}


void workout()
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
	int i, j, k, iter;
    	double foptold, fopt, *loga=NULL, *exprs=NULL, *vars=NULL, temp;
        /*FILE *pf=NULL;*/

/*        pf = fopen("lamda_noninfprior.txt","wt");*/
	
	loga = (double*)R_alloc(in_param.clusters, sizeof(double));
	exprs = (double*)R_alloc(in_param.chips, sizeof(double));
	vars = (double*)R_alloc(in_param.chips, sizeof(double));

	for (i=0; i<in_param.clusters; i++)
		in_param.pj[i] = 1.0/in_param.clusters;
	iter = 0;

	foptold = POSI_INF;
	fopt = 1.0e18;
	while (foptold-fopt > abs(in_param.eps*fopt))
	{
		/* E-step */
		for (i=0; i<in_param.genes; i++)
		{
			for (j=0; j<in_param.chips; j++)
			{
				exprs[j] = in_param.data_m[j*in_param.genes+i];
				vars[j] = in_param.data_var[j*in_param.genes+i];
			}
			for (k=0; k<in_param.clusters; k++)
			{
				loga[k] = log(in_param.pj[k]);
				for (j=0; j<in_param.chips; j++)
					loga[k] -= ((exprs[j]-in_param.centers[j*in_param.clusters+k])*(exprs[j]-in_param.centers[j*in_param.clusters+k])/(in_param.clsig[k]+vars[j])+log(in_param.clsig[k]+vars[j]))/2.0;
				loga[k] = exp(loga[k]);
			}
			
			temp = in_param.clusters*mean(loga,in_param.clusters);
			for (k=0; k<in_param.clusters; k++)
				in_param.pjPerData[i][k] = loga[k]/temp;
		}

		/* M-step */
		for (k=0; k<in_param.clusters; k++)
		{
			temp = 0.0;
			for (i=0; i<in_param.genes; i++) temp += in_param.pjPerData[i][k];
			in_param.pj[k] = temp/in_param.genes;
		}
		
		foptold = fopt;
		
		donlp2();
		fopt = fx;
		iter++;
		Rprintf(".");
		/*Rprintf("%d  %d %d %f %d\n", iter, icf, icgf, fx, (int)optite+11);*/
	}
	Rprintf("\n");

	for (i=0; i<in_param.genes; i++)
	{
		in_param.CIndex[i] = 1;
		for (k=1; k<in_param.clusters; k++)
			if (in_param.pjPerData[i][k]>in_param.pjPerData[i][in_param.CIndex[i]-1])
				in_param.CIndex[i] = k+1;
	}	

	for (k=0; k<in_param.clusters; k++)
	{
		for (j=0; j<in_param.chips; j++)
		{
			in_param.ClusterCenters[j*in_param.clusters+k] = in_param.centers[j*in_param.clusters+k];
		}
		in_param.ClusterSig[k] = in_param.clsig[k];
	} 

	for (i=0; i<in_param.genes; i++)
	{
		for (k=0; k<in_param.clusters; k++)
			in_param.LikeliPerGene[k*in_param.genes+i] = in_param.pjPerData[i][k];
	}	

	*(in_param.BIC) = 2.0*fopt+((2+in_param.chips)*in_param.clusters-1)*log(in_param.genes);
        /*fprintf(pf,"%f\n",sqrt(1.0/in_param.lamda_m));*/
	
        /*fclose(pf);*/
/*	if (loga != NULL) Free(loga);
	if (exprs != NULL) Free(exprs);
	if (vars != NULL) Free(vars);
*/}

SEXP pumaclust_c(SEXP Mmat, SEXP Stdmat, SEXP clusters, SEXP centers, SEXP clsig, SEXP eps, SEXP del0)
{
	SEXP dim=NULL;	
	SEXP res=NULL, CIndex=NULL, Centers=NULL, ClusterSig=NULL, pjPerData=NULL, BIC=NULL;
         
	initialparams_pumaclust();
	
	PROTECT(dim = getAttrib(Mmat, R_DimSymbol));
	in_param.genes = INTEGER(dim)[0];
	in_param.chips = INTEGER(dim)[1];
	
	in_param.clusters = INTEGER_POINTER(AS_INTEGER(clusters))[0];
	in_param.centers = NUMERIC_POINTER(AS_NUMERIC(centers));
	in_param.clsig = NUMERIC_POINTER(AS_NUMERIC(clsig));

	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];
	in_param.del0 = NUMERIC_POINTER(AS_NUMERIC(del0))[0];
	
	in_param.data_m = NUMERIC_POINTER(AS_NUMERIC(Mmat));
	in_param.data_var = NUMERIC_POINTER(AS_NUMERIC(Stdmat));
	
	allocatemem_pumaclust();
          
	PROTECT(CIndex = allocVector(INTSXP, in_param.genes));
	PROTECT(Centers = allocMatrix(REALSXP, in_param.clusters, in_param.chips));
	PROTECT(ClusterSig = allocVector(REALSXP, in_param.clusters));
	PROTECT(pjPerData = allocMatrix(REALSXP, in_param.genes, in_param.clusters));
	PROTECT(BIC = allocVector(REALSXP,1));
	PROTECT(res = allocVector(VECSXP,5));	
	in_param.CIndex = INTEGER_POINTER(AS_INTEGER(CIndex));
	in_param.ClusterCenters = NUMERIC_POINTER(AS_NUMERIC(Centers));
	in_param.ClusterSig = NUMERIC_POINTER(AS_NUMERIC(ClusterSig));
	in_param.LikeliPerGene = NUMERIC_POINTER(AS_NUMERIC(pjPerData));
	in_param.BIC = NUMERIC_POINTER(AS_NUMERIC(BIC));
	
	Rprintf("Clustering is performing ");
	workout();  

	SET_VECTOR_ELT(res, 0, CIndex);
	SET_VECTOR_ELT(res, 1, Centers);
	SET_VECTOR_ELT(res, 2, ClusterSig);
	SET_VECTOR_ELT(res, 3, pjPerData);
	SET_VECTOR_ELT(res, 4, BIC);
/*	freemem_pumaclust();
*/	Rprintf("Done.\n");
	UNPROTECT(7);
	
	return res;
}

/* **************************************************************************** */
/*                              donlp2-intv size initialization                           */
/* **************************************************************************** */
void user_init_size_pumaclust(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X


    /* problem dimension n = dim(donlp2_x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
      
    n = in_param.clusters;
    nstep = 20;
    
    nlin   =  0;
    nonlin =  0;
    iterma = 4000;
}

/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_pumaclust(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X
    
    static IINTEGER j;

    silent = TRUE;
/*    intakt = TRUE;*/

    big = 1.e20;
	/* initialise the parameters */
    for (j=1; j<=in_param.clusters; j++)
    {
	donlp2_x[j] = in_param.clsig[j-1];
		
	low[j] = LOWBOUND;
	up[j] = big;
    }
 
    analyt = TRUE;
    epsdif = 1.e-16;  
    
    nreset = n;
    
    del0 = in_param.del0;
    tau0 = 1.0e0;
    tau  = 0.1e0;

    return;
}

/* **************************************************************************** */
/*                                 special setup                                */
/* **************************************************************************** */
void setup_pumaclust(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
   
    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_pumaclust(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"

	int i, j, k;
	double t1, t2, **mujd=NULL, expr_temp, var_temp;
	
	mujd = (double**)R_alloc(in_param.clusters, sizeof(double*));
	for (i=0; i<in_param.clusters; i++)
		mujd[i] = (double*)R_alloc(in_param.chips, sizeof(double));
	
	for (k=0; k<in_param.clusters; k++) in_param.clsig[k] = donlp2_x[k+1];
	calMU(mujd, in_param.clsig, 0);
	
	for (k=0; k<in_param.clusters; k++)
	    for (j=0; j<in_param.chips; j++)
	        in_param.centers[j*in_param.clusters+k] = mujd[k][j];
    
/*	for (i=0; i<in_param.clusters; i++) Free(mujd[i]);
	Free(mujd);
*/
    return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_pumaclust(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	int i, j, k;
	double t1, t2, **mujd=NULL, expr_temp, var_temp;
	
	mujd = (double**)R_alloc(in_param.clusters, sizeof(double*));
	for (i=0; i<in_param.clusters; i++)
		mujd[i] = (double*)R_alloc(in_param.chips, sizeof(double));


	calMU(mujd, donlp2_x, 1);
	
	*fx = 0.0;
	for (i=0; i<in_param.genes; i++)
	{
		for (k=0; k<in_param.clusters; k++)
		{
			t1 = 0.0;
			t2 = 0.0;
			for (j=0; j<in_param.chips; j++)
			{
				expr_temp = in_param.data_m[j*in_param.genes+i];
				var_temp = in_param.data_var[j*in_param.genes+i];
				
				t1 += log(donlp2_x[k+1]+var_temp);
				t2 += (expr_temp-mujd[k][j])*(expr_temp-mujd[k][j])/(donlp2_x[k+1]+var_temp);
			}
		
			*fx += in_param.pjPerData[i][k]*(log(in_param.pj[k])-in_param.chips*log(2.0*M_PI)/2.0-t1/2.0-t2/2.0);
		}	
	}
	*fx = -*fx;

/*	for (i=0; i<in_param.clusters; i++) Free(mujd[i]);
	Free(mujd);
*/	
    return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_pumaclust(DDOUBLE donlp2_x[],DDOUBLE gradf[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	int i, j, k;
	double t1, t2, **mujd=NULL, expr_temp, var_temp;
	
	mujd = (double**)R_alloc(in_param.clusters, sizeof(double*));
	for (i=0; i<in_param.clusters; i++)
		mujd[i] = (double*)R_alloc(in_param.chips, sizeof(double));
	
	calMU(mujd, donlp2_x, 1);
	for (k=0; k<in_param.clusters; k++)
	{
		gradf[k+1] = 0.0;
		for (i=0; i<in_param.genes; i++)
		{
			t1 = 0.0;
			t2 = 0.0;
			for (j=0; j<in_param.chips; j++)
			{
				expr_temp = in_param.data_m[j*in_param.genes+i];
				var_temp = in_param.data_var[j*in_param.genes+i];
				
				t1 += 1.0/(donlp2_x[k+1]+var_temp);
				t2 += (expr_temp-mujd[k][j])*(expr_temp-mujd[k][j])/((donlp2_x[k+1]+var_temp)*(donlp2_x[k+1]+var_temp));
			}
			gradf[k+1] += -in_param.pjPerData[i][k]*(-t1/2.0+t2/2.0);
		}
	}

/*	for (i=0; i<in_param.clusters; i++) Free(mujd[i]);
	Free(mujd);
*/
    return;
}

/* **************************************************************************** */
/*  no nonlinear constraints */
/* **************************************************************************** */
void econ_pumaclust(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[], 
              LLOGICAL err[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_pumaclust(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_pumaclust(IINTEGER mode) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}
