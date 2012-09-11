/* **************************************************************************** */
/*        multi-mgMOS C implementation                                          */
/* **************************************************************************** */
#include <math.h>
#include "o8para.h"
#include "pplr_c.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

static pplrparam in_param;
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

void econ_pplr(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_pplr(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_pplr(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_pplr(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_pplr(IINTEGER mode);
void freemem_pplr();
void initialparams_pplr();
void setup_pplr();
void solchk_pplr();
void user_init_pplr(void);
void user_init_size_pplr(void);
void allocatemem_pplr();

void initialparams_pplr()
{
	in_param.replicates = NULL;
	in_param.data_m = NULL;
	in_param.data_var = NULL;
	in_param.expr = NULL;
	in_param.var = NULL;
	in_param.outp = NULL;
	in_param.mu1 = NULL;
	in_param.mu1sq = NULL;

	econ = econ_pplr;
	econgrad = econgrad_pplr;
	ef = ef_pplr;
	egradf = egradf_pplr;
	eval_extern = eval_extern_pplr;
	freemem = freemem_pplr;
	initialparams = initialparams_pplr;
	setup = setup_pplr;
	solchk = solchk_pplr;
	user_init = user_init_pplr;
	user_init_size = user_init_size_pplr;
	allocatemem = allocatemem_pplr;
}

void allocatemem_pplr()
{
	in_param.expr = (double*)R_alloc(in_param.chips, sizeof(double));
	in_param.var = (double*)R_alloc(in_param.chips, sizeof(double));
	if (in_param.method > 0)
	{
		in_param.mu1 = (double*)R_alloc(in_param.chips, sizeof(double));
		in_param.mu1sq = (double*)R_alloc(in_param.chips, sizeof(double));
	}
}

  
void freemem_pplr()
{
	if (in_param.expr != NULL) Free(in_param.expr);
	if (in_param.var != NULL) Free(in_param.var);
	if (in_param.method > 0)
	{
		if (in_param.mu1 != NULL) Free(in_param.mu1);
		if (in_param.mu1sq != NULL) Free(in_param.mu1sq);
	}
}

void findeforc(double *x, int c, double *e_c, int *num_c)
{
	int i,j;
	
	j = 0;
	for (i=0; i<in_param.chips; i++)
	{
		if (in_param.replicates[i] == c)
		{
			e_c[j] = x[i];
			j++;
		}
	}
	*num_c = j;
}

double mean(double *x, int len)
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

void workout_comb0()
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
	int i, j;
     
	for (i=0; i<in_param.genes; i++)
	{
		for (j=0; j<in_param.chips; j++)
		{
			in_param.expr[j] = in_param.data_m[j*in_param.genes+i];
			in_param.var[j] = in_param.data_var[j*in_param.genes+i];
		}
		in_param.cur_gene = i;
			
		donlp2();
			
		/*Rprintf("%d  %d %d %f %d\n", i+1, icf, icgf, fx, (int)optite+11);*/
		/*if ((int)i%200 == 0)
			Rprintf(".");*/
	}
}

void workout_comb1()
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
	int i, j, k,num_c;
    	double exp_c[MAX_NUM_REPLICATE], var_c[MAX_NUM_REPLICATE], x_temp[MAX_NUM_COND];
    	double v[MAX_NUM_COND], temp, tempx, foptold, fopt, lamda10[MAX_NUM_COND];
    	double mu_temp[MAX_NUM_COND], lam_temp[MAX_NUM_COND];
	double mu1old[MAX_NUM_COND], lamda1[MAX_NUM_COND], lamda1old[MAX_NUM_COND];
	double alpha1, alpha2, beta2, diff_mu1[MAX_NUM_COND], diff_lamda[MAX_NUM_COND];
	double *sample_x=NULL, *logy=NULL, *w_param=NULL;
	double beta1[MAX_NUM_CHIP], f_lamda[MAX_NUM_CHIP];
	double f_lamda_c[MAX_NUM_REPLICATE];
        /*FILE *pf=NULL;*/

/*        pf = fopen("lamda_noninfprior.txt","wt");*/
	
          sample_x = (double*)R_alloc(in_param.sample_num, sizeof(double));
	logy = (double*)R_alloc(in_param.sample_num, sizeof(double));
	w_param = (double*)R_alloc(in_param.sample_num, sizeof(double));
	
	GetRNGstate();
	for (i=0; i<in_param.genes; i++)
	{
		for (j=0; j<in_param.chips; j++)
		{
			in_param.expr[j] = in_param.data_m[j*in_param.genes+i];
			in_param.var[j] = in_param.data_var[j*in_param.genes+i];
		}
		in_param.cur_gene = i;
    		
		for (j=1; j<=in_param.conds; j++)
		{
			findeforc(in_param.expr, j, exp_c, &num_c);
			findeforc(in_param.var, j, var_c, &num_c);
			
			mu_temp[j-1] = mean(exp_c, num_c);
			if (num_c == 1)
				temp = 0.1;
			else
				temp = var(exp_c, num_c);
			if (temp == 0.0)
			    temp = 0.1;
			    
			lamda10[j-1] = 1/temp;
			tempx = 0.0;
			for (k=0; k<num_c; k++)
			{
				v[k] = fmax2(1.0e-4, temp-var_c[k]);
				if (v[k]>tempx)
					tempx = v[k];
			}
			lam_temp[j-1] = 1/mean(v, num_c);
		}
			
		in_param.mu0 = mean(mu_temp, in_param.conds);
		temp = var(mu_temp, in_param.conds);
		if (temp == 0.0)
		    temp = 0.1;
		in_param.eta0 = 1/temp;
                 
                   
		if (tempx<=1.0e-3)
		{
			in_param.alpha0 = 10.0;
			in_param.beta0 = 1.0;
		}
		else
		{
			temp = mean(lam_temp, in_param.conds);
			tempx = var(lam_temp, in_param.conds);
			if (tempx == 0.0)
			    tempx = 0.1;
			in_param.beta0 = temp/tempx;
			in_param.alpha0 = temp*temp/tempx;
		}
		
                   

		foptold = POSI_INF;
		fopt = 1.0e18;
		while (foptold-fopt > in_param.eps)
		{
			for (j=0; j<in_param.conds; j++)
			{
				in_param.mu1[j] = mu_temp[j];
				lamda1[j] = lamda10[j];
				mu1old[j] = POSI_INF;
				lamda1old[j] = POSI_INF;
				diff_mu1[j] = fabs(in_param.mu1[j]-mu1old[j]);
				diff_lamda[j] = fabs(lamda1[j]-lamda1old[j]);
			}
			
			alpha1 = 3.0/2.0;
			alpha2 = in_param.alpha0;
			beta2 = in_param.beta0;
			
			for (j=0; j<in_param.sample_num; j++)
				sample_x[j] = rgamma(alpha2, 1.0/beta2);
			
			while (fmaxn(diff_mu1, in_param.conds)>in_param.eps || fmaxn(diff_lamda, in_param.conds)>in_param.eps)
			{
				for (j=0; j<in_param.conds; j++)
					in_param.mu1sq[j] = in_param.mu1[j]*in_param.mu1[j]+1.0/lamda1[j];
				for (k=0; k<in_param.sample_num; k++)
					logy[k] = 0.0;
				
				for (j=0; j<in_param.chips; j++)
				{
					beta1[j] = (in_param.expr[j]*in_param.expr[j]-2.0*in_param.mu1[in_param.replicates[j]-1]*in_param.expr[j]+in_param.mu1sq[in_param.replicates[j]-1])/2.0;
					for (k=0; k<in_param.sample_num; k++)
						 logy[k] += dgamma(1.0/(1.0/sample_x[k]+in_param.var[j]), alpha1, 1.0/beta1[j], TRUE);
					
					f_lamda[j] = 0.0;
				}
				tempx = fmaxn(logy, in_param.sample_num);
				temp = 0.0;
				for (k=0; k<in_param.sample_num; k++)
					temp += exp(logy[k]-tempx);
				in_param.lamda_m = 0.0;
				for (k=0; k<in_param.sample_num; k++)
				{
					w_param[k] = exp(logy[k]-tempx-log(temp));
					in_param.lamda_m += w_param[k]*sample_x[k];
					for (j=0; j<in_param.chips; j++)
						f_lamda[j] += w_param[k]/(1.0/sample_x[k]+in_param.var[j]);
				}
				
				for (j=0; j<in_param.conds; j++)
				{
					lamda1old[j] = lamda1[j];
					mu1old[j] = in_param.mu1[j];
					findeforc(f_lamda, j+1, f_lamda_c, &num_c);
					findeforc(in_param.expr, j+1, exp_c, &num_c);
					
					lamda1[j] = in_param.eta0;
					for (k=0; k<num_c; k++)
						lamda1[j] += f_lamda_c[k];
					in_param.mu1[j] = in_param.eta0*in_param.mu0/lamda1[j];
					for (k=0; k<num_c; k++)
						in_param.mu1[j] += f_lamda_c[k]*exp_c[k]/lamda1[j];
					
					diff_mu1[j] = fabs(in_param.mu1[j]-mu1old[j]);
					diff_lamda[j] = fabs(lamda1[j]-lamda1old[j]);
				}
			}
			in_param.loglamda = 0.0;
			for (k=0; k<in_param.sample_num; k++)
				in_param.loglamda += w_param[k]*log(sample_x[k]);
				
			foptold = fopt;
			
			donlp2();
			
			fopt = fx;
                              
                             
                             
		
		}
		/*Rprintf("%d  %d %d %f %d\n", i+1, icf, icgf, fx, (int)optite+11);*/
                   	


		for (j=0; j<in_param.conds; j++)
		{
			/* mean of expression level */
			in_param.outp[j*in_param.genes+i] = in_param.mu1[j];
			/* standard deviation of expression level */
			in_param.outp[(in_param.conds+j)*in_param.genes+i] = sqrt(1.0/lamda1[j]);
		}	
                /*fprintf(pf,"%f\n",sqrt(1.0/in_param.lamda_m));*/
	
		/*if ((int)i%200 == 0)
			Rprintf(".");*/
	}
	PutRNGstate();
        /*fclose(pf);*/
/*	Free(sample_x);
	Free(logy);
	Free(w_param);
*/}

SEXP bcomb_c(SEXP Mmat, SEXP Stdmat, SEXP rep, SEXP method, SEXP conds, SEXP nsample, SEXP eps)
{
	SEXP dim=NULL;	
	SEXP res=NULL;

	initialparams_pplr();
	
	PROTECT(dim = getAttrib(Mmat, R_DimSymbol));
	in_param.genes = INTEGER(dim)[0];
	in_param.chips = INTEGER(dim)[1];
	
	in_param.conds = INTEGER_POINTER(AS_INTEGER(conds))[0];
	in_param.method = INTEGER_POINTER(AS_INTEGER(method))[0];
	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];
	
	in_param.data_m = NUMERIC_POINTER(AS_NUMERIC(Mmat));
	in_param.data_var = NUMERIC_POINTER(AS_NUMERIC(Stdmat));
	
	in_param.replicates = INTEGER_POINTER(AS_INTEGER(rep));
	in_param.sample_num = INTEGER_POINTER(AS_INTEGER(nsample))[0];

	allocatemem_pplr();
	PROTECT(res = allocMatrix(REALSXP, in_param.genes, 2*in_param.conds));
	in_param.outp = NUMERIC_POINTER(AS_NUMERIC(res));
	
	/*Rprintf("Combined signals are computing ");*/
	if (in_param.method == 4)
		workout_comb0();   /* non-informative prior, shared sig^2 */
	else
		workout_comb1();   /* conjugate prior, shared sig^2 */

/*	freemem_pplr();
*/	/*Rprintf("Done.\n");*/
	UNPROTECT(2);
	
	return res;
}

/* **************************************************************************** */
/*                              donlp2-intv size initialization                           */
/* **************************************************************************** */
void user_init_size_pplr(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X


    /* problem dimension n = dim(donlp2_x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
      
    if (in_param.method == 1)
    {
    	n = 4;
        nstep = 20;
    }
    else
    {
    	n = in_param.conds+5;
        nstep = 20;
    }
    
    nlin   =  0;
    nonlin =  0;
    iterma = 4000;
}

/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_pplr(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X
    
    static IINTEGER i,j, num_c;
    double exp_c[MAX_NUM_REPLICATE], var_c[MAX_NUM_REPLICATE], x_temp[MAX_NUM_COND];
    double v[MAX_NUM_COND];
    double mu_temp[MAX_NUM_COND], lam_temp[MAX_NUM_COND];

    silent = TRUE;
/*    intakt = TRUE;*/
    
    big = 1.e20;
    if (in_param.method==4)
    {
	/* initialise the parameters */
	double v_temp[MAX_NUM_COND], temp1, temp2;
	for (j=1; j<=in_param.conds; j++)
	{
		findeforc(in_param.expr, j, exp_c, &num_c);
		donlp2_x[j] = mean(exp_c, num_c);
		x_temp[j-1] = donlp2_x[j];
		if (num_c == 1)
			v[j-1] = 0.1;
		else
			v[j-1] = var(exp_c, num_c);
		if (v[j-1]==0.0)
		    v[j-1] = 0.1;
		v_temp[j-1] = 1/v[j-1];
	}
	donlp2_x[in_param.conds+1] = mean(v, in_param.conds);
	donlp2_x[in_param.conds+2] = mean(x_temp, in_param.conds);
	donlp2_x[in_param.conds+3] = var(x_temp, in_param.conds);
	temp1 = mean(v_temp, in_param.conds);
	temp2 = var(v_temp, in_param.conds);
	if (temp2==0.0)
	    temp2 = 0.1;
	donlp2_x[in_param.conds+4] = temp1*temp1/temp2;
	donlp2_x[in_param.conds+5] = temp1/temp2;
		
	for (i=1; i<=in_param.conds+5; i++)
	{
		low[i] = LOWBOUND;
		up[i] = big;
		/*Rprintf("%f ", donlp2_x[i]);*/
	}
	/*Rprintf("\n");*/
    }
    else
    {
    	donlp2_x[1] = in_param.eta0;
    	donlp2_x[2] = in_param.mu0;
    	donlp2_x[3] = in_param.alpha0;
    	donlp2_x[4] = in_param.beta0;
	
	for (i=1; i<=4; i++)
	{
		low[i] = LOWBOUND;
		up[i] = big;
	}
	low[2] = -big;
    	
    }

    analyt = TRUE;
    epsdif = 1.e-16;  
    
    nreset = n;
    
    del0 = 0.2e0;
    tau0 = 1.0e0;
    tau  = 0.1e0;
       
    return;
}

/* **************************************************************************** */
/*                                 special setup                                */
/* **************************************************************************** */
void setup_pplr(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
   
    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_pplr(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"

    IINTEGER i, j, num_c;
    double sig_ret, mu_ret, tau_ret;
    double t1,t2;
    double expr_c[MAX_NUM_REPLICATE], var_c[MAX_NUM_REPLICATE];
    
    if (in_param.method == 4)
    {
    	sig_ret = donlp2_x[in_param.conds+1];
	mu_ret = donlp2_x[in_param.conds+2];
	tau_ret = donlp2_x[in_param.conds+3];
	for (j=0; j<in_param.conds; j++)
	{
		findeforc(in_param.expr, j+1, expr_c, &num_c);
		findeforc(in_param.var, j+1, var_c, &num_c);
		t1 = mu_ret/tau_ret;
		t2 = 1/tau_ret;
		for (i=0; i<num_c; i++)
		{
			t1 += expr_c[i]/(sig_ret+var_c[i]);
			t2 += 1.0/(sig_ret+var_c[i]);
		}
		/* mean of expression level */
		in_param.outp[j*in_param.genes+in_param.cur_gene] = t1/t2;
		/* standard deviation of expression level */
		in_param.outp[(in_param.conds+j)*in_param.genes+in_param.cur_gene] = sqrt(1.0/t2);
	}
	/*Rprintf("%f %f %f %f %f\n", donlp2_x[1],donlp2_x[2],sig_ret,mu_ret,tau_ret);*/
    }
    else
    {
    	in_param.eta0 = donlp2_x[1];
    	in_param.mu0 = donlp2_x[2];
    	in_param.alpha0 = donlp2_x[3];
    	in_param.beta0 = donlp2_x[4];
    }
    return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_pplr(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	int i, j;
	
	*fx = 0.0;
	
        if (in_param.method == 1)
	{
		double eta_x, u_x, alpha_x, beta_x;
		
		eta_x = donlp2_x[1];
		u_x = donlp2_x[2];
		alpha_x = donlp2_x[3];
		beta_x = donlp2_x[4];
		
		*fx = in_param.conds*log(eta_x)/2.0+alpha_x*in_param.loglamda-beta_x*in_param.lamda_m+alpha_x*log(beta_x)-lgammafn(alpha_x);
		for (i=0; i<in_param.conds; i++)
			*fx -= eta_x*(in_param.mu1sq[i]-2.0*in_param.mu1[i]*u_x+u_x*u_x)/2.0;
		*fx = -*fx;
	}
	else
	{
		double sig_x, mu_x, tau_x, t1, t2, alpha_x, beta_x;
		double expr_c[MAX_NUM_REPLICATE], var_c[MAX_NUM_REPLICATE];
		int num_c;
	 	
		sig_x = donlp2_x[in_param.conds+1];
		mu_x = donlp2_x[in_param.conds+2];
		tau_x = donlp2_x[in_param.conds+3];
		alpha_x = donlp2_x[in_param.conds+4];
		beta_x = donlp2_x[in_param.conds+5];
		
		*fx = -alpha_x*log(beta_x)+(alpha_x-1.0)*log(sig_x)+beta_x/sig_x+lgammafn(alpha_x)+in_param.conds*log(tau_x)/2.0;
		t1 = 0.0;
		t2 = 0.0;
		for (j=0; j<in_param.conds; j++)
		{
			findeforc(in_param.expr, j+1, expr_c, &num_c);
			findeforc(in_param.var, j+1, var_c, &num_c);
			for (i=0; i<num_c; i++)
			{
				t1 += log(sig_x+var_c[i]);
				t2 += (expr_c[i]-donlp2_x[j+1])*(expr_c[i]-donlp2_x[j+1])
					/(2.0*(sig_x+var_c[i]));
			}
			*fx += (donlp2_x[j+1]-mu_x)*(donlp2_x[j+1]-mu_x)/(2.0*tau_x);
		}
		*fx += t1/2.0+t2;
		/*Rprintf("x: %f %f %f %f %f\n", donlp2_x[1], donlp2_x[2], sig_x, mu_x, tau_x);*/
		/*Rprintf("fx: %f\n", *fx);*/
	}
	
    return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_pplr(DDOUBLE donlp2_x[],DDOUBLE gradf[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	int i, j;
    
        if (in_param.method == 1)
	{
		double eta_x, u_x, alpha_x, beta_x;
		
		eta_x = donlp2_x[1];
		u_x = donlp2_x[2];
		alpha_x = donlp2_x[3];
		beta_x = donlp2_x[4];
		
		gradf[1] = -in_param.conds/(2.0*eta_x);
		gradf[2] = 0.0;
		for (i=0; i<in_param.conds; i++)
		{
			gradf[1] += (in_param.mu1sq[i]-2.0*in_param.mu1[i]*u_x+u_x*u_x)/2.0;
			gradf[2] -= (in_param.mu1[i]-u_x)*eta_x;
		}
		gradf[3] = -in_param.loglamda-log(beta_x)+digamma(alpha_x);
		gradf[4] = in_param.lamda_m-alpha_x/beta_x;
	}
	else
	{
		double sig_x, mu_x, tau_x, alpha_x, beta_x;
		double expr_c[MAX_NUM_REPLICATE], var_c[MAX_NUM_REPLICATE];
		int num_c;
	
		sig_x = donlp2_x[in_param.conds+1];
		mu_x = donlp2_x[in_param.conds+2];
		tau_x = donlp2_x[in_param.conds+3];
		alpha_x = donlp2_x[in_param.conds+4];
		beta_x = donlp2_x[in_param.conds+5];
		
		gradf[in_param.conds+1] = -(alpha_x-1.0)/sig_x+beta_x/(sig_x*sig_x);
		gradf[in_param.conds+2] = 0.0;
		gradf[in_param.conds+3] = in_param.conds/(2.0*tau_x);
		for (j=0; j<in_param.conds; j++)
		{
			gradf[j+1] = (donlp2_x[j+1]-mu_x)/tau_x;
			findeforc(in_param.expr, j+1, expr_c, &num_c);
			findeforc(in_param.var, j+1, var_c, &num_c);
			for (i=0; i<num_c; i++)
			{
				gradf[j+1] -= (expr_c[i]-donlp2_x[j+1])/(sig_x+var_c[i]);
				gradf[in_param.conds+1] += 1.0/(2.0*(sig_x+var_c[i]))-
						(expr_c[i]-donlp2_x[j+1])*(expr_c[i]-donlp2_x[j+1])/(2.0*(sig_x+var_c[i])*(sig_x+var_c[i]));
			}
			gradf[in_param.conds+2] -= (donlp2_x[j+1]-mu_x)/tau_x;
			gradf[in_param.conds+3] -= (donlp2_x[j+1]-mu_x)*(donlp2_x[j+1]-mu_x)/(2.0*tau_x*tau_x);
		}
		gradf[in_param.conds+4] = -log(beta_x)+log(sig_x)+digamma(alpha_x);
		gradf[in_param.conds+5] = -alpha_x/beta_x+1/sig_x;
		/*Rprintf("gradf: %f %f %f %f %f\n", gradf[1], gradf[2], gradf[3], gradf[4], gradf[5]);*/
	}

    return;
}

/* **************************************************************************** */
/*  no nonlinear constraints */
/* **************************************************************************** */
void econ_pplr(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[], 
              LLOGICAL err[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_pplr(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_pplr(IINTEGER mode) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}
