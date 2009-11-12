/* **************************************************************************** */
/*        pumaclustii C implementation                                          */
/* **************************************************************************** */
#include <math.h>
#include "o8para.h"
#include "pumaclustii_c.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

static clustiiparam in_param;
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

void econ_pumaclustii(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_pumaclustii(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_pumaclustii(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_pumaclustii(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_pumaclustii(IINTEGER mode);
void freemem_pumaclustii();
void initialparams_pumaclustii();
void setup_pumaclustii();
void solchk_pumaclustii();
void user_init_pumaclustii(void);
void user_init_size_pumaclustii(void);
void allocatemem_pumaclustii();

void initialparams_pumaclustii()
{
	in_param.cond_n = NULL;
	in_param.reps = NULL;
	in_param.cmu = NULL;
	in_param.csig = NULL;
	in_param.data_m = NULL;
	in_param.data_var = NULL;
	in_param.v = NULL;
	in_param.calpha = NULL;
	in_param.cbeta = NULL;
	in_param.q_kn = NULL;
	in_param.un = NULL;
	in_param.logun = NULL;
	in_param.etan = NULL;
	in_param.logetan = NULL;

	econ = econ_pumaclustii;
	econgrad = econgrad_pumaclustii;
	ef = ef_pumaclustii;
	egradf = egradf_pumaclustii;
	eval_extern = eval_extern_pumaclustii;
	freemem = freemem_pumaclustii;
	initialparams = initialparams_pumaclustii;
	setup = setup_pumaclustii;
	solchk = solchk_pumaclustii;
	user_init = user_init_pumaclustii;
	user_init_size = user_init_size_pumaclustii;
	allocatemem = allocatemem_pumaclustii;
}


double *CallocD(long s)
{
	return (double*)R_alloc(s, sizeof(double));
}

double **CallocDD(long s)
{
	return (double**)R_alloc(s, sizeof(double*));
}

double ***CallocDDD(long s)
{
	return (double***)R_alloc(s, sizeof(double**));
}

int *CallocI(long s)
{
	return (int*)R_alloc(s, sizeof(int));
}

void FreeD(double *p)
{
	if (p != NULL) Free(p);
}

void FreeDD(double **p)
{
	if (p != NULL) Free(p);
}

void FreeDDD(double ***p)
{
	if (p != NULL) Free(p);
}

void FreeI(int *p)
{
	if (p != NULL) Free(p);
}

void allocatemem_pumaclustii()
{
	int i;
//	in_param.pj = (double*)R_allo(in_param.clusters, sizeof(double));
	in_param.cond_n = CallocI(in_param.conds);
	in_param.v = CallocD(in_param.maxcls);
	in_param.calpha = CallocD(in_param.maxcls);
	in_param.cbeta = CallocD(in_param.maxcls);
	in_param.q_kn = (double**)R_alloc(in_param.genes, sizeof(double*));
	in_param.un = (double**)R_alloc(in_param.genes,  sizeof(double*));
	in_param.logun = (double**)R_alloc(in_param.genes,  sizeof(double*));
	in_param.etan = (double**)R_alloc(in_param.genes,  sizeof(double*));
	in_param.logetan = (double**)R_alloc(in_param.genes,  sizeof(double*));
	in_param.q_kn_best = CallocDD(in_param.genes);
	for (i=0; i<in_param.genes; i++)
	{
		in_param.q_kn[i] = (double*)R_alloc(in_param.maxcls,  sizeof(double));
		in_param.un[i] = (double*)R_alloc(in_param.maxcls,  sizeof(double));
		in_param.logun[i] = (double*)R_alloc(in_param.maxcls,  sizeof(double));
		in_param.etan[i] = (double*)R_alloc(in_param.maxcls,  sizeof(double));
		in_param.logetan[i] = (double*)R_alloc(in_param.maxcls,  sizeof(double));
		in_param.q_kn_best[i] = CallocD(in_param.maxcls);
	}
	in_param.pi_k_best = CallocD(in_param.maxcls);
	in_param.cmu = CallocDD(in_param.maxcls);
	in_param.csig = CallocDD(in_param.maxcls);
	in_param.cmu_best = CallocDD(in_param.maxcls);
	in_param.csig_best = CallocDD(in_param.maxcls);
	for (i=0; i<in_param.maxcls; i++)
	{
		in_param.cmu_best[i] = CallocD(in_param.conds);
		in_param.cmu[i] = CallocD(in_param.conds);
		in_param.csig[i] = CallocD(in_param.conds);
		in_param.csig_best[i] = CallocD(in_param.conds);
	}
}
  
void freemem_pumaclustii()
{
	int i;
	FreeI(in_param.cond_n);
	for (i=0; i<in_param.genes; i++)
	{
		FreeD(in_param.q_kn[i]);
		FreeD(in_param.un[i]);
		FreeD(in_param.logun[i]);
		FreeD(in_param.etan[i]);
		FreeD(in_param.logetan[i]);
		FreeD(in_param.q_kn_best[i]);
	}
	FreeDD(in_param.q_kn);
	FreeDD(in_param.un);
	FreeDD(in_param.logun);
	FreeDD(in_param.etan);
	FreeDD(in_param.logetan);
	FreeD(in_param.v);
	FreeD(in_param.calpha);
	FreeD(in_param.cbeta);
	FreeDD(in_param.q_kn_best);
	FreeD(in_param.pi_k_best);
	for (i=0; i<in_param.maxcls; i++)
	{
		FreeD(in_param.cmu_best[i]);
		FreeD(in_param.cmu[i]);
		FreeD(in_param.csig_best[i]);
	}
	FreeDD(in_param.cmu_best);
	FreeDD(in_param.cmu);
	FreeDD(in_param.csig_best);
	//if (in_param.pj != NULL) Free(in_param.pj);
}

double maxD(double x, double y)
{
	if (x>y)	return x;
	else return y;
}

double absD(double x)
{
	if (x<0)	return -x;
	else return x;
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
}*/

void workout_pumaclustii()
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
	int i, j, k, iter;
   	double *pi_k=NULL, **A=NULL, **temp1=NULL;
	double ***wn=NULL, ***wn2=NULL, **alpha_u=NULL, **beta_u=NULL, **alpha_eta=NULL, **beta_eta=NULL;
	double **t1=NULL, *t2=NULL, ***mu_w=NULL, ***sig_w=NULL, ***mu_t=NULL, ***sig_t=NULL;
	double  *exprs=NULL, *vars=NULL, temp, *temp2=NULL, temp3, **t3=NULL, *t4=NULL;
	double *a1=NULL, *a2=NULL, *a3=NULL, *a4=NULL, *a5=NULL, *a6=NULL, *a7=NULL, *a8=NULL, *a9=NULL;
	int K_nz, K_nz_t, n_hat, iit;
	double F_max, F0, F, Ft0, Ft, max_A, F_t;
	double **q_kn_t=NULL, *pi_k_t=NULL, **cmu_t=NULL, **csig_t=NULL;
	
	/* Allocate memory */
	exprs = CallocD(in_param.chips);
	vars = CallocD(in_param.chips);

	pi_k = CallocD(in_param.maxcls);
	t2 = CallocD(in_param.maxcls);

	A = CallocDD(in_param.genes);
	wn = CallocDDD(in_param.genes);
	wn2 = CallocDDD(in_param.genes);
	q_kn_t = CallocDD(in_param.genes);
	for (i=0; i<in_param.genes; i++)
	{
		A[i] = CallocD(in_param.maxcls);
		wn[i] = CallocDD(in_param.maxcls);
		wn2[i] = CallocDD(in_param.maxcls);
		for (j=0; j<in_param.maxcls; j++)
		{
			wn[i][j] = CallocD(in_param.conds);
			wn2[i][j] = CallocD(in_param.conds);
		}
		q_kn_t[i] = CallocD(in_param.maxcls);
	}	

	t1 = CallocDD(in_param.maxcls);
	for (i=0; i<in_param.maxcls; i++)
	{
		t1[i] = CallocD(in_param.conds);
	}
	
	mu_w = CallocDDD(in_param.genes);
	sig_w = CallocDDD(in_param.genes);
	mu_t = CallocDDD(in_param.genes);
	sig_t = CallocDDD(in_param.genes);
	alpha_u = CallocDD(in_param.genes);
	beta_u = CallocDD(in_param.genes);
	alpha_eta = CallocDD(in_param.genes);
	beta_eta = CallocDD(in_param.genes);
	
	for (i=0; i<in_param.genes; i++)
	{
		mu_w[i] = CallocDD(in_param.maxcls);
		sig_w[i] = CallocDD(in_param.maxcls);
		mu_t[i] = CallocDD(in_param.maxcls);
		sig_t[i] = CallocDD(in_param.maxcls);	
		for (j=0; j<in_param.maxcls; j++)
		{
			mu_w[i][j] = CallocD(in_param.conds);
			sig_w[i][j] = CallocD(in_param.conds);
			mu_t[i][j] = CallocD(in_param.chips);
			sig_t[i][j] = CallocD(in_param.chips);
		}
		alpha_u[i] = CallocD(in_param.maxcls);
		beta_u[i] = CallocD(in_param.maxcls);
		alpha_eta[i] = CallocD(in_param.maxcls);
		beta_eta[i] = CallocD(in_param.maxcls);
	}
	
	a1 = CallocD(in_param.maxcls);
	a2 = CallocD(in_param.maxcls);
	a3 = CallocD(in_param.maxcls);
	a4 = CallocD(in_param.maxcls);
	a5 = CallocD(in_param.maxcls);
	a6 = CallocD(in_param.maxcls);
	a7 = CallocD(in_param.maxcls);
	a8 = CallocD(in_param.maxcls);
	a9 = CallocD(in_param.maxcls);
	
	temp1 = CallocDD(in_param.maxcls);
	cmu_t = CallocDD(in_param.maxcls);
	csig_t = CallocDD(in_param.maxcls);
	t3 = CallocDD(in_param.maxcls);
	for (j=0; j<in_param.maxcls; j++)
	{
		temp1[j] = CallocD(in_param.conds);
		cmu_t[j] = CallocD(in_param.conds);
		csig_t[j] = CallocD(in_param.conds);
		t3[j] = CallocD(in_param.conds);
	}
		
	temp2 = CallocD(in_param.maxcls);
	
	t4 = CallocD(in_param.maxcls);

	pi_k_t = CallocD(in_param.maxcls);

	/* initialisation */
	for (i=0; i<in_param.maxcls; i++)
	{
		in_param.v[i] = 0.1;
		pi_k[i] = 1.0/in_param.maxcls;
		in_param.calpha[i] = 0.1;
		in_param.cbeta[i] = 0.1;
		
	}

	/* initialisation */
	for (i=0; i<in_param.genes; i++)
	{
		for (j=0; j<in_param.chips; j++)
		{
			exprs[j] = in_param.data_m[j*in_param.genes+i];
			vars[j] = in_param.data_var[j*in_param.genes+i];
		}
		for (j=0; j<in_param.maxcls; j++)
		{
			if (pi_k[j]>0.0)
			{
				alpha_u[i][j] = 0.0;
				beta_u[i][j] = 0.0;
				alpha_eta[i][j] = 0.0;
				beta_eta[i][j] = 0.0;
				for (k=0; k<in_param.conds; k++)
				{
					sig_w[i][j][k] = 0.5;
					mu_w[i][j][k] = 0.0;
				}
				for (k=0; k<in_param.chips; k++)
				{
					mu_t[i][j][k] = exprs[k];
					sig_t[i][j][k] = 1.0; 
					mu_w[i][j][in_param.reps[k]-1] += exprs[k];
				}
				for (k=0; k<in_param.conds; k++)
				{
					mu_w[i][j][k] /= in_param.cond_n[k];
//					Rprintf("mu_w: %f sig_w: %f\n", mu_w[i][j][k], sig_w[i][j][k]);
				}
			}
		}
	}
	
	K_nz = in_param.maxcls;
	F_max = R_NegInf;
	n_hat = 2*in_param.conds+3;
	
	while (K_nz >= in_param.mincls)
	{
		iter = 0;
		F0 = R_NegInf;

		while (1)
		{
			for (i=0; i<in_param.maxcls; i++)
			{
				t2[i] = 0.0;
				for (j=0; j<in_param.conds; j++)
					t1[i][j] = 0.0;
			}

			/* E-step */
			for (i=0; i<in_param.genes; i++)
			{
				for (j=0; j<in_param.chips; j++)
				{
					exprs[j] = in_param.data_m[j*in_param.genes+i];
					vars[j] = in_param.data_var[j*in_param.genes+i];
				}

				Ft0 = R_NegInf;
				iit = 1;

				/* iterations */
				while (1)
				{
					/* u ~ Ga(alpha_u, beta_u);  eta ~ Ga(alpha_eta, beta_eta) */
					for (j=0; j<in_param.maxcls; j++)
					{
						if (pi_k[j]>0.0)
						{
							alpha_u[i][j] = in_param.v[j]/2.0+in_param.conds/2.0;
							alpha_eta[i][j] = in_param.chips/2.0+in_param.calpha[j];
							temp = 0.0;
							for (k=0; k<in_param.conds; k++)
							{
								temp += (mu_w[i][j][k]*mu_w[i][j][k]+sig_w[i][j][k]-2.0*in_param.cmu[j][k]*mu_w[i][j][k]+in_param.cmu[j][k]*in_param.cmu[j][k])/in_param.csig[j][k];
						//		Rprintf("mu_w: %f sig_w: %f cmu: %f csig: %f\n",mu_w[i][j][k],sig_w[i][j][k],in_param.cmu[j][k],in_param.csig[j]);
							}
							beta_u[i][j] = temp/2.0+in_param.v[j]/2.0;
//							if (beta_u[i][j]<0.0) Rprintf("csig: %f\n", in_param.csig[j]);
						//	Rprintf("alpha_u: %f beta_u: %f\n", alpha_u[i][j],beta_u[i][j]);
							temp = 0.0;
							for (k=0; k<in_param.chips; k++)
								temp += (mu_t[i][j][k]*mu_t[i][j][k]+sig_t[i][j][k]+mu_w[i][j][in_param.reps[k]-1]*mu_w[i][j][in_param.reps[k]-1]+sig_w[i][j][in_param.reps[k]-1]-2.0*mu_t[i][j][k]*mu_w[i][j][in_param.reps[k]-1])/2.0;
							beta_eta[i][j] = in_param.cbeta[j]+temp;
							//Rprintf("alpha_eta: %f beta_eta: %f\n", alpha_eta[i][j],beta_eta[i][j]);
						}
					}

					/* A */
					max_A = R_NegInf;
					for (j=0; j<in_param.maxcls; j++)
					{
						if (pi_k[j]>0.0)
						{
							a1[j] = 0.0; a3[j] = 0.0; a6[j] = 0.0; a7[j] = 0.0;
							a2[j] = -lgammafn(in_param.calpha[j])+in_param.calpha[j]*log(in_param.cbeta[j])
								+(in_param.calpha[j]-1.0)*(digamma(alpha_eta[i][j])-log(beta_eta[i][j]))-in_param.cbeta[j]*alpha_eta[i][j]/beta_eta[i][j];
							a4[j] = -lgammafn(in_param.v[j]/2.0)+in_param.v[j]*log(in_param.v[j]/2.0)/2.0
								+(in_param.v[j]/2.0-1.0)*(digamma(alpha_u[i][j])-log(beta_u[i][j]))
								-in_param.v[j]*alpha_u[i][j]/(2.0*beta_u[i][j]);
							a5[j] = -lgammafn(alpha_eta[i][j])+alpha_eta[i][j]*log(beta_eta[i][j])
								+(alpha_eta[i][j]-1.0)*(digamma(alpha_eta[i][j])-log(beta_eta[i][j]))-alpha_eta[i][j];
							a8[j] = -lgammafn(alpha_u[i][j])+alpha_u[i][j]*log(beta_u[i][j])
								+(alpha_u[i][j]-1.0)*(digamma(alpha_u[i][j])-log(beta_u[i][j]))-alpha_u[i][j];
							for (k=0; k<in_param.chips; k++)
							{
							a1[j] += -log(2.0*M_PI)-log(vars[k])/2.0
								-(exprs[k]*exprs[k]-2.0*exprs[k]*mu_t[i][j][k]+mu_t[i][j][k]*mu_t[i][j][k]+sig_t[i][j][k])/(2.0*vars[k])
								+(digamma(alpha_eta[i][j])-log(beta_eta[i][j]))/2.0
								-(mu_t[i][j][k]*mu_t[i][j][k]+sig_t[i][j][k]-2.0*mu_t[i][j][k]*mu_w[i][j][in_param.reps[k]-1]
								+mu_w[i][j][in_param.reps[k]-1]*mu_w[i][j][in_param.reps[k]-1]+sig_w[i][j][in_param.reps[k]-1])
								*alpha_eta[i][j]/(2.0*beta_eta[i][j]);
							a6[j] += -log(2.0*M_PI)/2.0-log(sig_t[i][j][k])/2.0-1.0/2.0;
							}
							for (k=0; k<in_param.conds; k++)
							{
								a3[j] += (digamma(alpha_u[i][j])-log(beta_u[i][j]))/2.0-log(2.0*M_PI)/2.0-log(in_param.csig[j][k])/2.0
									-(mu_w[i][j][k]*mu_w[i][j][k]+sig_w[i][j][k]-2.0*in_param.cmu[j][k]*mu_w[i][j][k]
									+in_param.cmu[j][k]*in_param.cmu[j][k])*alpha_u[i][j]/(2.0*beta_u[i][j]*in_param.csig[j][k]);
								a7[j] += -log(2.0*M_PI)/2.0-log(sig_w[i][j][k])/2.0-1.0/2.0;
							}
							a9[j] = log(pi_k[j]);
						
							A[i][j] = a1[j]+a2[j]+a3[j]+a4[j]-a5[j]-a6[j]-a7[j]-a8[j]+a9[j];
							if (ISNAN(A[i][j]) && in_param.verbose)
							{
								Rprintf("A: %f a1: %f a2: %f a3: %f a4: %f a5: %f a6: %f a7: %f a8: %f a8: %f\n", A[i][j],a1[j],a2[j],a3[j],a4[j],a5[j],a6[j],a7[j],a8[j],a9[j]);
								Rprintf("alpha_u: %f  beta_u: %f csig: %f ", alpha_u[i][j], beta_u[i][j], in_param.csig[j][0]);
								for (k=0; k<in_param.conds; k++)
										Rprintf("cmu: %f ", in_param.cmu[j][k]);
								Rprintf("\n");
								
								return;
							}
							if (A[i][j]>max_A)	max_A = A[i][j];
						}
					}
					
					/* q_kn */
					temp = 0.0;
					Ft = 0.0;
					for (j=0; j<in_param.maxcls; j++) 
						if (pi_k[j]>0.0)
							temp += exp(A[i][j]-max_A);
					for (j=0; j<in_param.maxcls; j++) 
					{
						if (pi_k[j]>0.0)
						{
							in_param.q_kn[i][j] = exp(A[i][j]-max_A-log(temp));
//							if (i==0 && iter==0)
	//					Rprintf("q_kn: %f\n", in_param.q_kn[i][j]);
//						if (in_param.q_kn[i][j]>0)
							Ft += in_param.q_kn[i][j]*(A[i][j]-log(maxD(1.0e-18,in_param.q_kn[i][j])));
						}
					}
					
//					Rprintf("K_nz: %d iter: %d gene: %d Ft: %f\n", K_nz, iter, i, Ft);
					if (Ft-Ft0<absD(in_param.eps*Ft))
					{
//						Rprintf("K_nz: %d iter: %d gene: %d iit: %d Ft: %f\n", K_nz, iter, i, iit, Ft);
						break;
					}
					if (iit>3000 && in_param.verbose) Rprintf("K_nz: %d iter: %d gene: %d iit: %d Ft: %f\n", K_nz, iter, i, iit, Ft);
					Ft0 = Ft;
					
					/* t ~ N(mu_t, sig_t); w ~ N(mu_w, sig_w) */
					for (j=0; j<in_param.maxcls; j++)
						for (k=0; k<in_param.conds; k++)
							temp1[j][k] = 0.0;
					for (j=0; j<in_param.maxcls; j++)
					{
						if (pi_k[j]>0.0)
						{
							for (k=0; k<in_param.conds; k++)	temp1[j][k] = 0.0;
							for (k=0; k<in_param.chips; k++)
							{
								mu_t[i][j][k] = (exprs[k]+(alpha_eta[i][j]/beta_eta[i][j])*vars[k]*mu_w[i][j][in_param.reps[k]-1])
									/(1.0+alpha_eta[i][j]*vars[k]/beta_eta[i][j]);
								sig_t[i][j][k] = vars[k]/(1.0+alpha_eta[i][j]*vars[k]/beta_eta[i][j]);
							
								temp1[j][in_param.reps[k]-1] += mu_t[i][j][k];
							}
							for (k=0; k<in_param.conds; k++)
							{
								mu_w[i][j][k] = (alpha_eta[i][j]*temp1[j][k]*in_param.csig[j][k]/beta_eta[i][j]+alpha_u[i][j]*in_param.cmu[j][k]/beta_u[i][j])
									/(in_param.cond_n[k]*alpha_eta[i][j]*in_param.csig[j][k]/beta_eta[i][j]+alpha_u[i][j]/beta_u[i][j]);
								sig_w[i][j][k] = in_param.csig[j][k]/(in_param.cond_n[k]*alpha_eta[i][j]*in_param.csig[j][k]/beta_eta[i][j]+alpha_u[i][j]/beta_u[i][j]);
							}
						}
					}
					iit++;
					
				} /* while 1 */
				
				for (j=0; j<in_param.maxcls; j++)
				{
					if (pi_k[j]>0.0)
					{
						in_param.un[i][j] = alpha_u[i][j]/beta_u[i][j];
						in_param.logun[i][j] = digamma(alpha_u[i][j])-log(beta_u[i][j]);
						in_param.etan[i][j] = alpha_eta[i][j]/beta_eta[i][j];
						in_param.logetan[i][j] = digamma(alpha_eta[i][j])-log(beta_eta[i][j]);
						for (k=0; k<in_param.conds; k++)
						{
							wn[i][j][k] = mu_w[i][j][k];
							wn2[i][j][k] = mu_w[i][j][k]*mu_w[i][j][k]+sig_w[i][j][k];
							t1[j][k] += in_param.q_kn[i][j]*alpha_u[i][j]*mu_w[i][j][k]/beta_u[i][j];
						}
						t2[j] += in_param.q_kn[i][j]*alpha_u[i][j]/beta_u[i][j];
					}
//					if (i==0)
	//				Rprintf("%f %f %f %f %f %f\n", in_param.un[i][j],in_param.logun[i][j],in_param.etan[i][j],in_param.logetan[i][j],wn[i][j][0],wn2[i][j][0]);
				}
				
			} /* for genes */
			
			K_nz = 0;
			temp = 0.0;
			
			for (j=0; j<in_param.maxcls; j++)
			{
				if (pi_k[j]>0.0)
				{
					K_nz++;
					temp += log(in_param.genes*pi_k[j]/12.0);
				}
			}
				
			F = 0.0;
			for (i=0; i<in_param.genes; i++)
				for (j=0; j<in_param.maxcls; j++)
					if (pi_k[j]>0.0)
						F += in_param.q_kn[i][j]*(A[i][j]-log(maxD(1.0e-18,in_param.q_kn[i][j])));
				F = F-n_hat*temp/2.0-K_nz*log(in_param.genes/12.0)/2.0-K_nz*(n_hat+1)/2.0;
			
			iter++;
            if(in_param.verbose)
    			Rprintf("K_nz: %d iter: %d F: %f\n", K_nz, iter, F);
			
			if (F-F0<absD(in_param.eps*F) || K_nz<in_param.mincls)
			{
				break;
			}
			F0 = F;
			F_t = F;
			K_nz_t = K_nz;
			for (j=0; j<in_param.maxcls; j++)
			{
				pi_k_t[j] = pi_k[j];

				for (k=0; k<in_param.conds; k++)	
				{
					cmu_t[j][k] = in_param.cmu[j][k];
					csig_t[j][k] = in_param.csig[j][k];
				}
				for (k=0; k<in_param.genes; k++)
					q_kn_t[k][j] = in_param.q_kn[k][j];
			}

			/* M-step */
			temp3 = 0.0;
			for (k=0; k<in_param.maxcls; k++)
			{
				if (pi_k[k]>0.0)
				{
					temp = 0.0;
					t4[k] = 0.0;
					for (j=0; j<in_param.conds; j++)
					{
						in_param.cmu[k][j] = t1[k][j]/t2[k];
						t3[k][j] = 0.0;
					}

					for (i=0; i<in_param.genes; i++) 
					{
//						temp4 = 0.0;
						temp += in_param.q_kn[i][k];
						for (j=0; j<in_param.conds; j++)
						{
							t3[k][j] += in_param.q_kn[i][k]*in_param.un[i][k]*(wn2[i][k][j]
								-2.0*in_param.cmu[k][j]*wn[i][k][j]+in_param.cmu[k][j]*in_param.cmu[k][j]);
//							if (in_param.un[i][j]*(wn2[i][k][j]-2.0*in_param.cmu[k][j]*wn[i][k][j]
	//							+in_param.cmu[k][j]*in_param.cmu[k][j])<0)
//								Rprintf("%d mu_w: %f sig_w: %f cmu: %f csig: %f 1/eta: %f\n", i, mu_w[i][k][j],sig_w[i][k][j], in_param.cmu[k][j], in_param.csig[k]/in_param.un[i][k], beta_eta[i][k]/alpha_eta[i][k]);
						}
//						t3[k] += in_param.q_kn[i][k]*temp4;
						t4[k] += in_param.q_kn[i][k];
					}
					pi_k[k] = maxD(0.0,temp-n_hat/2.0);
					temp3 += pi_k[k];
				}
			}
			for (k=0; k<in_param.maxcls; k++)	
			{
				if (pi_k[k]>0.0)
				{
					pi_k[k] = pi_k[k]/temp3;
			//		Rprintf("pi_k: %f ",pi_k[k]);
					for (j=0; j<in_param.conds; j++)
					{
						in_param.csig[k][j] = t3[k][j]/t4[k];
						if (in_param.csig[k][j]<0 && in_param.verbose)
							Rprintf("t3: %f t4: %f\n", t3[k][j], t4[k]);
					}
				}
			}
	//		Rprintf("\n");			
			for (in_param.i_opt=0; in_param.i_opt<in_param.maxcls ; in_param.i_opt++)
			{
				if (pi_k[in_param.i_opt]>0.0)
				{
					in_param.var_flag = 1;
					donlp2();
					//Rprintf("%d  %d %d %f %d\n", iter, icf, icgf, fx, (int)optite+11);					
					in_param.var_flag = 2;
					//Rprintf("%d  %d %d %f %d\n", iter, icf, icgf, fx, (int)optite+11);	
					donlp2();
				}
			}
			
		} /* while 1 */
		
		if (F_t>=F_max)
		{
			F_max = F_t;
			in_param.K_best = K_nz_t;
			in_param.F_value = F_t;
			for (j=0; j<in_param.maxcls; j++)
			{
				in_param.pi_k_best[j] = pi_k_t[j];

				for (k=0; k<in_param.conds; k++)	
				{
					in_param.cmu_best[j][k] = cmu_t[j][k];
					in_param.csig_best[j][k] = csig_t[j][k];
				}
				for (k=0; k<in_param.genes; k++)
					in_param.q_kn_best[k][j] = q_kn_t[k][j];
			}
		}
		
		temp = 2.0;
		i = -1;
		for (j=0; j<in_param.maxcls; j++)
		{
			if (pi_k[j]<temp && pi_k[j]>0.0)
			{
				temp = pi_k[j];
				i = j;
			}
		}
		pi_k[i] = 0.0;
		
		temp3 = 0.0;
		for (j=1; j<in_param.maxcls; j++)
			temp3 += pi_k[j];
		for (j=1; j<in_param.maxcls; j++)
			pi_k[j] /= temp3;
		K_nz--;
		
	} /* while K_nz>=K_min */
		

    /*fprintf(pf,"%f\n",sqrt(1.0/in_param.lamda_m));*/
	
    /*************** Free memory ******************************/
/*	FreeD(exprs);
	FreeD(vars);
	FreeD(pi_k);
	FreeD(t2);

	for (i=0; i<in_param.genes; i++)
	{
		for (j=0; j<in_param.maxcls; j++)
		{
			FreeD(wn[i][j]);
			FreeD(wn2[i][j]);
			FreeD(mu_w[i][j]);
			FreeD(sig_w[i][j]);
			FreeD(mu_t[i][j]);
			FreeD(sig_t[i][j]);
		}
		FreeD(A[i]);
		FreeDD(wn[i]);
		FreeDD(wn2[i]);
		FreeDD(mu_w[i]);
		FreeDD(sig_w[i]);
		FreeDD(mu_t[i]);
		FreeDD(sig_t[i]);
		
		FreeD(alpha_u[i]);
		FreeD(beta_u[i]);
		FreeD(alpha_eta[i]);
		FreeD(beta_eta[i]);
		
		FreeD(q_kn_t[i]);
	}	
	FreeDD(A);
	FreeDDD(wn);
	FreeDDD(wn2);
	FreeDDD(mu_w);
	FreeDDD(sig_w);
	FreeDDD(mu_t);
	FreeDDD(sig_t);
	FreeDD(alpha_u);
	FreeDD(beta_u);
	FreeDD(alpha_eta);
	FreeDD(beta_eta);
	
	FreeDD(q_kn_t);

	for (i=0; i<in_param.maxcls; i++)
	{
		FreeD(t1[i]);
		FreeD(cmu_t[i]);
	}
	FreeDD(t1);
	FreeDD(cmu_t);
	
	FreeD(a1);
	FreeD(a2);
	FreeD(a3);
	FreeD(a4);
	FreeD(a5);
	FreeD(a6);
	FreeD(a7);
	FreeD(a8);
	FreeD(a9);

	for (j=0; j<in_param.maxcls; j++)
	{
		FreeD(temp1[j]);
		FreeD(t3[j]);
		FreeD(csig_t[j]);
	}
	
	FreeDD(temp1);
	FreeD(temp2);
	
	FreeDD(t3);
	FreeD(t4);
	
	FreeDD(csig_t);
	FreeD(pi_k_t);*/
	/***********************************************************/
}

SEXP pumaclustii_c(SEXP Mmat, SEXP Stdmat, SEXP conds, SEXP reps, SEXP mincls, SEXP maxcls, SEXP centers, SEXP clsig, SEXP verbose, SEXP eps, SEXP del0)
{
	SEXP dim=NULL;	
	SEXP res=NULL, CIndex=NULL, Centers=NULL, ClusterSig=NULL, pjPerData=NULL, BestK=NULL, Fvalue=NULL;
	int *outCIndex=NULL;    
	double *outClusterCenters=NULL;
	double *outClusterSig=NULL, maxv;
	double *outLikeliPerGene=NULL, *cmu=NULL, *csig=NULL;
	int *outBestK=NULL;
	double *outFvalue=NULL;
	int i,j,k;

	initialparams_pumaclustii();

	PROTECT(dim = getAttrib(Mmat, R_DimSymbol));
	in_param.genes = INTEGER(dim)[0];
	in_param.chips = INTEGER(dim)[1];
	
	in_param.conds = INTEGER_POINTER(AS_INTEGER(conds))[0];
	in_param.mincls = INTEGER_POINTER(AS_INTEGER(mincls))[0];
	in_param.maxcls = INTEGER_POINTER(AS_INTEGER(maxcls))[0];
	in_param.reps = INTEGER_POINTER(AS_INTEGER(reps));
	cmu = NUMERIC_POINTER(AS_NUMERIC(centers));
	csig = NUMERIC_POINTER(AS_NUMERIC(clsig));
    in_param.verbose = INTEGER_POINTER(AS_INTEGER(verbose))[0];
	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];
	in_param.del0 = NUMERIC_POINTER(AS_NUMERIC(del0))[0];
	in_param.data_m = NUMERIC_POINTER(AS_NUMERIC(Mmat));
	in_param.data_var = NUMERIC_POINTER(AS_NUMERIC(Stdmat));

	allocatemem_pumaclustii();

	for (i=0; i<in_param.maxcls; i++)
	{
		for (j=0; j<in_param.conds; j++)
		{
			in_param.cmu[i][j] = cmu[j*in_param.maxcls+i];
			in_param.csig[i][j] = csig[j*in_param.maxcls+i];
		}
	}

	for (i=0; i<in_param.conds; i++)
		in_param.cond_n[i] = 0;
	for (i=0; i<in_param.chips; i++)
	{
		in_param.cond_n[in_param.reps[i]-1]++;
	}

	Rprintf("Clustering is performing ...\n");
	workout_pumaclustii();  
	Rprintf("Done.\n");

	PROTECT(CIndex = allocVector(INTSXP, in_param.genes));
	PROTECT(Centers = allocMatrix(REALSXP, in_param.K_best, in_param.conds));
	PROTECT(ClusterSig = allocMatrix(REALSXP, in_param.K_best, in_param.conds));
	PROTECT(pjPerData = allocMatrix(REALSXP, in_param.genes, in_param.K_best));
	PROTECT(BestK = allocVector(INTSXP,1));
	PROTECT(Fvalue = allocVector(REALSXP,1));
	PROTECT(res = allocVector(VECSXP,6));	
	outCIndex = INTEGER_POINTER(AS_INTEGER(CIndex));
	outClusterCenters = NUMERIC_POINTER(AS_NUMERIC(Centers));
	outClusterSig = NUMERIC_POINTER(AS_NUMERIC(ClusterSig));
	outLikeliPerGene = NUMERIC_POINTER(AS_NUMERIC(pjPerData));
	outBestK = INTEGER_POINTER(AS_INTEGER(BestK));
	outFvalue = NUMERIC_POINTER(AS_NUMERIC(Fvalue));
	
	*(outBestK) = in_param.K_best;
	*(outFvalue) = in_param.F_value;
	for (i=0; i<in_param.genes; i++)
	{
		outCIndex[i] = 1;
		if (in_param.pi_k_best[0]>0.0)
			j = 1;
		else
			j = 0;
		maxv = in_param.q_kn_best[i][0];
		for (k=1; k<in_param.maxcls; k++)
			if (in_param.pi_k_best[k]>0.0)
			{
				j++;
				if (in_param.q_kn_best[i][k]>maxv)
				{
					outCIndex[i] = j;
					maxv = in_param.q_kn_best[i][k];
				}
			}
	}	

	for (j=0; j<in_param.conds; j++)
	{
		i = 0;
		for (k=0; k<in_param.maxcls; k++)
		{
			if (in_param.pi_k_best[k]>0.0)
			{
				outClusterCenters[j*in_param.K_best+i] = in_param.cmu_best[k][j];
				outClusterSig[j*in_param.K_best+i] = in_param.csig_best[k][j];
				i++;
			}
		}
	} 

	for (i=0; i<in_param.genes; i++)
	{
		j = 0;
		for (k=0; k<in_param.maxcls; k++)
			if (in_param.pi_k_best[k]>0.0)
			{
				outLikeliPerGene[j*in_param.genes+i] = in_param.q_kn_best[i][k];
				j++;
			}
	}	


	SET_VECTOR_ELT(res, 0, CIndex);
	SET_VECTOR_ELT(res, 1, Centers);
	SET_VECTOR_ELT(res, 2, ClusterSig);
	SET_VECTOR_ELT(res, 3, pjPerData);
	SET_VECTOR_ELT(res, 4, BestK);
	SET_VECTOR_ELT(res, 5, Fvalue);
	
	/*freemem_pumaclustii();*/
	UNPROTECT(8);
	
	return res;
}

/* **************************************************************************** */
/*                              donlp2-intv size initialization                           */
/* **************************************************************************** */
void user_init_size_pumaclustii(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X


    /* problem dimension n = dim(donlp2_x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
      
    if (in_param.var_flag == 1)
    	n = 1;
    else
    	n = 2;
    	
    nstep = 20;
    
    nlin   =  0;
    nonlin =  0;
    iterma = 4000;
}

/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_pumaclustii(void) {
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
	if (in_param.var_flag == 1)
	{
		donlp2_x[1] = in_param.v[in_param.i_opt];
		
		low[1] = LOWBOUND;
		up[1] = big;
    }
    else
    {
    	donlp2_x[1] = in_param.calpha[in_param.i_opt];
    	donlp2_x[2] = in_param.cbeta[in_param.i_opt];
 		low[1] = LOWBOUND;
		up[1] = big;
		low[2] = LOWBOUND;
		up[2] = big;
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
void setup_pumaclustii(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
   
    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_pumaclustii(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"

	if (in_param.var_flag == 1)
	{
		in_param.v[in_param.i_opt] = donlp2_x[1];
//		Rprintf("v: %f\n", donlp2_x[1]);
	}
	else
	{
		in_param.calpha[in_param.i_opt] = donlp2_x[1];
		in_param.cbeta[in_param.i_opt] = donlp2_x[2];
//		Rprintf("alpha: %f beta: %f\n", donlp2_x[1], donlp2_x[2]);
	}
    return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_pumaclustii(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	int i, j;
	double t1, t2, t3;
	
	j = in_param.i_opt;
	*fx = 0.0;
	
	if (in_param.var_flag == 1 )
	{
		for (i=0; i<in_param.genes; i++)
		{
			*fx += in_param.q_kn[i][j]*(-lgammafn(donlp2_x[1]/2.0)+donlp2_x[1]*log(donlp2_x[1]/2.0)/2.0
				+donlp2_x[1]*in_param.logun[i][j]/2.0-donlp2_x[1]*in_param.un[i][j]/2.0);
		}
	}
	else
	{
		t1 = 0.0;
		t2 = 0.0;
		t3 = 0.0;
		for (i=0; i<in_param.genes; i++)
		{
			t1 += in_param.q_kn[i][j];
			t2 += in_param.q_kn[i][j]*in_param.logetan[i][j];
			t3 += in_param.q_kn[i][j]*in_param.etan[i][j];
		}
		*fx = t1*(-lgammafn(donlp2_x[1])+donlp2_x[1]*log(donlp2_x[2]))+donlp2_x[1]*t2-donlp2_x[2]*t3;
	}
	
	*fx = -*fx;

    return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_pumaclustii(DDOUBLE donlp2_x[],DDOUBLE gradf[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	int i, j;
	double t1, t2, t3;
	
	j = in_param.i_opt;
	if (in_param.var_flag == 1)
	{
		gradf[1] = 0.0;
		for (i=0; i<in_param.genes; i++)
		{
			gradf[1] += in_param.q_kn[i][j]*(-digamma(donlp2_x[1]/2.0)/2.0+log(donlp2_x[1]/2.0)/2.0+1.0/2.0-
				in_param.un[i][j]/2.0+in_param.logun[i][j]/2.0);
		}
		gradf[1] = -gradf[1];
	}
	else
	{
		t1 = 0.0;
		t2 = 0.0;
		t3 = 0.0;
		for (i=0; i<in_param.genes; i++)
		{
			t1 += in_param.q_kn[i][j];
			t2 += in_param.q_kn[i][j]*in_param.logetan[i][j];
			t3 += in_param.q_kn[i][j]*in_param.etan[i][j];
		}
		gradf[1] = t1*(-digamma(donlp2_x[1])+log(donlp2_x[2]))+t2;
		gradf[2] = t1*donlp2_x[1]/donlp2_x[2]-t3;
		gradf[1] = -gradf[1];
		gradf[2] = -gradf[2];
	}
    return;
}

/* **************************************************************************** */
/*  no nonlinear constraints */
/* **************************************************************************** */
void econ_pumaclustii(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[], 
              LLOGICAL err[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_pumaclustii(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_pumaclustii(IINTEGER mode) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}
