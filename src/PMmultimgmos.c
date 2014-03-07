/* **************************************************************************** */
/*        PMmulti-mgMOS                                                         */
/* **************************************************************************** */
#include <math.h>
#include "o8para.h"
#include "PMmultimgmoshead.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

static pmexpparam in_param;
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

void econ_mmgmospm(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_mmgmospm(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_mmgmospm(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_mmgmospm(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_mmgmospm(IINTEGER mode);
void freemem_mmgmospm();
void initialparams_mmgmospm();
void setup_mmgmospm();
void solchk_mmgmospm();
void user_init_mmgmospm(void);
void user_init_size_mmgmospm(void);
void allocatemem_mmgmospm();

/*econ
econgrad
ef
egradf
eval_extern
freemem
initialparams
setup
solchk
user_init
user_init_size
allocatemem
*/

double pmdierfc(double x)  /* inverse error function */
{
        return -qnorm((1-x)/2,0,1,TRUE,FALSE)/sqrt(2);
}

double pmerfc(double x)    /* complementary error function */
{
        return 2*pnorm(sqrt(2)*x,0,1,FALSE,FALSE);
}

void initialparams_mmgmospm()
{
	in_param.data_pm = NULL;
	in_param.parameters = NULL;
	in_param.prctiles = NULL;
	econ = econ_mmgmospm;
	econgrad = econgrad_mmgmospm;
	ef = ef_mmgmospm;
	egradf = egradf_mmgmospm;
	eval_extern = eval_extern_mmgmospm;
	freemem = freemem_mmgmospm;
	initialparams = initialparams_mmgmospm;
	setup = setup_mmgmospm;
	solchk = solchk_mmgmospm;
	user_init = user_init_mmgmospm;
	user_init_size = user_init_size_mmgmospm;
	allocatemem = allocatemem_mmgmospm;
}
void pmgetgenedata(int g)
{
	int i,j;
	
	in_param.num_probe = in_param.probesets[g];
 
	
	for (i=0; i<in_param.num_probe; i++)
	{
		in_param.totalprobe ++;
		for (j=0; j<in_param.chips; j++)
		{
			in_param.pm[i][j] = in_param.data_pm[j*in_param.probes+in_param.totalprobe];
			
		}
	}
}

void allocatemem_mmgmospm()
{
	int i;
	
	
	in_param.probesets = (int*)R_alloc(in_param.genes, sizeof(int));
	
	in_param.parameters = (double**)R_alloc(in_param.genes, sizeof(double*));
	for (i=0; i<in_param.genes; i++)
	{
		in_param.parameters[i] = (double*)R_alloc((in_param.conds+2), sizeof(double));  
		in_param.probesets[i] = 0;
	}
		
	
}

  
void freemem_mmgmospm()
{
	int j;
	
	for (j=0; j<in_param.genes; j++)
		if (in_param.parameters[j] != NULL) Free(in_param.parameters[j]);
        
	if (in_param.parameters != NULL) Free(in_param.parameters);

	if (in_param.probesets != NULL) Free(in_param.probesets);
}


void pmcalparameters()
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
	int niter = 1, nx;
	double fstart;
	int p, i, j;
	int finishflag = 0;
     
  	fstart = HUGE_VAL;
	
	nx = in_param.conds+2;
	
	
	while (1)
	{
		void R_CheckUserInterrupt(void);
		
		in_param.totalprobe = -1;
		for (p=0; p<in_param.genes; p++)   //in_param.genes
		{
			void R_CheckUserInterrupt(void);
			
			/* do optimisation for each gene*/
			in_param.cur_gene = p;
			
			pmgetgenedata(p);
			
			in_param.flag = 0;
			
			if (in_param.num_probe > 1)
				donlp2();
			
/*			in_param.parameters[p][nx] = fx;*/
/*			in_param.parameters[p][nx+1] = optite+11;*/

/*			Rprintf("%d_mmgmos %d_mmgmos  %d_mmgmos %d_mmgmos %d_mmgmos\n", niter, p+1, icf, icgf, (int)optite+11);*/
			if ((int)p%500 == 0)
				Rprintf(".");
  		}
	
		
		
	
		/* check the termination criteria*/
		if (fstart-fx < in_param.eps*fx || finishflag == 1) 
		{
			if (in_param.saveflag == TRUE)
			{
				FILE *df = fopen("par_pmmmgmos.txt", "wt");
		
				if (!df)
				{ 
					Rprintf("Cannot open file for saving parameters\n");
					break;
				}
		
				for (i=0; i<in_param.genes; i++)
				{
					for (j=0; j<in_param.conds+2; j++)
						fprintf(df," %f", in_param.parameters[i][j]);
					fprintf(df,"\n");
				}
				fclose(df);
				
				
			}
			Rprintf("\n");

			break;
		}
		
		
		niter++;
		fstart = fx;
	}
}

void pmcalexpression()
{
	int p, i, j, k, q;
	double  alphai, c, d;
	double mu_Gauss, var_Gauss, mu_truncGauss, var_truncGauss, kk;
	double secdrv, firstdrv,temp;
	double t1, ym, aym, yym;
	
	//FILE *pf=NULL;
	
	/*pf=fopen("result.txt", "wt");*/
	
	in_param.totalprobe = -1;
	for (p=0; p<in_param.genes; p++)
	{
		c = in_param.parameters[p][in_param.conds];
		d = in_param.parameters[p][in_param.conds+1];
				
		pmgetgenedata(p);
		
		if (in_param.num_probe > 1)
		{
			for (i=0; i<in_param.conds; i++)
			{
				alphai = in_param.parameters[p][i];
				
				t1 = 0.0;
				for (j=0; j<in_param.conds; j++)
					t1 += in_param.parameters[p][j];
				t1 += c;	

				secdrv = in_param.num_probe*(trigamma(t1)-trigamma(alphai));
				var_Gauss = -1.0/secdrv;

				if (alphai<1.0e-6)
				{
					yym = 0.0;
					aym = 0.0;
					for (k=0; k<in_param.num_probe; k++)
					{
						ym = 0.0;
						for (q=0; q<in_param.conds; q++)
							ym += in_param.pm[k][q];
						ym += d;
						yym += log(ym);
						aym += log(in_param.pm[k][i]);
					}

					firstdrv = in_param.num_probe*(digamma(t1)-digamma(alphai))-yym+aym;
					mu_Gauss = var_Gauss*firstdrv;
				}
				else
					mu_Gauss = alphai;

				kk=2.0/pmerfc(-mu_Gauss/sqrt(2.0*var_Gauss));

				mu_truncGauss = kk*(sqrt(var_Gauss)*exp(-mu_Gauss*mu_Gauss/(2.0*var_Gauss))/sqrt(2.0*M_PI)
						+mu_Gauss*pmerfc(-mu_Gauss/sqrt(2.0*var_Gauss))/2.0);
				var_truncGauss = kk*((var_Gauss+(mu_Gauss-mu_truncGauss)*(mu_Gauss-mu_truncGauss))*pmerfc(-mu_Gauss/sqrt(2.0*var_Gauss))/2.0
						+sqrt(var_Gauss/(2.0*M_PI))*exp(-mu_Gauss*mu_Gauss/(2.0*var_Gauss))*(mu_Gauss-2.0*mu_truncGauss));

				/* calculate expression value -- the approximate mean of <log(s)> */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+p] = (digamma(mu_truncGauss)+log(d)-digamma(c))/log(2.0)+
										tetragamma(mu_truncGauss)*var_truncGauss/(2.0*log(2.0)*log(2.0));

				/* calculate standard deviation */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+in_param.genes+p] = sqrt(pow(trigamma(mu_truncGauss),2)*var_truncGauss/(log(2.0)*log(2.0)));
                                /* calculate percentiles */
				for (j=0; j<in_param.num_prctile; j++)
				{
					temp = mu_Gauss+sqrt(2.0*var_Gauss)*pmdierfc(1-2.0*(1.0-in_param.prctiles[j])/kk);
					in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+(j+2)*in_param.genes+p] = (digamma(temp)+log(d)-digamma(c))/log(2.0);
				}

				
			}	
				
		if  ((int)p%500 == 0)
			Rprintf(".");
		}
	} /* end of for genes */
	Rprintf("\n");
}


SEXP pmmmgmos_c(SEXP PMmat, SEXP ngenes, SEXP probeNames,SEXP prctiles, SEXP nprc,  SEXP saveflag, SEXP eps)
{
	void R_CheckUserInterrupt(void);
         SEXP dim=NULL;	
	SEXP res=NULL;
	
	int i, j;
	const char *geneName=NULL;

/*	Rprintf("initialparams_mmgmospm ");
*/	initialparams_mmgmospm();
	
	PROTECT(dim = getAttrib(PMmat, R_DimSymbol));
	
	in_param.chips = INTEGER(dim)[1];
	in_param.conds = in_param.chips;
	in_param.num_prctile = INTEGER(nprc)[0];
	in_param.genes = INTEGER(ngenes)[0];
	//in_param.genes = 10;//
	in_param.probes = INTEGER(dim)[0];

	//in_param.probes = 160;//	
	in_param.data_pm = NUMERIC_POINTER(AS_NUMERIC(PMmat));
	in_param.prctiles = NUMERIC_POINTER(AS_NUMERIC(prctiles));
	
	
	in_param.saveflag = LOGICAL_POINTER(AS_LOGICAL(saveflag))[0];
        
	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];
	 allocatemem_mmgmospm();
	
	
	
	geneName = CHAR(STRING_ELT(probeNames,0));
	j = 0;
	for (i=0; i<in_param.probes; i++)
	{
		if (!strcmp(geneName, CHAR(STRING_ELT(probeNames,i))))
			in_param.probesets[j]++;
		else
		{
			geneName = CHAR(STRING_ELT(probeNames,i));
			j++;
			in_param.probesets[j] = 1;
		}
	}
	
	Rprintf("Model optimising ");
 // Rprintf("%d\n",in_param.chips);
        // Rprintf("%d\n",in_param.genes);//
	pmcalparameters();

	//PROTECT(res = allocMatrix(REALSXP, in_param.genes*2, in_param.conds));
          PROTECT(res = allocMatrix(REALSXP, in_param.genes*(2+in_param.num_prctile), in_param.conds));
	in_param.outp = NUMERIC_POINTER(AS_NUMERIC(res));
	
	Rprintf("Expression values calculating ");
	pmcalexpression();
	//freemem_mmgmospm();
	Rprintf("Done.\n");
	UNPROTECT(2);
         return res;


}



/* **************************************************************************** */
/*                              donlp2-intv size initialization                           */
/* **************************************************************************** */
void user_init_size_mmgmospm(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X

  
    /* problem dimension n = dim(donlp2_x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
       
    if (in_param.flag == 0)
    {
    	n = in_param.chips+2;
    	nstep = 20;
    }
   /* else if (in_param.flag == 1)
    {
    	n = 1;
        nstep = 40;
    }
    else
    {
    	n = 4;
	nstep = 20;
    }*/
    
    nlin   =  0;
    nonlin =  0;
    iterma = 4000;
  
}

/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_mmgmospm(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X
    
    static IINTEGER i,j;

    silent = TRUE;
/*    intakt = TRUE;*/
    
    big = 1.e20;
   
    if (in_param.flag==0)
    {
	
	/* initialise the parameters */
	for (i=1; i<=in_param.conds; i++)
	{
		donlp2_x[i] = 2.0;
		low[i] = ALOW;
		up[i] = big;
	}
	/*for (j=1; j<=2; j++)
	{
		donlp2_x[in_param.conds+j] = 10.0;
		low[in_param.conds+j] = ALOW;
		up[in_param.conds+j] = big;
	}*/
	donlp2_x[in_param.conds+1] = 10.0;
	low[in_param.conds+1] = CLOW;
	up[in_param.conds+1] = big;
	donlp2_x[in_param.conds+2] = 10.0;
	low[in_param.conds+2] = DLOW;
	up[in_param.conds+2] = big;
    }
   
   
    analyt = TRUE;
    epsdif = 1.e-16;  
    
    nreset = n;
    
    del0 = 1.0e0;
    tau0 = 1.0e1;
    tau  = 0.1e0;
    
    return;
}

/* **************************************************************************** */
/*                                 special setup                                */
/* **************************************************************************** */
void setup_mmgmospm(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
   
    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_mmgmospm(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"

    IINTEGER i;
    
    if (in_param.flag == 0)
    {
    	for (i=0; i<n; i++)
		in_param.parameters[in_param.cur_gene][i] = donlp2_x[i+1];
    }
   

    return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_mmgmospm(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	double alphaii[MAX_NUM_COND]={0.0}, c, d, t11, t22;
	double ym[MAX_NUM_PROBE_PM]={0.0}, aym[MAX_NUM_PROBE_PM]={0.0};
	
	int  j, k;
	
	*fx = 0.0;
	
	

	if (in_param.flag==0)
	{
		
		c = donlp2_x[in_param.conds+1];
		d = donlp2_x[in_param.conds+2];
	
		t11 = 0.0;
		t22 = 0.0;
		for (j=0; j<in_param.conds; j++)
		{
			alphaii[j] = donlp2_x[j+1];
			
			
			t11 +=alphaii[j];
			t22 += lgammafn(alphaii[j]);
		}
		t11+= c;
		
		for (k=0; k<in_param.num_probe; k++)
		{
			for (j=0; j<in_param.conds; j++)
			{
				ym[k] += in_param.pm[k][j];
				aym[k] += (alphaii[j]-1.0)*log(in_param.pm[k][j]);
			}
			ym[k] += d;
		
			*fx += c*log(d) + lgammafn(t11) - lgammafn(c) - t22 - t11*log(ym[k]) + aym[k];
		}
	
		*fx = -*fx;
		
	}
	
    return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_mmgmospm(DDOUBLE donlp2_x[],DDOUBLE gradf[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	double alphaii[MAX_NUM_COND]={0.0}, c, d, t1;
	double ym[MAX_NUM_PROBE_PM]={0.0};
	int  j, k;
    
	if (in_param.flag==0)
	{
		for (j=0; j<in_param.conds+2; j++)
			gradf[j+1] = 0.0;
		
		c = donlp2_x[in_param.conds+1];
		d = donlp2_x[in_param.conds+2];
	
		t1 = 0.0;
		for (j=0; j<in_param.conds; j++)
		{
			alphaii[j] = donlp2_x[j+1];
					
			t1 += alphaii[j];
		}
		t1 += c;
		
		for (k=0; k<in_param.num_probe; k++)
		{
			for (j=0; j<in_param.conds; j++)
			{
					
				ym[k] += in_param.pm[k][j];
				
			}
			ym[k] += d;
		
			for (j=0; j<in_param.conds; j++)
			{
				gradf[j+1] += digamma(t1)+log(in_param.pm[k][j])-log(ym[k])-digamma(alphaii[j]);
				
			}
			gradf[in_param.conds+1] += log(d)+digamma(t1)-digamma(c)-log(ym[k]);
			gradf[in_param.conds+2] += c/d-t1/ym[k];
		}
	
		for (j=0; j<in_param.conds+2; j++)
			gradf[j+1] = -gradf[j+1];
	}
	
		
		
	
    return;
}

/* **************************************************************************** */
/*  no nonlinear constraints */
/* **************************************************************************** */
void econ_mmgmospm(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[], 
              LLOGICAL err[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_mmgmospm(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_mmgmospm(IINTEGER mode) {
    #define  X extern
    #include "o8comm.h"


    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}

