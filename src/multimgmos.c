/* **************************************************************************** */
/*        multi-mgMOS C implementation                                          */
/* **************************************************************************** */
#include <math.h>
#include "o8para.h"
#include "multimgmos.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

static expparam in_param;
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

void econ_mmgmos(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_mmgmos(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_mmgmos(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_mmgmos(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_mmgmos(IINTEGER mode);
void freemem_mmgmos();
void initialparams_mmgmos();
void setup_mmgmos();
void solchk_mmgmos();
void user_init_mmgmos(void);
void user_init_size_mmgmos(void);
void allocatemem_mmgmos();

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

double dierfc(double x)  /* inverse error function */
{
        return -qnorm((1-x)/2,0,1,TRUE,FALSE)/sqrt(2);
}

double erfc(double x)    /* complementary error function */
{
        return 2*pnorm(sqrt(2)*x,0,1,FALSE,FALSE);
}

void initialparams_mmgmos()
{
	in_param.data_pm = NULL;
	in_param.data_mm = NULL;
	in_param.probesets = NULL;
	in_param.replicates = NULL;
	in_param.paramphi = NULL;
	in_param.prctiles = NULL;
	in_param.parameters = NULL;
	
	econ = econ_mmgmos;
	econgrad = econgrad_mmgmos;
	ef = ef_mmgmos;
	egradf = egradf_mmgmos;
	eval_extern = eval_extern_mmgmos;
	freemem = freemem_mmgmos;
	initialparams = initialparams_mmgmos;
	setup = setup_mmgmos;
	solchk = solchk_mmgmos;
	user_init = user_init_mmgmos;
	user_init_size = user_init_size_mmgmos;
	allocatemem = allocatemem_mmgmos;
}
void getgenedata(int g)
{
	int i,j;
	
	in_param.num_probe = in_param.probesets[g];
	
	for (i=0; i<in_param.num_probe; i++)
	{
		in_param.totalprobe ++;
		for (j=0; j<in_param.chips; j++)
		{
			in_param.pm[i][j] = in_param.data_pm[j*in_param.probes+in_param.totalprobe];
			in_param.mm[i][j] = in_param.data_mm[j*in_param.probes+in_param.totalprobe];
		}
	}
}

void allocatemem_mmgmos()
{
	int i;
	
	in_param.replicates = (int*)R_alloc(in_param.conds, sizeof(int));
	in_param.probesets = (int*)R_alloc(in_param.genes, sizeof(int));
	
	in_param.parameters = (double**)R_alloc(in_param.genes, sizeof(double*));
	for (i=0; i<in_param.genes; i++)
	{
		in_param.parameters[i] = (double*)R_alloc((in_param.conds+in_param.chips+2), sizeof(double));  
		in_param.probesets[i] = 0;
	}
		
	in_param.paramphi = (double*)R_alloc(1, sizeof(double));
}

  
void freemem_mmgmos()
{
	int j;
	
	for (j=0; j<in_param.genes; j++)
		if (in_param.parameters[j] != NULL) Free(in_param.parameters[j]);
	if (in_param.parameters != NULL) Free(in_param.parameters);

	if (in_param.paramphi != NULL) Free(in_param.paramphi);
	
	
	if (in_param.replicates != NULL) Free(in_param.replicates);
	if (in_param.probesets != NULL) Free(in_param.probesets);
}


void calparameters()
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
	
	nx = in_param.conds+in_param.chips+2;
	
	
	while (1)
	{
		void R_CheckUserInterrupt(void);
		
		in_param.totalprobe = -1;
		for (p=0; p<in_param.genes; p++)
		{
			void R_CheckUserInterrupt(void);
			
			/* do optimisation for each gene*/
			in_param.cur_gene = p;
			
			getgenedata(p);
			
			in_param.flag = 0;
			
			if (in_param.num_probe > 1)
				donlp2();
			
/*			in_param.parameters[p][nx] = fx;*/
/*			in_param.parameters[p][nx+1] = optite+11;*/

/*			Rprintf("%d_mmgmos %d_mmgmos  %d_mmgmos %d_mmgmos %d_mmgmos\n", niter, p+1, icf, icgf, (int)optite+11);*/
			if ((int)p%500 == 0)
				Rprintf(".");
  		}
	
		if (in_param.phi == 0)
			finishflag = 1;
		else
		{
			in_param.flag = 1;
			
			/* do optimisation to estimate phi*/
			donlp2();
		}
	
		/* check the termination criteria*/
		if (fstart-fx < in_param.eps*fx || finishflag == 1) 
		{
			if (in_param.saveflag == TRUE)
			{
				FILE *df = fopen("par_mmgmos.txt", "wt");
		
				if (!df)
				{ 
					Rprintf("Cannot open file for saving parameters\n");
					break;
				}
		
				for (i=0; i<in_param.genes; i++)
				{
					for (j=0; j<in_param.conds+in_param.chips+2; j++)
						fprintf(df," %f", in_param.parameters[i][j]);
					fprintf(df,"\n");
				}
				fclose(df);
				
				df = fopen("phi_mmgmos.txt", "wt");
				if (!df)
				{ 
					Rprintf("Cannot open file for saving phi\n");
					break;
				}
				fprintf(df, "%f\n", in_param.phi);
				fclose(df);
			}
			Rprintf("\n");

			break;
		}
		
		in_param.phi = in_param.paramphi[0];
		
/*		in_param.paramphi[0] = (double)niter;
		in_param.paramphi[2] = fx;
		in_param.paramphi[3] = optite+11;
		
		Rprintf("%.15f %.15f %d_mmgmos\n", in_param.phi, fx, (int)optite+11);
*/		
		
		
		niter++;
		fstart = fx;
	}
}

void calexpression()
{
	int p, i, j, k, q;
	double ai, alphai, c, d_mmgmos, phi;
	double mu_Gauss, var_Gauss, mu_truncGauss, var_truncGauss, kk;
	double secdrv, firstdrv;
	double t1, ym, aym, yym;
	double temp;
	FILE *pf=NULL;
	
	/*pf=fopen("result.txt", "wt");*/
	
	phi = in_param.phi;
	in_param.totalprobe = -1;
	for (p=0; p<in_param.genes; p++)
	{
		c = in_param.parameters[p][in_param.conds+in_param.chips];
		d_mmgmos = in_param.parameters[p][in_param.conds+in_param.chips+1];
				
		getgenedata(p);
		
		if (in_param.num_probe > 1)
		{
			for (i=0; i<in_param.conds; i++)
			{
				alphai = in_param.parameters[p][i];
				ai = in_param.parameters[p][in_param.conds+i];

				t1 = 0.0;
				for (j=0; j<in_param.conds; j++)
					t1 += 2.0*in_param.parameters[p][in_param.conds+j]
					+(1.0+phi)*in_param.parameters[p][j];
				t1 += c;	

				secdrv = in_param.num_probe*((1.0+phi)*(1.0+phi)*trigamma(t1)
					-trigamma(ai+alphai)-phi*phi*trigamma(ai+phi*alphai));
				var_Gauss = -1.0/secdrv;

				if (alphai<1.0e-6)
				{
					yym = 0.0;
					aym = 0.0;
					for (k=0; k<in_param.num_probe; k++)
					{
						ym = 0.0;
						for (q=0; q<in_param.conds; q++)
							ym += in_param.pm[k][q]+in_param.mm[k][q];
						ym += d_mmgmos;
						yym += log(ym);
						aym += log(in_param.pm[k][i])+phi*log(in_param.mm[k][i]);
					}

					firstdrv = in_param.num_probe*((1.0+phi)*digamma(t1)-digamma(ai+alphai)-
						phi*digamma(ai+phi*alphai))-(1.0+phi)*yym+aym;
					mu_Gauss = var_Gauss*firstdrv;
				}
				else
					mu_Gauss = alphai;

				kk=2.0/erfc(-mu_Gauss/sqrt(2.0*var_Gauss));

				mu_truncGauss = kk*(sqrt(var_Gauss)*exp(-mu_Gauss*mu_Gauss/(2.0*var_Gauss))/sqrt(2.0*M_PI)
						+mu_Gauss*erfc(-mu_Gauss/sqrt(2.0*var_Gauss))/2.0);
				var_truncGauss = kk*((var_Gauss+(mu_Gauss-mu_truncGauss)*(mu_Gauss-mu_truncGauss))*erfc(-mu_Gauss/sqrt(2.0*var_Gauss))/2.0
						+sqrt(var_Gauss/(2.0*M_PI))*exp(-mu_Gauss*mu_Gauss/(2.0*var_Gauss))*(mu_Gauss-2.0*mu_truncGauss));

				/* calculate expression value -- the approximate mean of <log(s)> */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+p] = (digamma(mu_truncGauss)+log(d_mmgmos)-digamma(c))/log(2.0)+
										tetragamma(mu_truncGauss)*var_truncGauss/(2.0*log(2.0)*log(2.0));

				/* calculate standard deviation */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+in_param.genes+p] = sqrt(pow(trigamma(mu_truncGauss),2)*var_truncGauss/(log(2.0)*log(2.0)));

				/* calculate percentiles */
				for (j=0; j<in_param.num_prctile; j++)
				{
					temp = mu_Gauss+sqrt(2.0*var_Gauss)*dierfc(1-2.0*(1.0-in_param.prctiles[j])/kk);
					in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+(j+2)*in_param.genes+p] = (digamma(temp)+log(d_mmgmos)-digamma(c))/log(2.0);
				}

			} /* end of for conds */
		}
		else
		{
			for (i=0; i<in_param.conds; i++)
			{
				/* using simple method to calculate the expression level */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+p] = log(0.01 > ((in_param.pm[0][i]-in_param.mm[0][i])/(1.0-phi)) ? 0.01 : ((in_param.pm[0][i]-in_param.mm[0][i])/(1.0-phi)))/log(2.0);

				/* set standard deviation 0 */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+in_param.genes+p] = 0.0;

				/* no percentiles, set all to expression value */
				for (j=0; j<in_param.num_prctile; j++)
					in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+(j+2)*in_param.genes+p] = in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+p];
			}
		}
		if  ((int)p%500 == 0)
			Rprintf(".");
		
	} /* end of for genes */
	Rprintf("\n");
}

void workout_mgmos()
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
	int p, i, j;
	double aa, alphaii, c, d_mgmos;
	double mu_Gauss, var_Gauss, mu_truncGauss, var_truncGauss, kk;
	double secdrv, firstdrv, temp;
	FILE *df=NULL;

	in_param.totalprobe = -1;
	
	if (in_param.saveflag == TRUE)
	{
		df = fopen("par_mgmos.txt", "a");
		
		if (!df)
			Rprintf("Cannot open file for saving parameters\n");
	}
	
	for (p=0; p<in_param.genes; p++)
	{
		/* do optimisation for each gene */
		in_param.cur_gene = p;
			
		getgenedata(p);
	
		for (i=0; i<in_param.chips; i++)
		{	
			in_param.cur_cond = i;
			in_param.flag = 2;
			
			if (in_param.num_probe > 1)
				donlp2();
			
/*			Rprintf("%d_mmgmos %d_mmgmos  %d_mmgmos %d_mmgmos %d_mmgmos\n", niter, p+1, icf, icgf, (int)optite+11);*/
			
			if (in_param.num_probe > 1)
 			{
				aa = in_param.par_mgmos[0];
				alphaii = in_param.par_mgmos[1];
				c = in_param.par_mgmos[2];
				d_mgmos = in_param.par_mgmos[3];

			        secdrv = in_param.num_probe*(trigamma(2*aa+alphaii+c)-trigamma(aa+alphaii));
			        var_Gauss=1/(-secdrv);

			        if (alphaii < 1.0e-6)
				{
				    double s1=0.0, s2=0.0;
				    for (j=0; j<in_param.num_probe; j++)
				    {
				    	s1 += log(in_param.pm[j][i]);
					s2 += log(in_param.pm[j][i]+in_param.mm[j][i]+d_mgmos);
				    }
			            firstdrv = in_param.num_probe*(digamma(2*aa+alphaii+c)-digamma(aa+alphaii))+s1-s2;
			            mu_Gauss = var_Gauss*firstdrv;
        			}
				else
            				mu_Gauss = alphaii;
        
        			kk=2.0/erfc(-mu_Gauss/sqrt(2.0*var_Gauss));
				
				mu_truncGauss = kk*(sqrt(var_Gauss)*exp(-mu_Gauss*mu_Gauss/(2.0*var_Gauss))/sqrt(2.0*M_PI)
						+mu_Gauss*erfc(-mu_Gauss/sqrt(2.0*var_Gauss))/2.0);
				var_truncGauss = kk*((var_Gauss+(mu_Gauss-mu_truncGauss)*(mu_Gauss-mu_truncGauss))*erfc(-mu_Gauss/sqrt(2.0*var_Gauss))/2.0
						+sqrt(var_Gauss/(2.0*M_PI))*exp(-mu_Gauss*mu_Gauss/(2.0*var_Gauss))*(mu_Gauss-2.0*mu_truncGauss));

				/* calculate expression value -- the approximate mean of <log(s)> */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+p] = (digamma(mu_truncGauss)+log(d_mgmos)-digamma(c))/log(2.0)+
										     tetragamma(mu_truncGauss)*var_truncGauss/(2.0*log(2.0)*log(2.0));
										     
				/* calculate standard deviation */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+in_param.genes+p] = sqrt(pow(trigamma(mu_truncGauss),2)*var_truncGauss/(log(2.0)*log(2.0)));
				
				/* calculate percentiles */
				for (j=0; j<in_param.num_prctile; j++)
				{
					temp = mu_Gauss+sqrt(2.0*var_Gauss)*dierfc(1-2.0*(1.0-in_param.prctiles[j])/kk);
					in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+(j+2)*in_param.genes+p] = (digamma(temp)+log(d_mgmos)-digamma(c))/log(2.0);
				}
				if (in_param.saveflag == TRUE)
					fprintf(df," %f %f %f %f", alphaii, aa, c, d_mgmos);

			}
			else
			{
				/* using simple method to calculate the expression level */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+p] = log(0.01 > (in_param.pm[0][i]-in_param.mm[0][i]) ? 0.01 : (in_param.pm[0][i]-in_param.mm[0][i]))/log(2.0);

				/* set standard deviation 0 */
				in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+in_param.genes+p] = 0.0;

				/* no percentiles, set all to expression value */
				for (j=0; j<in_param.num_prctile; j++)
					in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+(j+2)*in_param.genes+p] = in_param.outp[i*(2+in_param.num_prctile)*in_param.genes+p];
				if (in_param.saveflag == TRUE)
					fprintf(df," %f %f %f %f", 0.0, 0.0, 0.0, 0.0);

			}
 		}
		if ((int)p%500 == 0)
			Rprintf(".");
		if (in_param.saveflag == TRUE)
			fprintf(df,"\n");

	}
	Rprintf("\n");	
	if (in_param.saveflag == TRUE)
		fclose(df);
	
}


SEXP mmgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps)
{
	void R_CheckUserInterrupt(void);

	SEXP dim=NULL;	
	SEXP res=NULL;
	double *p_phis=NULL;
	int i, j;
	const char *geneName=NULL;

/*	Rprintf("initialparams_mmgmos ");
*/	initialparams_mmgmos();
	
	PROTECT(dim = getAttrib(PMmat, R_DimSymbol));
	
	in_param.chips = INTEGER(dim)[1];
	in_param.conds = in_param.chips;
	in_param.probes = INTEGER(dim)[0];
/*	in_param.probes = 160;*/
	
	in_param.genes = INTEGER(ngenes)[0];
/*	in_param.genes = 10;*/
	in_param.num_prctile = INTEGER(nprc)[0];
	p_phis = NUMERIC_POINTER(AS_NUMERIC(phis));
	in_param.phi = p_phis[0];
	in_param.mu = p_phis[1];
	in_param.sigma = p_phis[2];
	
	in_param.data_pm = NUMERIC_POINTER(AS_NUMERIC(PMmat));
	in_param.data_mm = NUMERIC_POINTER(AS_NUMERIC(MMmat));
	in_param.prctiles = NUMERIC_POINTER(AS_NUMERIC(prctiles));
	
	in_param.saveflag = LOGICAL_POINTER(AS_LOGICAL(saveflag))[0];
	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];
	
	/*in_param.eps = 1.0e-6;*/
	
/*	Rprintf("allocatemem_mmgmos ");
*/	allocatemem_mmgmos();
	
	for (i=0; i<in_param.conds; i++)
		in_param.replicates[i] = 1;
	
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
	calparameters();

	PROTECT(res = allocMatrix(REALSXP, in_param.genes*(2+in_param.num_prctile), in_param.conds));
	in_param.outp = NUMERIC_POINTER(AS_NUMERIC(res));
	
	Rprintf("Expression values calculating ");
	calexpression();

/*	freemem_mmgmos();*/
	Rprintf("Done.\n");
	UNPROTECT(2);
	
	return res;
}

SEXP mgmos_c(SEXP PMmat, SEXP MMmat, SEXP ngenes, SEXP probeNames, SEXP phis, SEXP prctiles, SEXP nprc, SEXP saveflag, SEXP eps)
{
	SEXP dim=NULL;	
	SEXP res=NULL;
	double *p_phis=NULL;
	int i, j;
	const char *geneName=NULL;

	initialparams_mmgmos();
	
	PROTECT(dim = getAttrib(PMmat, R_DimSymbol));
	
	in_param.chips = INTEGER(dim)[1];
	in_param.conds = in_param.chips;
	in_param.probes = INTEGER(dim)[0];
	
	in_param.genes = INTEGER(ngenes)[0];
	in_param.num_prctile = INTEGER(nprc)[0];
	p_phis = NUMERIC_POINTER(AS_NUMERIC(phis));
	in_param.phi = p_phis[0];
	in_param.mu = p_phis[1];
	in_param.sigma = p_phis[2];
	
	in_param.data_pm = NUMERIC_POINTER(AS_NUMERIC(PMmat));
	in_param.data_mm = NUMERIC_POINTER(AS_NUMERIC(MMmat));
	in_param.prctiles = NUMERIC_POINTER(AS_NUMERIC(prctiles));
	
	in_param.saveflag = LOGICAL_POINTER(AS_LOGICAL(saveflag))[0];
	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];
	
	/*in_param.eps = 1.0e-6;*/
	
	in_param.probesets = (int*)R_alloc(in_param.genes, sizeof(int));
	for (i=0; i<in_param.genes; i++)
		in_param.probesets[i] = 0;
	
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
	
	PROTECT(res = allocMatrix(REALSXP, in_param.genes*(2+in_param.num_prctile), in_param.conds));
	in_param.outp = NUMERIC_POINTER(AS_NUMERIC(res));
/*	in_param.probes = 160;
	in_param.genes = 10;
*/	
	Rprintf("Model optimising ");
	workout_mgmos();

/*	if (in_param.probesets != NULL) Free(in_param.probesets);
*/	
	Rprintf("Done.\n");
	UNPROTECT(2);
	
	return res;
}

/* **************************************************************************** */
/*                              donlp2-intv size initialization                           */
/* **************************************************************************** */
void user_init_size_mmgmos(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X


    /* problem dimension n = dim(donlp2_x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
       
    if (in_param.flag == 0)
    {
    	n = in_param.conds+in_param.chips+2;
    	nstep = 20;
    }
    else if (in_param.flag == 1)
    {
    	n = 1;
        nstep = 40;
    }
    else
    {
    	n = 4;
	nstep = 20;
    }
    
    nlin   =  0;
    nonlin =  0;
    iterma = 4000;
}

/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_mmgmos(void) {
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
		donlp2_x[i] = 1.0;
		low[i] = ALPHALOW;
		up[i] = big;
	}
	for (i=in_param.conds+1; i<=in_param.conds+in_param.chips; i++)
	{
		donlp2_x[i] = 10.0;
		low[i] = ALOW;
		up[i] = big;
	}
	donlp2_x[in_param.conds+in_param.chips+1] = 1.0;
	low[in_param.conds+in_param.chips+1] = CLOW;
	up[in_param.conds+in_param.chips+1] = big;
	donlp2_x[in_param.conds+in_param.chips+2] = 10.0;
	low[in_param.conds+in_param.chips+2] = DLOW;
	up[in_param.conds+in_param.chips+2] = big;
    }
    else if (in_param.flag == 1)
    {
/*        strcpy(name, "phiopt");
	silent = FALSE;
*/	donlp2_x[1] = in_param.phi;
	low[1] = PHILOW;
	up[1] = PHIUP;
    }
    else
    {
 	donlp2_x[1] = 10.0;
	low[1] = ALOW;
	up[1] = big;
	donlp2_x[2] = 1.0;
	low[2] = ALPHALOW;
	up[2] = big;
	donlp2_x[3] = 1.0;
	low[3] = CLOW;
	up[3] = big;
	donlp2_x[4] = 10.0;
	low[4] = DLOW;
	up[4] = big;
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
void setup_mmgmos(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
   
    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_mmgmos(void) {
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
    else if (in_param.flag == 1)
    	in_param.paramphi[0] = donlp2_x[1];
    else
    {
    	for (i=0; i<n; i++)
		in_param.par_mgmos[i] = donlp2_x[i+1];
    }

    return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_mmgmos(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	double alphaii[MAX_NUM_COND]={0.0}, c, d_mmgmos, t1, t2;
	double ym[MAX_NUM_PROBE]={0.0}, aym[MAX_NUM_PROBE]={0.0}, alphas[MAX_NUM_CHIP]={0.0}, a_mos[MAX_NUM_CHIP]={0.0};
	double s5[MAX_NUM_PROBE]={0.0};
	int i, j, k, p;
	
	*fx = 0.0;
	
	

	if (in_param.flag==0)
	{
		
		c = donlp2_x[in_param.conds+in_param.chips+1];
		d_mmgmos = donlp2_x[in_param.conds+in_param.chips+2];
	
		t1 = 0.0;
		t2 = 0.0;
		for (j=0; j<in_param.conds; j++)
		{
			alphaii[j] = donlp2_x[j+1];
			alphas[j] = donlp2_x[j+1];
			a_mos[j] = donlp2_x[in_param.conds+j+1];
			
			t1 += 2.0*a_mos[j]+(1.0+in_param.phi)*alphas[j];
			t2 += lgammafn(a_mos[j]+alphas[j])+lgammafn(a_mos[j]+in_param.phi*alphas[j]);
		}
		t1 += c;
	
		for (k=0; k<in_param.num_probe; k++)
		{
			for (j=0; j<in_param.conds; j++)
			{
				ym[k] += in_param.pm[k][j]+in_param.mm[k][j];
				aym[k] += (a_mos[j]+alphas[j]-1.0)*log(in_param.pm[k][j])+(a_mos[j]+in_param.phi*alphas[j]-1.0)*log(in_param.mm[k][j]);
			}
			ym[k] += d_mmgmos;
		
			*fx += c*log(d_mmgmos) + lgammafn(t1) - lgammafn(c) - t2 - t1*log(ym[k]) + aym[k];
		}
	
		*fx = -*fx;
	}
	else if (in_param.flag == 1)
	{
		in_param.totalprobe = -1;
		for (p=0; p<in_param.genes; p++)
		{
			getgenedata(p);
			if (in_param.num_probe > 1)
			{
				c = in_param.parameters[p][in_param.conds+in_param.chips];
				d_mmgmos = in_param.parameters[p][in_param.conds+in_param.chips+1];
		
				t1 = 0.0;
				t2 = 0.0;
				for (j=0; j<in_param.conds; j++)
				{
					alphaii[j] = in_param.parameters[p][j];
					alphas[j] = in_param.parameters[p][j];
		
					a_mos[j] = in_param.parameters[p][in_param.conds+j];
					t1 += 2.0*a_mos[j]+(1.0+donlp2_x[1])*alphas[j];
					t2 += lgammafn(a_mos[j]+alphas[j])+lgammafn(a_mos[j]+donlp2_x[1]*alphas[j]);
				}
				t1 += c;
		
				for (k=0; k<in_param.num_probe; k++)
				{
					ym[k] = 0.0;
					aym[k] = 0.0;
				}
				for (k=0; k<in_param.num_probe; k++)
				{
					for (i=0; i<in_param.chips; i++)
					{
						ym[k] += in_param.pm[k][i]+in_param.mm[k][i];
						aym[k] += (a_mos[i]+alphas[i]-1.0)*log(in_param.pm[k][i])+(a_mos[i]+donlp2_x[1]*alphas[i]-1.0)*log(in_param.mm[k][i]);
					}
					ym[k] += d_mmgmos;
			
					*fx += c*log(d_mmgmos) + lgammafn(t1) - lgammafn(c) - t2 - t1*log(ym[k]) + aym[k];
				}
			}
		}
	
		*fx += -pow(log(donlp2_x[1])-in_param.mu,2.0)/(2.0*pow(in_param.sigma,2.0))-log(donlp2_x[1]*in_param.sigma*sqrt(2.0*M_PI));
	
		*fx = -*fx;
	}
	else
	{
		double aa, alphaii, c, d_mgmos, s1=0.0, s2=0.0, s3=0.0;
		
		aa = donlp2_x[1];
		alphaii = donlp2_x[2];
		c = donlp2_x[3];
		d_mgmos = donlp2_x[4];
		for (i=0; i<in_param.num_probe; i++)
		{
			s1 += log(in_param.pm[i][in_param.cur_cond]);
			s2 += log(in_param.mm[i][in_param.cur_cond]);
			s3 += log(in_param.pm[i][in_param.cur_cond]+in_param.mm[i][in_param.cur_cond]+d_mgmos);
		}
		
		*fx = in_param.num_probe*(c*log(d_mgmos)+lgammafn(2*aa+alphaii+c)-lgammafn(aa+alphaii)-lgammafn(aa)-lgammafn(c))
		      + (aa+alphaii-1)*s1+(aa-1)*s2-(2*aa+alphaii+c)*s3;
		*fx = -*fx;
		
/*		printf("%.10f %.10f %.10f %.10f %d %d\n",*fx, s1, s2, s3, in_param.num_probe,in_param.cur_cond);*/
	}
	
    return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_mmgmos(DDOUBLE donlp2_x[],DDOUBLE gradf[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	double alphaii[MAX_NUM_COND]={0.0}, s1[MAX_NUM_COND]={0.0}, s2[MAX_NUM_PROBE][MAX_NUM_COND]={0.0}, c, d_mmgmos, t1, s3, s4;
	double ym[MAX_NUM_PROBE]={0.0}, aym[MAX_NUM_PROBE]={0.0}, alphas[MAX_NUM_CHIP]={0.0}, a_mos[MAX_NUM_CHIP]={0.0};
	double s5[MAX_NUM_PROBE]={0.0};
	int i, j, k, p;
    
	if (in_param.flag==0)
	{
		for (j=0; j<in_param.conds+in_param.chips+2; j++)
			gradf[j+1] = 0.0;
		
		c = donlp2_x[in_param.conds+in_param.chips+1];
		d_mmgmos = donlp2_x[in_param.conds+in_param.chips+2];
	
		t1 = 0.0;
		for (j=0; j<in_param.conds; j++)
		{
			alphaii[j] = donlp2_x[j+1];
			s1[j] = 0.0;
			alphas[j] = donlp2_x[j+1];
			a_mos[j] = donlp2_x[in_param.conds+j+1];
			
			s1[j] = digamma(a_mos[j]+alphas[j])+in_param.phi*digamma(a_mos[j]+in_param.phi*alphas[j]);
		
			t1 += 2.0*a_mos[j]+(1.0+in_param.phi)*alphas[j];
		}
		t1 += c;
		
		for (k=0; k<in_param.num_probe; k++)
		{
			for (j=0; j<in_param.conds; j++)
			{
				s2[k][j] = log(in_param.pm[k][j])+in_param.phi*log(in_param.mm[k][j]);
			
				ym[k] += in_param.pm[k][j]+in_param.mm[k][j];
				aym[k] += (a_mos[j]+alphas[j]-1.0)*log(in_param.pm[k][j])+(a_mos[j]+in_param.phi*alphas[j]-1.0)*log(in_param.mm[k][j]);
			}
			ym[k] += d_mmgmos;
		
			for (j=0; j<in_param.conds; j++)
			{
				gradf[j+1] += in_param.replicates[j]*(1.0+in_param.phi)*digamma(t1)-s1[j]-in_param.replicates[j]*(1.0+in_param.phi)*log(ym[k])+s2[k][j];
				gradf[in_param.conds+j+1] += 2.0*digamma(t1)-digamma(a_mos[j]+alphas[j])-digamma(a_mos[j]+in_param.phi*alphas[j])-2.0*log(ym[k])+log(in_param.pm[k][j])+log(in_param.mm[k][j]);
			}
			gradf[in_param.conds+in_param.chips+1] += log(d_mmgmos)+digamma(t1)-digamma(c)-log(ym[k]);
			gradf[in_param.conds+in_param.chips+2] += c/d_mmgmos-t1/ym[k];
		}
	
		for (j=0; j<in_param.conds+in_param.chips+2; j++)
			gradf[j+1] = -gradf[j+1];
	}
	else if (in_param.flag == 1)
	{
		gradf[1] = 0.0;
		
		in_param.totalprobe = -1;
		for (p=0; p<in_param.genes; p++)
		{
			getgenedata(p);
			if (in_param.num_probe > 1)
			{
				c = in_param.parameters[p][in_param.conds+in_param.chips];
				d_mmgmos = in_param.parameters[p][in_param.conds+in_param.chips+1];
		
				t1 = 0.0;
				s4 = 0.0;
				s3 = 0.0;
				for (j=0; j<in_param.conds; j++)
				{
					alphaii[j] = in_param.parameters[p][j];
					alphas[j] = in_param.parameters[p][j];
		
					a_mos[j] = in_param.parameters[p][in_param.conds+j];
					t1 += 2.0*a_mos[j]+(1.0+donlp2_x[1])*alphas[j];
					
					s4 += alphas[j]*digamma(a_mos[j]+donlp2_x[1]*alphas[j]);
					s3 += alphas[j];
				}
				t1 += c;
/*				s3 = vec_sum(alphas,in_param.chips);*/
			
				for (k=0; k<in_param.num_probe; k++)
				{
					ym[k] = 0;
					s5[k] = 0;
				}
				for (k=0; k<in_param.num_probe; k++)
				{
					for (i=0; i<in_param.chips; i++)
					{
						ym[k] += in_param.pm[k][i]+in_param.mm[k][i];
						s5[k] += alphas[i]*log(in_param.mm[k][i]); 
					}
					ym[k] += d_mmgmos;
		
					gradf[1] += s3*digamma(t1)-s4-s3*log(ym[k])+s5[k];
				}
			}
		}
	
		gradf[1] += -(log(donlp2_x[1])-in_param.mu)/(pow(in_param.sigma,2.0)*donlp2_x[1])-1.0/donlp2_x[1];
	
		gradf[1] = -gradf[1];
	}
	else
	{
		double aa, alphaii, c, d_mgmos, s1=0.0, s2=0.0, s3=0.0, s4=0.0;
		
		aa = donlp2_x[1];
		alphaii = donlp2_x[2];
		c = donlp2_x[3];
		d_mgmos = donlp2_x[4];
		for (i=0; i<in_param.num_probe; i++)
		{
			s1 += log(in_param.pm[i][in_param.cur_cond]);
			s2 += log(in_param.mm[i][in_param.cur_cond]);
			s3 += log(in_param.pm[i][in_param.cur_cond]+in_param.mm[i][in_param.cur_cond]+d_mgmos);
			s4 += 1/(in_param.pm[i][in_param.cur_cond]+in_param.mm[i][in_param.cur_cond]+d_mgmos);
		}
		gradf[1] = in_param.num_probe*(2.0*digamma(2.0*aa+alphaii+c)-digamma(aa+alphaii)-digamma(aa))+s1+s2-2.0*s3; 
		gradf[2] = in_param.num_probe*(digamma(2.0*aa+alphaii+c)-digamma(aa+alphaii))+s1-s3;
		gradf[3] = in_param.num_probe*(log(d_mgmos)+digamma(2.0*aa+alphaii+c)-digamma(c))-s3;
		gradf[4] = in_param.num_probe*c/d_mgmos-(2.0*aa+alphaii+c)*s4;
		
		for (i=1; i<5; i++)
			gradf[i] = -gradf[i];
	}
    return;
}

/* **************************************************************************** */
/*  no nonlinear constraints */
/* **************************************************************************** */
void econ_mmgmos(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[], 
              LLOGICAL err[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_mmgmos(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_mmgmos(IINTEGER mode) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}
