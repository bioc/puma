#include <math.h>
#include "o8para.h"
#include "ipplr_c.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

static ipplrparam in_param;
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

void econ_ipplr(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_ipplr(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_ipplr(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_ipplr(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_ipplr(IINTEGER mode);
void freemem_ipplr();
void initialparams_ipplr();
void setup_ipplr();
void solchk_ipplr();
void user_init_ipplr(void);
void user_init_size_ipplr(void);
void allocatemem_ipplr();

void initialparams_ipplr()
{
	in_param.replicates = NULL;
	in_param.data_m = NULL;
	in_param.data_var = NULL;
	in_param.expr = NULL;
	in_param.var = NULL;
	in_param.outp = NULL;
	in_param.mu1 = NULL;
	in_param.mu1sq = NULL;

	econ = econ_ipplr;
	econgrad = econgrad_ipplr;
	ef = ef_ipplr;
	egradf = egradf_ipplr;
	eval_extern = eval_extern_ipplr;
	freemem = freemem_ipplr;
	initialparams = initialparams_ipplr;
	setup = setup_ipplr;
	solchk = solchk_ipplr;
	user_init = user_init_ipplr;
	user_init_size = user_init_size_ipplr;
	allocatemem = allocatemem_ipplr;
}


void allocatemem_ipplr()
{
	in_param.expr = (double*)R_alloc(in_param.chips, sizeof(double));
	in_param.var = (double*)R_alloc(in_param.chips, sizeof(double));
	in_param.mu1 = (double*)R_alloc(in_param.conds, sizeof(double));
	in_param.mu1sq = (double*)R_alloc(in_param.conds, sizeof(double));
}

void freemem_ipplr()
{
	if (in_param.expr != NULL) Free(in_param.expr);
	if (in_param.var != NULL) Free(in_param.var);
	if (in_param.mu1 != NULL) Free(in_param.mu1);
	if (in_param.mu1sq != NULL) Free(in_param.mu1sq);

}

void findeforc_ipplr(double *x, int c, double *e_c, int *num_c)
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

double mean_ipplr(double *x, int len)
{
	int i;
	double s=0.0;

	for (i=0; i<len; i++)
		s += x[i];
	return s/len;
}

double var_ipplr(double *x, int len)
{
	int i;
	double m,s=0.0;

	m = mean_ipplr(x,len);
	for (i=0; i<len; i++)
		s += (x[i]-m)*(x[i]-m);

	return s/(len-1);
}

double fmaxn_ipplr(double *x, int n)
{
	int i;
	double m;

	m = x[0];
	for (i=1; i<n; i++)
		if (x[i]>m)
			m = x[i];

	return m;
}


void workout0()
{
       #define  X extern
       #include "o8comm.h"
       #undef   X
       #include "o8cons.h"

        int i, j,k,  num_c;
        double exp_c[MAX_NUM_REPLICATE], var_c[MAX_NUM_REPLICATE];
        double mu_temp[MAX_NUM_COND],v[MAX_NUM_COND],lam_temp[MAX_NUM_COND],lamda10[MAX_NUM_COND];		
        double mu1old[MAX_NUM_COND],diff_mu1[MAX_NUM_COND],lamda1[MAX_NUM_COND],lamda1old[MAX_NUM_COND],diff_lamda[MAX_NUM_COND];
	    double x_exp[MAX_NUM_CHIP],x_var[MAX_NUM_CHIP],x_sq[MAX_NUM_CHIP];
        double temp,tempx,alpha1,beta1,fopt,foptold,Fx1,Fx2;


           

               /*FILE *temp_f=NULL;
	        temp_f=fopen("temp_folder.txt","wt");*/

	for (i=0; i<in_param.genes; i++)
	 {
		for (j=0; j<in_param.chips; j++)
		{
		        in_param.expr[j] = in_param.data_m[j*in_param.genes+i];
			in_param.var[j] = in_param.data_var[j*in_param.genes+i];

		}
		

		for (j=1; j<=in_param.conds; j++)
		{
			findeforc_ipplr(in_param.expr, j, exp_c, &num_c);
			findeforc_ipplr(in_param.var, j, var_c, &num_c);
			
      
			mu_temp[j-1] = mean_ipplr(exp_c, num_c);
			if (num_c == 1)
			    temp = 0.1;
			else
		            temp = var_ipplr(exp_c, num_c);
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
			lam_temp[j-1] = 1/mean_ipplr(v, num_c);

		}
		

		in_param.mu0 = mean_ipplr(mu_temp, in_param.conds);
		temp = var_ipplr(mu_temp, in_param.conds);
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
			temp = mean_ipplr(lam_temp, in_param.conds);
			tempx = var_ipplr(lam_temp, in_param.conds);
			if (tempx == 0.0)
			    tempx = 0.1;
			in_param.alpha0 = temp*temp/tempx;
                        in_param.beta0 = temp/tempx;   

               if(in_param.beta0<0.001||in_param.beta0>100){
                           in_param.alpha0 = 10.0;
                           in_param.beta0 = 1.0;
                       }                         
                }
            //  fprintf(temp_f,"%d\t%f\t%f\t%f\t%f\n",i+1,in_param.mu0,in_param.eta0,in_param.alpha0,in_param.beta0);

                int m=0,n=0;
                foptold = POSI_INF;
		fopt = 1.0e18;
   
		while (foptold-fopt > in_param.eps && m<in_param.max_num)
		{
                      m++;
                      n=0;

		      for(j=0;j<in_param.conds;j++)
                         {
			 in_param.mu1[j]=mu_temp[j];
			 lamda1[j]=lamda10[j];
			 mu1old[j]=POSI_INF;
			 lamda1old[j]=POSI_INF;
			 diff_mu1[j]=fabs(in_param.mu1[j]-mu1old[j]);
			 diff_lamda[j]=fabs(lamda1[j]-lamda1old[j]);
                         
			 }

                      alpha1 = in_param.alpha0+in_param.chips/2.0;
                      beta1 = in_param.beta0+0.15;
                      in_param.lamda_m = alpha1/beta1;
		

		   while(fmaxn_ipplr(diff_mu1,in_param.conds)>in_param.eps||fmaxn_ipplr(diff_lamda,in_param.conds)>in_param.eps&&n<2000)
		       {
                         n++;
			      
                         for(j=0;j<in_param.conds;j++)
                            in_param.mu1sq[j]=in_param.mu1[j]*in_param.mu1[j]+1.0/lamda1[j];
                             
                         for(j=0;j<in_param.chips;j++)
                         {                       
                            temp=1+in_param.var[j]*in_param.lamda_m;
                            x_exp[j]=(in_param.expr[j]+in_param.var[j]*in_param.mu1[in_param.replicates[j]-1]*in_param.lamda_m)/temp;
			    x_var[j]=in_param.var[j]/temp;
                            x_sq[j]=x_exp[j]*x_exp[j]+x_var[j];                       
                         }


                            temp=0.0;                           
                         for(j=0;j<in_param.chips;j++)
			    temp+=(x_sq[j]+in_param.mu1sq[in_param.replicates[j]-1]-2.0*x_exp[j]*in_param.mu1[in_param.replicates[j]-1])/2.0;
                            
                            alpha1=in_param.alpha0+in_param.chips/2.0;
                            beta1=in_param.beta0+temp;
                            in_param.lamda_m=alpha1/beta1;

                         for(j=0;j<in_param.conds;j++)
                         {
                            lamda1old[j]=lamda1[j];
                            mu1old[j]=in_param.mu1[j];
                                       
                            findeforc_ipplr(x_exp,j+1,exp_c,&num_c);

                            lamda1[j]=in_param.eta0+num_c*in_param.lamda_m;
                                       
                            in_param.mu1[j]=in_param.eta0*in_param.mu0/lamda1[j];
                            for(k=0;k<num_c;k++)
                                in_param.mu1[j]+=in_param.lamda_m*exp_c[k]/lamda1[j];
                                        
                            diff_mu1[j]=fabs(in_param.mu1[j]-mu1old[j]);
                            diff_lamda[j]=fabs(lamda1[j]-lamda1old[j]);
  
                          }  
                                                           

                    }     
                           in_param.loglamda = digamma(alpha1)-log(beta1);

                           for(j=0;j<in_param.conds;j++)
                               in_param.mu1sq[j]=in_param.mu1[j]*in_param.mu1[j]+1.0/lamda1[j];
                              
                        

                  Fx1=0.0; 
                  Fx2=0.0;  

                  for(j=0;j<in_param.chips;j++)
                  {
                        Fx1-=((in_param.expr[j]*in_param.expr[j]-2.0*in_param.expr[j]*x_exp[j]+x_sq[j])/in_param.var[j]+log(2.0*M_PI)+log(in_param.var[j]))/2.0;
                        Fx1-=(in_param.lamda_m*(x_sq[j]-2.0*x_exp[j]*in_param.mu1[in_param.replicates[j]-1]+in_param.mu1sq[in_param.replicates[j]-1])-in_param.loglamda+log(2.0*M_PI))/2.0;

                  }

                 for(j=0;j<in_param.chips;j++)
                        Fx2-=(1+log(2.0*M_PI)+log(x_var[j]))/2.0;

                 for(j=0;j<in_param.conds;j++)               
                        Fx2-=(1+log(2.0*M_PI)-log(lamda1[j]))/2.0;

                 Fx2-=alpha1*log(beta1)-lgammafn(alpha1)+(alpha1-1)*in_param.loglamda-beta1*in_param.lamda_m;


                 foptold = fopt;

	         donlp2();

                 fopt = fx-Fx1+Fx2;                
              
         }       
               // Rprintf("Gene->%d: %d\t%d\t%d\t%f\n",i+1,(int)optite+11,m,n,fopt);

           for(j=0;j<in_param.conds;j++)
	       { 
	          in_param.outp[j*in_param.genes+i] = in_param.mu1[j];
	          in_param.outp[(in_param.conds+j)*in_param.genes+i] = sqrt(1.0/lamda1[j]);
	       }

//	   if ((int)i%200 == 0)
//		   Rprintf(".");

    }
       // fclose(temp_f);
 }

/***********************************************************************************************************************/
/*----------------------------------------interfac    ----------------------------------------------  */
/***********************************************************************************************************************/

SEXP hcomb_c(SEXP Mmat, SEXP Stdmat, SEXP rep,  SEXP conds,  SEXP max_num,  SEXP eps )
{
	SEXP dim=NULL;
	SEXP res=NULL;

	initialparams_ipplr();

	PROTECT(dim = getAttrib(Mmat, R_DimSymbol));
	in_param.genes = INTEGER(dim)[0];
	in_param.chips = INTEGER(dim)[1];

	in_param.conds = INTEGER_POINTER(AS_INTEGER(conds))[0];

	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];

	in_param.data_m = NUMERIC_POINTER(AS_NUMERIC(Mmat));
	in_param.data_var = NUMERIC_POINTER(AS_NUMERIC(Stdmat));

	in_param.replicates = INTEGER_POINTER(AS_INTEGER(rep));
             in_param.max_num = INTEGER_POINTER(AS_INTEGER(max_num))[0];

	allocatemem_ipplr();
	PROTECT(res = allocMatrix(REALSXP, in_param.genes, 2*in_param.conds));
	in_param.outp = NUMERIC_POINTER(AS_NUMERIC(res));

	//Rprintf("Combined signals are computing\n ");
    
		workout0();
           
	/*freemem_ipplr();*/
	/*Rprintf("Done.\n");*/
	UNPROTECT(2);

	return res;
}

/* **************************************************************************** */
/*                              donlp2-intv size initialization                           */
/* **************************************************************************** */
void user_init_size_ipplr(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X

   /* problem dimension n = dim(donlp2_x),
          nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */

     n=4;
     nstep=20;

     nlin   =  0;
     nonlin =  0;
     iterma = 4000;
}


/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_ipplr(void)
{
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X

    static IINTEGER i;
    silent=TRUE;

    big=1.0e20;


    donlp2_x[1] = in_param.mu0;
    donlp2_x[2] = in_param.eta0;
    donlp2_x[3] = in_param.alpha0;
    donlp2_x[4] = in_param.beta0;


    for(i=1;i<=4;i++)
     {
       low[i]=LOWBOUND;
       up[i]=big;
     }
    low[1]=-big;

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
void setup_ipplr(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_ipplr(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"

    in_param.mu0 = donlp2_x[1];
    in_param.eta0 = donlp2_x[2];
    in_param.alpha0 = donlp2_x[3];
    in_param.beta0 = donlp2_x[4];

    return;

 }

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_ipplr(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    int i;
    double alpha_x,beta_x,mu_x,eta_x;


    *fx=0.0;

    mu_x = donlp2_x[1];
    eta_x = donlp2_x[2];
    alpha_x = donlp2_x[3];
    beta_x = donlp2_x[4];


    *fx=alpha_x*log(beta_x)+(alpha_x-1)*in_param.loglamda-beta_x*in_param.lamda_m-lgammafn(alpha_x)+in_param.conds*(log(eta_x)-log(2.0*M_PI))/2.0;

    for(i=0;i<in_param.conds;i++)     
         *fx-=eta_x*(in_param.mu1sq[i]-2.0*in_param.mu1[i]*mu_x+mu_x*mu_x)/2.0;
  
     *fx=-*fx;

    return;

}
/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_ipplr(DDOUBLE donlp2_x[],DDOUBLE gradf[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    int i;
    double alpha_x,beta_x,mu_x,eta_x;


    mu_x = donlp2_x[1];
    eta_x = donlp2_x[2];
    alpha_x = donlp2_x[3];
    beta_x = donlp2_x[4];

    gradf[1] = 0.0;
    gradf[2] = -in_param.conds/(2.0*eta_x);             
	
    for(i=0;i<in_param.conds;i++)
       {
	 gradf[1]+=(mu_x-in_param.mu1[i])*eta_x;
         gradf[2]+=(in_param.mu1sq[i]-2.0*in_param.mu1[i]*mu_x+mu_x*mu_x)/2.0;
       }

    gradf[3] = digamma(alpha_x)-log(beta_x)-in_param.loglamda;
    gradf[4] = in_param.lamda_m-alpha_x/beta_x;

    return;

}



/* **************************************************************************** */
/*  no nonlinear constraints */
/* **************************************************************************** */
void econ_ipplr(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[],
              LLOGICAL err[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_ipplr(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_ipplr(IINTEGER mode) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}
