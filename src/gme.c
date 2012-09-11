/* **************************************************************************** */
/*       gme C implementation                                          */
/* **************************************************************************** */
#include <math.h>
#include "o8para.h"
#include "gmehead.h"
#include "global_constants.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Memory.h>

static gme_expparam in_param;
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

void econ_gme(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_gme(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_gme(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_gme(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_gme(IINTEGER mode);
void freemem_gme();
void initialparams_gme();
void setup_gme();
void solchk_gme();
void user_init_gme(void);
void user_init_size_gme(void);
void allocatemem_gme();

int pmst=0;
int gtst=0;
double **PM;
double **GTM;
int  numofalpha=0;
int loca=0;


void gme_getgenedata(int g)  
{
	int i,j,k;
       
           
	in_param.cur_num_gene=in_param.num_gene[g];
	in_param.cur_num_alpha=in_param.num_alpha[g];
          in_param.cur_num_probe=in_param.num_probe[g];
	in_param.conds=in_param.cur_num_alpha*in_param.chips;
	
         
	for( i=pmst; i<pmst+in_param.cur_num_probe; i++)  //pm of gene
	{
		
		for(j=0; j<in_param.chips+1; j++)
		{
			in_param.pm[i-pmst][j]=PM[i][j];
                             Rprintf("%f",in_param.pm[i-pmst][j]);
		}
                   
                   Rprintf("\n");
           
	}
          
          pmst=pmst+in_param.cur_num_probe;
          
           
 
	for(k=gtst; k<gtst+in_param.cur_num_gene;k++) //corresponding of gene
	{	
		
		for(j=0;j<in_param.gtdim;j++)
		{
			in_param.gt[k-gtst][j]=GTM[k][j];
                             
                              
		}
             
	}
       
	gtst=gtst+in_param.cur_num_gene;
                 
	for(i=0; i<in_param.cur_num_probe; i++)  
	{
		for(j=0; j<in_param.cur_num_gene; j++)
		{     
                         
	             if(in_param.gt[j][1]==in_param.pm[i][0])
                          { 
			   in_param.MB[i][ (int)in_param.gt[j][0]-1]=1.0;
                               for(k=1;k<in_param.chips; k++)
               		     {
			            in_param.MB[i][(int)in_param.gt[j][0]-1+in_param.cur_num_alpha*k]=1.0;
                               
			      }
	               }
                  
                    }
                 
      
	}
   


}


void initialparams_gme()
{
  in_param.data_pm=NULL;
  in_param.data_gt=NULL;
  in_param.data_GT=NULL;
  in_param.data_GN=NULL;
  in_param.data_PN=NULL;
  in_param.data_AN=NULL;
  in_param.parameters=NULL;
 
       econ = econ_gme;
	econgrad = econgrad_gme;
	ef = ef_gme;
	egradf = egradf_gme;
	eval_extern = eval_extern_gme;
	freemem = freemem_gme;
	initialparams = initialparams_gme;
	setup = setup_gme;
	solchk = solchk_gme;
	user_init = user_init_gme;
	user_init_size = user_init_size_gme;
	allocatemem = allocatemem_gme;
  
}



void allocatemem_gme()
{ 
   int j;
  
  
    in_param.parameters=(double**)R_alloc(in_param.numofgenes,sizeof(double*));
    
  
    for(j=0; j<in_param.numofgenes;j++)
    {
     
             in_param.parameters[j]=(double *)R_alloc(in_param.num_alpha[j]*in_param.chips+2,sizeof(double));
    }
   
}


void freemem_gme()
{
    int j;
    for(j=0;j<in_param.numofgenes;j++)
    {
         if(in_param.parameters[j]!=NULL) Free(in_param.parameters[j]);
    }
    if(in_param.parameters!=NULL) Free(in_param.parameters);
}

double  ** allocate_matrix(int m,int n) 
{
	int i,j;
	double **c=NULL;
	c=(double **)Calloc(m,double);
       for(i=0;i<m;i++)
              c[i]=(double *)Calloc(n,double);
	for(i=0;i<m;i++)
	{	
		for(j=0;j<n;j++)
		{ c[i][j]=0;}
	}
	return c;

}


void gme_calparameters()
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
	int niter = 1, nx;
	double fstart;
	int p, i, j;
	int finishflag = 1;
     
  	fstart = HUGE_VAL;
	
	nx = in_param.conds+2;
    
	
	
	while(1)
	{
		void R_CheckUserInterrupt(void);
		
		
             
		for (p=0;p<in_param.numofgenes; p++)
		{
			void R_CheckUserInterrupt(void);
			
			/* do optimisation for each gene*/
			in_param.cur_gene = p;
			
			gme_getgenedata(p);
                    	donlp2();
                             
			if ((int)p%1000 == 0)
			   Rprintf(".");
                              
                                
                              
  		}
		/* check the termination criteria*/
                    
		  if (in_param.saveflag ==1)
			{
				FILE *df = fopen("par_gmoexon.txt", "wt");
		
				if (!df)
				{ 
					Rprintf("Cannot open file for saving parameters\n");
					break;
				}
		
				for (i=0; i<in_param.numofgenes; i++)
				{
					for (j=0; j<in_param.num_alpha[i]*in_param.chips+2; j++)
						fprintf(df," %f", in_param.parameters[i][j]);
                                                       
					
				}
				fclose(df);
					
		      }
                     
			Rprintf("\n");
                       
			break;
	
	 }
}





/* **************************************************************************** */
/*                              donlp2-intv size initialization                 */
/* **************************************************************************** */
void user_init_size_gme(void){
   
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X

    
    /* problem dimension n = dim(donlp2_x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
       
    
  
    n = in_param.conds+2;
    nstep = 20;
   
    
    
    nlin   =  0;
    nonlin =  0;
    iterma = 4000;//
}

/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_gme(void) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X
    
    static IINTEGER i;

    silent = TRUE;

    
    big = 1.e20;
 	
	/* initialise the parameters */
	for (i=1; i<=in_param.conds; i++)
	{
		donlp2_x[i] = 2.0;
		low[i] = ALOW;
		up[i] = big;
	}
	
	donlp2_x[in_param.conds+1] = 6.0;
	low[in_param.conds+1] = CLOW;
	up[in_param.conds+1] = big;
	donlp2_x[in_param.conds+2] = 10.0;
	low[in_param.conds+2] = DLOW;
	up[in_param.conds+2] = big;
   

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
void setup_gme(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
   
    return;
}
/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_gme(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
    
   int i;
   for (i=0; i<n; i++)
    {    in_param.parameters[in_param.cur_gene][i] = donlp2_x[i+1];}
       
  
    return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_gme(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

	double  c,d,q=0.0,w=0.0, t1=0.0, t2=0.0,t3=0.0;
	int i, j, k,k1;
	
	*fx = 0.0;
	for (j=0; j<in_param.conds; j++)
	{	
	   in_param.alphaii[j] = donlp2_x[j+1];
                                         
	 }
     
         c = donlp2_x[in_param.conds+1];
         d = donlp2_x[in_param.conds+2];
         
	  
	for(k1=0;k1<in_param.cur_num_probe; k1++)
	{        
               q=0.0;
	     w=0.0;
	     t3=0.0;	
               t2=0.0; 
        
           for(i=0; i<in_param.conds; i++)
	       in_param.alpha[i]=in_param.alphaii[i]*in_param.MB[k1][i];
	                        	
           for(i=0; i<in_param.conds; i++)
	       { 
                 
                  q+=in_param.alpha[i];
                 }
                   
	  q+=c;
		
	  for(j=0; j<in_param.chips; j++)
	      {
                 
                  t1=0.0;
                  t3+=log(in_param.pm[k1][j+1]);
	           
                  for(k=j*in_param.cur_num_alpha; k<(j+1)*in_param.cur_num_alpha; k++)
	          {
                        t1+=in_param.alpha[k];
                    }
			
                  t2+=(t1)*log(in_param.pm[k1][j+1])-lgammafn(t1);
	        w+=in_param.pm[k1][j+1];
	      }
                
              w+=d;
              
		
	  *fx += c*log(d) + lgammafn(q) + t2 - lgammafn(c)-q*log(w)-t3;
		
		
      }
       
         
                  
	*fx=-*fx;
   

}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_gme(DDOUBLE donlp2_x[],DDOUBLE gradf[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
       
    int i, j, k,k2;
    double c, d,q=0.0,w=0.0;
    double t4=0.0,t5=0.0;

   
	
     for (j=0; j<in_param.conds+2; j++)
	    gradf[j+1] = 0.0;

     for (j=0; j<in_param.conds; j++)
	 {	
	       in_param.alphaii2[j] = donlp2_x[j+1];
	 }
	
           c = donlp2_x[in_param.conds+1];
	 d = donlp2_x[in_param.conds+2];

      for(k2=0; k2<in_param.cur_num_probe; k2++)
	   {
		   
	           q=0.0;
		 w=0.0;
		
	         for(i=0; i<in_param.conds; i++)
			   in_param.alpha2[i]=in_param.alphaii2[i]*in_param.MB[k2][i];


	         for(i=0; i<in_param.conds; i++)
		     {  
              
                           q+=in_param.alpha2[i];
                         }
		   
	          q+=c;
		    
                    for(j=0; j<in_param.chips; j++)
	             
                       {
                          w+=in_param.pm[k2][j+1];
                        }
		  
                    w+=d;
		

		 for(j=0; j<in_param.chips; j++)
		    {
              
		       t4=0.0;
	                 t5=0.0;
				
                        for(k=j*in_param.cur_num_alpha; k<(j+1)*in_param.cur_num_alpha; k++)
			 
                               t4+=in_param.alpha2[k];
		           t5=log(in_param.pm[k2][1+j])-log(w)-digamma(t4)+digamma(q);
		           
                        for(k=j*in_param.cur_num_alpha; k<(j+1)*in_param.cur_num_alpha; k++)
		                   
                                in_param.Dif[k2][k]=t5;
		     }
            
		   gradf[in_param.conds+1] += log(d)+digamma(q)-digamma(c)-log(w);
		
                       gradf[in_param.conds+2] += c/d-q/w;
      
	}
            
                
    for(i=0; i<in_param.cur_num_probe; i++)
	{
              for(j=0;j<in_param.conds;j++)
			  
                    in_param.G1[j]=in_param.Dif[i][j]*in_param.MB[i][j];
              
              for(j=0;j<in_param.conds; j++)
                 {
                 	   in_param.G[i][j]=in_param.G1[j];
                 }
           }
             


	for(i=0; i<in_param.conds; i++)
	{
	       for(j=0;j<in_param.cur_num_probe; j++)
		   gradf[i+1]+=in_param.G[j][i];
	}
	
        for (j=0; j<in_param.conds+2; j++)
	     gradf[j+1] = -gradf[j+1];
               
    	
	
     return;
  }

/* **************************************************************************** */
/*  no nonlinear constraints */
/* **************************************************************************** */
void econ_gme(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[], 
              LLOGICAL err[]) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_gme(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_gme(IINTEGER mode) {
    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}

SEXP gme_c(SEXP PMmat, SEXP GTmat, SEXP PNmat, SEXP GNmat, SEXP ANmat, SEXP totalgene, SEXP saveflag, SEXP eps)
{
	void R_CheckUserInterrupt(void);
       

          SEXP dim=NULL;
	SEXP dim11=NULL;
	SEXP res=NULL;
        
	int i, j;
                    
	initialparams_gme();
	
	PROTECT(dim = getAttrib(PMmat, R_DimSymbol));
       
	PROTECT(dim11 =getAttrib(GTmat, R_DimSymbol));
	in_param.chips =INTEGER(dim)[1]-1;
          in_param.probes=INTEGER(dim)[0];
	in_param.genes = INTEGER(dim11)[0];
          in_param.gtdim= INTEGER(dim11)[1];        
	in_param.data_pm = NUMERIC_POINTER(AS_NUMERIC(PMmat));
       
          in_param.data_GT = NUMERIC_POINTER(AS_NUMERIC(GTmat));
          in_param.data_GN = INTEGER_POINTER(AS_INTEGER(GNmat)); 
	       
	in_param.data_PN = INTEGER_POINTER(AS_INTEGER(PNmat));
	in_param.data_AN = INTEGER_POINTER(AS_INTEGER(ANmat));
             
         
	in_param.saveflag = LOGICAL_POINTER(AS_LOGICAL(saveflag))[0];
     
	in_param.eps = NUMERIC_POINTER(AS_NUMERIC(eps))[0];
          in_param.numofgenes= INTEGER_POINTER(AS_INTEGER(totalgene))[0];
          
             
        
          in_param.totalprobe=-1;
          in_param.totalgene=-1;
          numofalpha=0;
     
       for(j=0;j<in_param.numofgenes;j++)
	{
          
            in_param.num_gene[j]=in_param.data_GN[j];
            in_param.num_alpha[j]=in_param.data_AN[j];
            numofalpha=numofalpha+in_param.num_alpha[j]*in_param.chips;
            in_param.num_probe[j]=in_param.data_PN[j];
           
           
         }
       

      // double **PM;
      PM=allocate_matrix(in_param.probes+1,in_param.chips+1);
 
      for( i=0; i<in_param.probes; i++)  //pm of gene
	{
	    in_param.totalprobe++;
	    for(j=0; j<in_param.chips+1; j++)
	     {
		PM[i][j]=in_param.data_pm[j*in_param.probes+in_param.totalprobe];
                    
	     }
    
	}
       // double **GT;
      
       GTM=allocate_matrix(in_param.genes+1,in_param.gtdim);

       for(i=0; i<in_param.genes;i++) //corresponding of gene
	{	
	     in_param.totalgene++;
	     for(j=0;j<in_param.gtdim;j++)
	     {
		GTM[i][j]=in_param.data_GT[j*in_param.genes+in_param.totalgene];
                   
	      }
               
	}
	 

     
       allocatemem_gme();

     
		
       Rprintf("Model optimising now ");
       gme_calparameters();
      
       PROTECT(res = allocMatrix(REALSXP,numofalpha+2*in_param.numofgenes,1));
       in_param.outp = NUMERIC_POINTER(AS_NUMERIC(res));
       	
       loca=0;
       for(i=0;i<in_param.numofgenes;i++)
       {
             for(j=0;j<in_param.num_alpha[i]*in_param.chips+2;j++)
              {  
                                             
		in_param.outp[loca] = in_param.parameters[i][j];
                    loca=loca+1;   
                          
              }

					
        }
         
	
      	
      pmst=0;
      gtst=0;

	//freemem_gme();
      for(i=0;i<in_param.probes+1; i++)
           Free(PM[i]);
      for(i=0;i<in_param.genes+1; i++)
          Free(GTM[i]);
      Rprintf("Done.\n");
      UNPROTECT(3);
	
      return res;
}








