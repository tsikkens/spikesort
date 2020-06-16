
/*

fsmem_mvgm : Free Split and Merge Expectation-Maximization algorithm for Multivariate Gaussian Mixtures

Usage  
-------

[logl , M , S , P]  = fsmem_mvgm(Z , [M0] , [S0] , [P0] , [option]);


Inputs
-------

Z                 Measurements (d x N), where d denotes the measure's dimension and N number of samples
M0                Initial mean vector. M0 can be (d x 1 x K), where K denotes the a priori number of clusters (default [Kini random elements from Z])
S0                Initial covariance matrix. S0 can be (d x d x K)  (default [cov(Z)/40])
P0                Initial mixture probablities (1 x 1 x K) :        (default [1/Kini])
options       
                  Kmin              Minimum number of compounds (default [1])
                  Kini              Initial number of compounds (default [5])
                  Kmax              Maximum number of compounds (default [15])
                  maxite_fsmem      Number of maximum iteration for the main loop of the fsmem (default [100])
                  maxite_fullem     Number of maximum iteration for the full EM inside the main loop (default [100])
                  maxite_partialem  Number of maximum iteration for the partial EM inside the main loop (default [100])
                  epsi_fullem       Tolerance in loglikelihood improuvement of the Full EM (default [1e-6])
                  epsi_partialem    Tolerance in loglikelihood improuvement of the Partial EM (default [1e-6])
                  lambda            Covariance regularization parameter (default [0.01])
                  maxcands_split    Maximum number of split candidate (default [5])
                  splitinit_epsi    Split Initialisation parameter for the mean of splitted cluster (default [1])
                  maxcands_merge    Maximum number of merge candidate (default [5])
                  covtype           Covariance type : 0 = full , 1 = elliptical , 2 = spherical (default [0])
                  fail_exit         Number of tentatives of split/merge operations before exit. If fail_exit = 0, then FSMEM = EM 


Ouputs
-------

logl              Loglikelihood & information matrix (6 x Tfinal) where Tfinal is the total number of fsmem iterations
logl(1,i)         Loglikelihood at iteration i=1,...,Tfinal
logl(2,i)         Loglikelihood of the partial EM at iteration i=1,...,Tfinal
logl(3,i)         Number of iterations of the partial EM
logl(4,i)         Loglikelihood of the full EM at iteration i=1,...,Tfinal
logl(5,i)         Number of iterations of the full EM
logl(6,i)         Flag of the split/merge moves : {-1 , -0.5} = success/failed merge move, {1, 0.5} sucess/failed split move


M                 Estimated mean vector (d x 1 x Kest), where Kest is the number of estimated coupounds
S                 Estimated covariance vector (d x d x Kest)
P                 Estimated initial probabilities (1 x 1 x Kest)


To compile
----------

mex  -DranSHR3 -output fsmem_mvgm.dll fsmem_mvgm.c

or

mex -DranSHR3 -f mexopts_intel10.bat -output fsmem_mvgm.dll fsmem_mvgm.c



Example 1
----------


d                                   = 2;
m                                   = 2;
L                                   = 1;
K                                   = 10000;
nbite                               = 100;

P                                   = cat(3 , [0.4] , [0.6]);
M                                   = cat(3 , [-1 ; -1] , [1 ; 1]);
S                                   = cat(3 , [1 0.3 ; 0.3 0.8] , [0.7 0.6; 0.6 1]);


[Z ,  X]                            = sample_mvgm(K , M , S , P);

[logl , Mest , Sest , Pest]         = fsmem_mvgm(Z);



[x , y]                             = ndellipse(M , S);
[xest , yest]                       = ndellipse(Mest , Sest);

L                                   = likelihood_mvgm(Z , Mest , Sest , Pest);


[val , Xest]                        = max(L);

ind1                                = (Xest == 1);
ind2                                = (Xest == 2);

sum(X~=Xest , 2)/K

figure(1) , plot(Z(1 , ind1) , Z(2 , ind1) , 'k+' , Z(1 , ind2) , Z(2 , ind2) , 'g+' , x , y , 'b' ,xest   , yest  ,'r', 'linewidth' , 2);


[Znew ,  Xnew]                      = sample_mvgm(K , M , S , P);

Lnew                                = likelihood_mvgm(Znew , Mest , Sest , Pest);


[val , Xnewest]                     = max(Lnew);

ind1                                = (Xnewest == 1);
ind2                                = (Xnewest == 2);

sum(Xnew~=Xnewest)/K

figure(2), plot(Znew(1 , ind1) , Znew(2 , ind1) , 'k+' , Znew(1 , ind2) , Znew(2 , ind2) , 'g+' , x , y , 'b' , xest  , yest ,'r', 'linewidth' , 2);



Example 2
----------

clear , clc , close all hidden

options.Kmin                           = 1;
options.Kini                           = 10;
options.Kmax                           = 20;
options.maxite_fullem                  = 500;
options.maxite_partialem               = 500;
options.epsi_fullem                    = 1e-7;
options.epsi_partialem                 = 1e-7;
options.maxcands_split                 = 12;
options.maxcands_merge                 = 12;
options.fail_exit                      = 2;


Z                                      = spiral2d(3000);
tic,[logl , Mest , Sest , Pest]        = fsmem_mvgm(Z , [] , [] , [] , options);,toc
K_est                                  = size(Pest , 3);

[xest , yest]                          = ndellipse(Mest , Sest);
[xaxes , yaxes]                        = ndellipse(Mest , Sest , [] , [] , 3);
xaxes(end , :)                         = [];
yaxes(end , :)                         = [];
h                                      = plot(Z(1 , :) , Z(2 , :) , '+' , xest , yest , 'r' , xaxes , yaxes , 'k', 'linewidth' , 2);
legend([h(1), h(2) , h(2+K_est)] , 'Z' , sprintf('Estimated, K=%d' , K_est) , 'PC');



Author : Sébastien PARIS  © (sebastien.paris@lsis.org)
--------


Ver 1.1 (19/01/09)
------- 


References
----------

[1] Wagenaar, D. A. : FSMEM for MoG, Term project for CS/CNS/EE 156b: Learning systems, 
class by P. Perona and M. Welling, Caltech, 2000

[2] Ueda, N., Nakano, R., Ghahramani, Z., Hinton, G.: SMEM algorithm for mixture
models. Neural Computation 12(9), pp. 2109-2128 (2000)

[3] Zhang B., Zhang C, Yi X. : Competitive EM algorithm for finite mixture models,
Pattern Recognition, Volume 37, Number 1, (January 2004) , pp. 131-144(14)

[4] Blekas K., Lagaris I. E.: Split-Merge Incremental LEarning (SMILE) of Mixture Models. ICANN (2) 2007: pp. 291-300



*/



#include <math.h>
#include <time.h>
#include "mex.h"

#define NUMERICS_FLOAT_MIN 1.0E-37
#define M_PI 3.14159265358979323846

#define tiny 1e-300
#define huge 1e300
#define likely_ini -1e300
#define fudge 1e-6
#define maxite_eigen 100
#define tol_eigen 1e-6


struct opts
{

	int    Kmin;

	int    Kini;

	int    Kmax;

	int    maxite_fsmem;

	int    maxite_fullem;

	int    maxite_partialem;

	double epsi_fullem;

	double epsi_partialem;

	double lambda;

	int    maxcands_split;

	double splitinit_epsi;

	int    maxcands_merge;

	int    covtype;

	int    fail_exit;

};



#define MAX(A,B)   (((A) > (B)) ? (A) : (B) )
#define MIN(A,B)   (((A) < (B)) ? (A) : (B) ) 


#define mix(a , b , c) \
{ \
	a -= b; a -= c; a ^= (c>>13); \
	b -= c; b -= a; b ^= (a<<8); \
	c -= a; c -= b; c ^= (b>>13); \
	a -= b; a -= c; a ^= (c>>12);  \
	b -= c; b -= a; b ^= (a<<16); \
	c -= a; c -= b; c ^= (b>>5); \
	a -= b; a -= c; a ^= (c>>3);  \
	b -= c; b -= a; b ^= (a<<10); \
	c -= a; c -= b; c ^= (b>>15); \
}

#define zigstep 128 /* Number of Ziggurat'Steps */


#define znew   (z = 36969*(z&65535) + (z>>16) )
#define wnew   (w = 18000*(w&65535) + (w>>16) )
#define MWC    ((znew<<16) + wnew )
#define SHR3   ( jsr ^= (jsr<<17), jsr ^= (jsr>>13), jsr ^= (jsr<<5) )
#define CONG   (jcong = 69069*jcong + 1234567)
#define KISS   ((MWC^CONG) + SHR3)


#ifdef ranKISS

#define randint KISS

#define rand() (randint*2.328306e-10)

#endif



#ifdef ranSHR3

#define randint SHR3

#define rand() (0.5 + (signed)randint*2.328306e-10)

#endif

/*--------------------------------------------------------------- */


#ifdef __x86_64__
    typedef int UL;
#else
    typedef unsigned long UL;
#endif


/*--------------------------------------------------------------- */


static UL jsrseed = 31340134 , jsr;

#ifdef ranKISS

static UL z=362436069, w=521288629, jcong=380116160;

#endif


static UL iz , kn[zigstep];
static long hz;
static double wn[zigstep] , fn[zigstep];



/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/

static double ng;


/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/


void randini(void);
void randnini(void);
double nfix(void);
double randn(void);
void lubksb(double *, int , int *, double *);
int ludcmp(double *, int , int *, double * , double *);
double inv(double * , double * , double *  , double * , double * , int * , int , int );
double maxeigen(double * , int , double * , double *);
void copy_GMM(double * , double * , double * , double * , double * , double * , int , int );
void   RZt(double * , double * , double * , int , int , int );
void qsindex( double * , int * , int , int  );

void Jsplit_GMM(double * , double * ,
				double * ,
				int , int , int , 
				double * , int *);

void Jmerge_GMM(double * ,
				double * ,
				int , int , int , 
				double *, int * , int *);

void responsibility_GMM(double * , double * , double * , double * , 
						double * , double * , 
						int , int  , int * , int  , int , 
						double * , 
						double * , double * , double * , double * , int *);

void mergeinit_GMM(double * , double * , double * , int , int , 
				   int);

void splitinit_GMM(double * , double * , double * , int , int , double , 
				   int , double * , double * , double *);



void truncate_GMM(double * , double * , double * , int , int , int );


void extend_GMM(double * , double * , double * , int , int);




void em_mvgm(double * , double * , double * , double * , double , double , int * , int , int , 
			 double * , int * ,
			 int , int , int , int , 
			 double *, double * , 
			 int *, double * , double *, 
			 double * , double * , 
			 double * , double * , double * , double * , int *);



void fsmem_mvgm(double * ,  double *  , double * , double * , struct opts  ,
				int , int , 
				double * , double * , double * , double * ,
				int * , int *);


/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/


void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )

{


	/* inputs */

	double *Z , *M0 , *S0 , *P0;

	int numdimsZ  , numdimsM0 , numdimsS0 , numdimsP0;

	const int *dimsZ , *dimsM0 , *dimsS0 , *dimsP0;

	mxArray *mxtemp;


	/* outputs */

	double *loglfinal , *Mfinal , *Sfinal , *Pfinal;

	int numdimsloglfinal , numdimsMfinal , numdimsSfinal , numdimsPfinal;

	int *dimsloglfinal  , *dimsMfinal , *dimsSfinal , *dimsPfinal;


	struct opts options = {1 , 5 , 20 , 100 , 100 , 100 , 1e-6 , 1e-6 , 0.01 , 5 , 1 , 5 , 0 , 2};


	int i , l ,  d , N , d2 , idx, ld;

	int Kfinal=5 , Tfinal=1 , Kini=5 , Kmax , tempint;

	double temp;

	double *logl , *M , *S , *P;     /* Local GMM parameters */

	double  *tmp , *moy , *var;





	/*---------------------------------------------------------------*/
	/*---------------------- PARSE INPUTS----------------------------*/	
	/*---------------------------------------------------------------*/

	if( (nrhs < 1) )

	{

		mexPrintf("                                                                                                                                \n");
		mexPrintf("fsmem_mvgm : Free Split and Merge Expectation-Maximization algorithm for Multivariate Gaussian Mixtures                         \n");
		mexPrintf("                                                                                                                                \n");
		mexPrintf("Usage                                                                                                                           \n");        
		mexPrintf(" -------                                                                                                                        \n");
		mexPrintf("                                                                                                                                \n");
		mexPrintf("[logl , M , S , P]  = fsmem_mvgm(Z , [M0] , [S0] , [P0] , [option]);                                                            \n");
		mexPrintf("                                                                                                                                \n");
		mexPrintf("                                                                                                                                \n");
		mexPrintf("Inputs                                                                                                                          \n");
		mexPrintf("-------                                                                                                                         \n");
		mexPrintf("                                                                                                                                \n");
		mexPrintf("Z             Measurements (d x N)                                                                                              \n");
		mexPrintf("M0            Initial mean vector. M0 can be (d x 1 x K)        (default [Kini random elements from Z])                         \n");
		mexPrintf("S0            Initial covariance matrix. S0 can be (d x d x K)  (default [cov(Z)/40])                                           \n");
		mexPrintf("P0            Initial mixture probablities (1 x 1 x K) :        (default [1/Kini])                                              \n");         
		mexPrintf("options                                                                                                                         \n");
		mexPrintf("	             Kmin              Minimum number of compounds (default [1])                                                       \n");
		mexPrintf("	             Kini              Initial number of compounds (default [5])                                                       \n");
		mexPrintf("	             Kmax              Maximum number of compounds (default [15])                                                      \n");
		mexPrintf("                     maxite_fsmem      Number of maximum iteration for the main loop of the fsmem (default [100])                      \n");
		mexPrintf("	             maxite_fullem     Number of maximum iteration for the full EM inside the main loop (default [100])                \n");
		mexPrintf("	             maxite_partialem  Number of maximum iteration for the partial EM inside the main loop (default [100])             \n");
		mexPrintf("	             epsi_fullem       Tolerance in loglikelihood improuvement of the Full EM (default [1e-6])                         \n");
		mexPrintf("	             epsi_partialem    Tolerance in loglikelihood improuvement of the Partial EM (default [1e-6])                      \n");
		mexPrintf("	             lambda            Covariance regularization parameter (default [0.01])                                            \n");
		mexPrintf("	             maxcands_split    Maximum number of split candidate (default [5])                                                 \n");
		mexPrintf("	             splitinit_epsi    Split Initialisation parameter for the mean of splitted cluster (default [1])                   \n");
		mexPrintf("	             maxcands_merge    Maximum number of merge candidate (default [5])                                                 \n");
		mexPrintf("	             covtype           Covariance type : 0 = full , 1 = elliptical , 2 = spherical (default [0])                       \n");
		mexPrintf("	             fail_exit         Number of tentatives of split/merge operations before exit. If fail_exit = 0, then FSMEM = EM   \n");
		mexPrintf("                                                                                                                                \n");
		mexPrintf("Ouputs                                                                                                                          \n");
		mexPrintf("-------                                                                                                                         \n");
		mexPrintf("                                                                                                                                \n");
		mexPrintf("logl          Loglikelihood & information matrix (6 x Tfinal) where Tfinal is the total number of fsmem iterations              \n");
		mexPrintf("              logl(1,i) loglikelihood at iteration i=1,...,Tfinal                                                               \n");
		mexPrintf("              logl(2,i) loglikelihood of the partial EM at iteration i=1,...,Tfinal                                             \n");
		mexPrintf("              logl(3,i) Number of iterations of the partial EM                                                                  \n");
		mexPrintf("              logl(4,i) loglikelihood of the full EM at iteration i=1,...,Tfinal                                                \n");
		mexPrintf("              logl(5,i) Number of iterations of the full EM                                                                     \n");
		mexPrintf("              logl(6,i) Flag of the split/merge moves : {-1,-0.5} success/failed merge move, {1,0.5} sucess/failed split move   \n");
		mexPrintf("M             Estimated mean vector (d x 1 x Kest), where Kest is the number of estimated coupounds                             \n");
		mexPrintf("S             Estimated covariance vector (d x d x Kest)                                                                        \n");
		mexPrintf("P             Estimated initial probabilities (1 x 1 x Kest)                                                                    \n");
		mexPrintf("                                                                                                                                \n");
		mexErrMsgTxt("At least 1 input for fsmem_mvgm\n");
	}




	/* Inputs 1 */


	Z          = mxGetPr(prhs[0]);

	numdimsZ   = mxGetNumberOfDimensions(prhs[0]);

	dimsZ      = mxGetDimensions(prhs[0]);

	d          = dimsZ[0];

	N          = dimsZ[1];

	d2         = d*d;

	ng         = 1.0/pow(2*M_PI , d/2);


	/* Inputs 2 */


	if ( (nrhs> 1) && !mxIsEmpty(prhs[1]))
	{


		M0         = mxGetPr(prhs[1]);

		numdimsM0  = mxGetNumberOfDimensions(prhs[1]);

		dimsM0     = mxGetDimensions(prhs[1]);

		if (  (dimsM0[1] != 1)  )
		{

			mexErrMsgTxt("M0 must be (d x 1 x Kini)");			 

		}

		if(numdimsM0 == 3)
		{

			Kini           = dimsM0[2];

			options.Kini   = Kini;   /* If Initial GMM is inputed, force options's field update */

		}
		else
		{

			options.Kini   = 1;   /* If Initial GMM is inputed, force options's field update */

		}


	}


	/* Inputs 3 */


	if ( (nrhs> 2) && !mxIsEmpty(prhs[2]))
	{


		S0         = mxGetPr(prhs[2]);

		numdimsS0  = mxGetNumberOfDimensions(prhs[2]);

		dimsS0     = mxGetDimensions(prhs[2]);

		if ( (dimsS0[0] != d) && (dimsS0[1] != d) && (dimsS0[2] != Kini)  )
		{

			mexErrMsgTxt("S0 must be (d x d x Kini)");			 

		}	

	}


	/* Inputs 4 */


	if ( (nrhs> 3) && !mxIsEmpty(prhs[3]))
	{


		P0         = mxGetPr(prhs[3]);

		numdimsP0  = mxGetNumberOfDimensions(prhs[3]);

		dimsP0     = mxGetDimensions(prhs[3]);

		if (  dimsP0[0] != 1 && dimsP0[1] != 1 && dimsP0[2] != Kini )
		{

			mexErrMsgTxt("P must be (1 x 1 x Kini ) ");			 

		}

	}	


	/* Inputs 5 */


	if((nrhs > 4)  && !mxIsEmpty(prhs[4]))

	{


		mxtemp                            = mxGetField(prhs[4] , 0, "Kmin");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < 1)
			{

				mexPrintf("Kmin >= 1, force default value, Kmin = 1");	

				options.Kmin              = 1;

			}
			else
			{

				options.Kmin              = tempint;

			}

		}


		mxtemp                            = mxGetField(prhs[4] , 0 , "Kini");

		if((mxtemp != NULL) && mxIsEmpty(prhs[1]))
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < options.Kmin)
			{

				mexPrintf("Kini >= Kmin, force default value");	

				options.Kini              = MAX(5 , options.Kmin);

			}
			else
			{

				options.Kini              = tempint;

			}

		}



		mxtemp                            = mxGetField(prhs[4] , 0, "Kmax");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if (tempint < options.Kini)
			{

				mexPrintf("Kmax >= Kini, force default value");	

				options.Kmax              = options.Kini;

			}

			else
			{

				options.Kmax               = tempint;

			}

		}

		mxtemp                            = mxGetField(prhs[4] , 0, "maxite_fsmem");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < 0)
			{

				mexPrintf("options.maxite_fsmem > 0");

				options.maxite_fsmem      = 1000;

			}
			else
			{

				options.maxite_fsmem      = tempint;

			}

		}


		mxtemp                            = mxGetField(prhs[4] , 0, "maxite_fullem");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < 0)
			{

				mexPrintf("options.maxite_fullem > 0");	

				options.maxite_fullem     = 1000;

			}

			else
			{

				options.maxite_fullem         = tempint;

			}

		}

		mxtemp                            = mxGetField(prhs[4] , 0, "maxite_partialem");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < 0)
			{

				mexPrintf("options.maxite_partialem > 0");	

				options.maxite_partialem  = 1000;

			}

			else
			{

				options.maxite_partialem  = tempint;

			}

		}


		mxtemp                            = mxGetField(prhs[4] , 0, "epsi_fullem");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			if(*tmp < 0)
			{

				mexPrintf("options.epsi_fullem >= 0");

				options.epsi_fullem       = 1e-6;



			}
			else
			{

				options.epsi_fullem       = tmp[0];

			}

		}


		mxtemp                            = mxGetField(prhs[4] , 0, "epsi_partialem");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			if(*tmp < 0)
			{

				mexPrintf("options.epsi_partialem >= 0");

				options.epsi_partialem   = 1e-6;




			}
			else

			{

				options.epsi_partialem    = tmp[0];

			}

		}


		mxtemp                            = mxGetField(prhs[4] , 0, "lambda");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			if(*tmp < 0)
			{

				mexPrintf("options.lambda >= 0");	

				options.lambda            = 0.01;
			}

			else
			{

				options.lambda            = tmp[0];

			}

		}


		mxtemp                            = mxGetField(prhs[4] , 0, "maxcands_split");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < 0)				
			{

				mexPrintf("options.maxcands_split >= 0");	

				options.maxcands_split    = 5;

			}
			else

			{
				options.maxcands_split    = tempint;


			}

		}


		mxtemp                            = mxGetField(prhs[4] , 0, "splitinit_epsi");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);


			if(*tmp < 0)

			{
				mexPrintf("options.splitinit_epsi >= 0");	

				options.splitinit_epsi        = 1.0;

			}

			else
			{

				options.splitinit_epsi        = tmp[0];

			}

		}

		mxtemp                            = mxGetField(prhs[4] , 0, "maxcands_merge");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < 0)				
			{

				mexPrintf("options.maxcands_merge >= 0");	

				options.maxcands_merge    = 5;


			}
			else

			{

				options.maxcands_merge    = tempint;

			}


		}


		mxtemp                            = mxGetField(prhs[4] , 0, "covtype");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if( (tempint <0) || (tempint > 2) )
			{

				mexPrintf("options.covtype = {0,1,2} for {full, elliptical, spherical}, force to 0");

				options.covtype           = 0;

			}

			options.covtype               = tempint;

		}


		mxtemp                            = mxGetField(prhs[4] , 0, "fail_exit");

		if(mxtemp != NULL)
		{

			tmp                           = mxGetPr(mxtemp);

			tempint                       = (int) tmp[0];

			if(tempint < 0)
			{

				mexPrintf("options.fail_exit must be > 0, force to 2");

				options.fail_exit         = 2;


			}
			else
			{

				options.fail_exit         = tempint;

			}

		}

	}


	Kini                                  = options.Kini;

	Kmax                                  = options.Kmax;


	/* Random Generator Initialization */

	randini();

	randnini();



	if ( (nrhs < 2) || mxIsEmpty(prhs[1]))

	{

		/*		Kini        = options.Kini; */

		M0          = (double *)mxMalloc(d*Kini*sizeof(double));

		for( l = 0 ; l < Kini ; l++)
		{

			idx         = ( (int)(ceil(N*rand()))-1)*d;

			ld          = l*d;

			for (i = 0 ; i < d ; i++)
			{

				M0[i + ld]  = Z[i + idx];

			}

		}

	}


	if ( (nrhs < 3) || mxIsEmpty(prhs[2]))

	{

		S0          = (double *)mxMalloc(d2*Kini*sizeof(double));

		moy         = (double *)mxMalloc(d*sizeof(double));

		var         = (double *)mxMalloc(d*sizeof(double));



		for (i = 0 ; i < d ; i++)
		{

			moy[i]  = 0.0;

			var[i]  = 0.0;

		}

		for (i = 0 ; i < d2*Kini ; i++)
		{

			S0[i]   = 0.0;

		}

		/* Mean */

		for (l = 0 ; l < N ; l++)
		{

			ld    = l*d;

			for(i = 0 ; i < d ; i++)
			{

				moy[i] += Z[i + ld];

			}

		}


		temp        = 1.0/N;

		for (i = 0 ; i < d ; i++)
		{

			moy[i]  *= temp;

		}

		/* Variance  */

		for (l = 0 ; l < N ; l++)
		{

			ld    = l*d;

			for(i = 0 ; i < d ; i++)
			{

				temp        = (Z[i + ld] - moy[i]);

				var[i]     += temp*temp;

			}

		}


		temp        = 1.0/(40*(N-1));

		for (i = 0 ; i < d ; i++)
		{

			var[i]  *= temp;

		}


		for( l = 0 ; l < Kini ; l++)
		{

			ld          = l*d2;

			for (i = 0 ; i < d ; i++)
			{

				S0[i*(d+1) + ld]  = var[i];

			}

		}

	}



	if ( (nrhs < 4) || mxIsEmpty(prhs[3]))

	{


		P0          = (double *)mxMalloc(Kini*sizeof(double));


		temp        = 1.0/Kini;

		for (i = 0 ; i < Kini ; i++)
		{

			P0[i]  = temp;

		}

	}



	logl       = (double *)mxMalloc(6*Kmax*sizeof(double));

	M          = (double *)mxMalloc(d*Kmax*sizeof(double));

	S          = (double *)mxMalloc(d2*Kmax*sizeof(double));

	P          = (double *)mxMalloc(Kmax*sizeof(double));


	/*---------------------------------------------------------------*/
	/*------------------------ MAIN CALL ----------------------------*/
	/*---------------------------------------------------------------*/



	fsmem_mvgm(Z ,  M0  , S0 , P0 , options ,
		d , N , 
		logl , M , S , P , 
		&Kfinal , &Tfinal);


	/*---------------------------------------------------------------*/
	/*---------------------- PARSE OUTPUT ---------------------------*/	
	/*---------------------------------------------------------------*/

	numdimsloglfinal    = 2;

	dimsloglfinal       = (int *)mxMalloc(numdimsloglfinal*sizeof(int));

	dimsloglfinal[0]    = 6;

	dimsloglfinal[1]    = Tfinal;




	numdimsMfinal       = 3;

	dimsMfinal          = (int *)mxMalloc(numdimsMfinal*sizeof(int));

	dimsMfinal[0]       = d;

	dimsMfinal[1]       = 1;

	dimsMfinal[2]       = Kfinal;



	numdimsSfinal       = 3;

	dimsSfinal          = (int *)mxMalloc(numdimsSfinal*sizeof(int));

	dimsSfinal[0]       = d;

	dimsSfinal[1]       = d;

	dimsSfinal[2]       = Kfinal;



	numdimsPfinal       = 3;

	dimsPfinal          = (int *)mxMalloc(numdimsPfinal*sizeof(int));

	dimsPfinal[0]       = 1;

	dimsPfinal[1]       = 1;

	dimsPfinal[2]       = Kfinal;




	plhs[0]             = mxCreateNumericArray(numdimsloglfinal , dimsloglfinal , mxDOUBLE_CLASS, mxREAL);

	loglfinal           = mxGetPr(plhs[0]);


	plhs[1]             = mxCreateNumericArray(numdimsMfinal , dimsMfinal , mxDOUBLE_CLASS, mxREAL);

	Mfinal              = mxGetPr(plhs[1]);


	plhs[2]             = mxCreateNumericArray(numdimsSfinal , dimsSfinal , mxDOUBLE_CLASS, mxREAL);

	Sfinal              = mxGetPr(plhs[2]);


	plhs[3]             = mxCreateNumericArray(numdimsPfinal , dimsPfinal , mxDOUBLE_CLASS, mxREAL);

	Pfinal              = mxGetPr(plhs[3]);



	/* Copy FSMEM results */


	for (i = 0 ; i < 6*Tfinal ; i++)
	{

		loglfinal[i]  = logl[i];

	}


	copy_GMM(Mfinal , Sfinal , Pfinal , M , S , P , d , Kfinal);



	/*---------------------------------------------------------------*/
	/*------------------------ FREE MEMORY --------------------------*/
	/*---------------------------------------------------------------*/



	mxFree(logl);

	mxFree(M);

	mxFree(S);

	mxFree(P);


	if ( (nrhs < 2) || mxIsEmpty(prhs[1]))

	{

		mxFree(M0);

	}

	if ( (nrhs < 3) || mxIsEmpty(prhs[2]))

	{

		mxFree(S0);

		mxFree(moy);

		mxFree(var);

	}


	if ( (nrhs < 4) || mxIsEmpty(prhs[3]))

	{

		mxFree(P0);

	}

}

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/



void fsmem_mvgm(double *Z ,  double *M0  , double *S0 , double *P0 , struct opts options ,
				int d , int N , 
				double *logl , double *M , double *S , double *P ,
				int *K_final , int *T_final)

{


	int i , t , c , d2 = d*d , tlog;

	int maxit_fsmem = options.maxite_fsmem, Kmin = options.Kmin, Kini = options.Kini , Kmax = options.Kmax , Kpair = Kmax*(Kmax-1)/2;

	int maxit_partial = options.maxite_partialem , maxit_full = options.maxite_fullem;

	int Kcurrent = Kini , Knew = -1;

	int next = 1 , maxcands_merge = options.maxcands_merge , nbmerge_cands , maxcands_split = options.maxcands_split , nbsplit_cands;

	int covtype = options.covtype , fail = 0, fail_exit = options.fail_exit;

	int m1 , m2 , s , s2;

	int T_current , T_partial , T_new;

	double lambda = options.lambda , epsi_fullem = options.epsi_fullem , epsi_partialem = options.epsi_partialem , splitinit_epsi = options.splitinit_epsi;

	double likely_current , likely_partial , likely_new ;


	int *indx , *index , *indexsort , *idxcluster_partial , *tempJmergeindex;


	double *R , *G ;

	double *dZ , *invSkdZ;

	double *px , *primalR;

	double *invSk , *tempSk, *vect , *vv;

	double *Jsplit , *Jmerge , *tempJsplit , *tempJmerge;

	double *Mnew , *Snew , *Pnew ; 



	R                  = (double *)mxMalloc(N*Kmax*sizeof(double));

	G                  = (double *)mxMalloc(N*Kmax*sizeof(double));

	dZ                 = (double *)mxMalloc(d*N*sizeof(double));

	invSkdZ            = (double *)mxMalloc(d*N*sizeof(double));

	px                 = (double *)mxMalloc(N*sizeof(double));

	primalR            = (double *)mxMalloc(N*sizeof(double));


	invSk              = (double *)mxMalloc(d2*sizeof(double));

	tempSk             = (double *)mxMalloc(d2*sizeof(double));

	vect               = (double *)mxMalloc(d*sizeof(double));

	vv                 = (double *)mxMalloc(d*sizeof(double));	

	Jsplit             = (double *)mxMalloc(Kmax*2*sizeof(double));

	tempJsplit         = (double *)mxMalloc(Kmax*sizeof(double));

	Jmerge             = (double *)mxMalloc(Kpair*3*sizeof(double));

	tempJmerge         = (double *)mxMalloc(Kpair*sizeof(double));

	Mnew               = (double *)mxMalloc(d*Kmax*sizeof(double));

	Snew               = (double *)mxMalloc(d2*Kmax*sizeof(double));

	Pnew               = (double *)mxMalloc(Kmax*sizeof(double));



	index              = (int *)mxMalloc(Kmax*sizeof(int));

	idxcluster_partial = (int *)mxMalloc(Kmax*sizeof(int));

	indexsort          = (int *)mxMalloc(Kpair*sizeof(int));

	tempJmergeindex    = (int *)mxMalloc(Kpair*2*sizeof(int));

	indx               = (int *)mxMalloc(d*sizeof(int));



	/* Copy Initial parameters */


	copy_GMM(M , S , P , M0 , S0 , P0 , d , Kcurrent);



	/* First Run of Full-EM */

	for (i = 0 ; i < Kcurrent ; i++)
	{

		idxcluster_partial[i] = i;

	}

	em_mvgm(Z , M , S , P , epsi_fullem , lambda , idxcluster_partial , Kcurrent , maxit_full , 
		&likely_current , &T_current ,
		d , Kcurrent , N , covtype , 
		R , G , 
		index , dZ , invSkdZ , 
		px , primalR , 
		invSk , tempSk , vect , vv , indx);


	logl[0]     = likely_current;

	logl[1]     = 0.0;

	logl[2]     = 0.0;

	logl[3]     = likely_current;

	logl[4]     = (double) T_current;

	logl[5]     = 0.0;

	(*T_final)  = 1;

	(*K_final)  = Kcurrent;





	if(fail_exit > 0)
	{

		/* Main Loop of FSMEM */

		if(Kcurrent == 1)
		{

			next = -1;

		}

		tlog             = 6;

		for (t = 0 ; t < maxit_fsmem - 1 ; t++)

		{

			if(fail == 0)     /* A previous S & M move was successfull, update R and J matrices */
			{

				for (i = 0 ; i < Kcurrent ; i++)
				{

					index[i] = i;

				}


				responsibility_GMM(Z , M , S , P , 
					R , G , 
					d , N , index , Kcurrent , covtype ,
					px , 
					invSk , tempSk , vect , vv , indx);


				Jsplit_GMM(R , G ,
					Jsplit ,
					N , Kcurrent , Kmax, 
					tempJsplit , indexsort);


				Jmerge_GMM(R ,
					Jmerge ,
					N , Kcurrent , Kpair , 
					tempJmerge , tempJmergeindex , indexsort);	

			}

			if((next > 0) && (Kcurrent > Kmin - 1)) /* Merge Move	*/	
			{


				Knew           = Kcurrent - 1;    /* New hypothetic number of clusters */

				nbmerge_cands  = Kcurrent*(Kcurrent-1)/2;



				if(nbmerge_cands > maxcands_merge)

				{

					nbmerge_cands = maxcands_merge;

				}

				likely_new   = likely_ini;

				for(c =  0 ; c < nbmerge_cands ; c++)
				{


					copy_GMM(Mnew , Snew , Pnew , M , S , P , d , Kcurrent);


					m1    = (int)(Jmerge[c + Kpair]);

					m2    = (int)(Jmerge[c + 2*Kpair]);

					mergeinit_GMM(Mnew , Snew , Pnew , m1 , m2 , d);

					truncate_GMM(Mnew , Snew , Pnew , m2 , Kcurrent , d);



					/* Partial EM on merged cluster m1 */

					idxcluster_partial[0] = m1;

					em_mvgm(Z , Mnew , Snew , Pnew , epsi_partialem , lambda , idxcluster_partial , 1 , maxit_partial , 
						&likely_partial , &T_partial ,
						d , Knew , N , covtype , 
						R , G , 
						index , dZ , invSkdZ , 
						px , primalR , 
						invSk , tempSk , vect , vv , indx);


					/* Full EM */


					for (i = 0 ; i < Knew ; i++)
					{

						idxcluster_partial[i] = i;

					}

					em_mvgm(Z , Mnew , Snew , Pnew , epsi_fullem , lambda , idxcluster_partial , Knew , maxit_full , 
						&likely_new , &T_new ,
						d , Knew , N , covtype , 
						R , G , 
						index , dZ , invSkdZ  , 
						px , primalR , 
						invSk , tempSk , vect , vv , indx);

					if(likely_new > likely_current)
					{

						break;

					}

				}

				if(likely_new > likely_current)
				{

					copy_GMM(M , S , P , Mnew , Snew , Pnew , d , Knew);

					Kcurrent       = Knew;

					likely_current = likely_new;

					fail           = 0;

					logl[0 + tlog] = likely_new;

					logl[1 + tlog] = likely_partial;

					logl[2 + tlog] = (double) T_partial;

					logl[3 + tlog] = likely_new;

					logl[4 + tlog] = (double) T_new;

					logl[5 + tlog] = -1.0;

					if(Kcurrent == Kmin)
					{

						next       = 1;
					}

				}

				else

				{

					fail          += 1;

					next           = -1;

					logl[0 + tlog] = likely_current;

					logl[1 + tlog] = likely_partial;

					logl[2 + tlog] = (double) T_partial;

					logl[3 + tlog] = likely_new;

					logl[4 + tlog] = (double) T_new;

					logl[5 + tlog] = -0.5;

				}		
			}

			if((next < 0) && (Kcurrent < Kmax))   /* Split Move  */

			{


				nbsplit_cands  = Kcurrent;

				if(nbsplit_cands > maxcands_split)

				{

					nbsplit_cands = maxcands_split;

				}

				Knew          = Kcurrent + 1;

				likely_new    = likely_ini;


				for(c =  0 ; c < nbsplit_cands ; c++)
				{


					copy_GMM(Mnew , Snew , Pnew , M , S , P , d , Kcurrent);


					s                     = (int) (Jsplit[c + Kmax]);

					s2                    = Kcurrent; /* Knew - 1, an index */

					extend_GMM(Mnew , Snew , Pnew , d , s2);

					splitinit_GMM(Mnew , Snew , Pnew , s , s2 , splitinit_epsi , 
						d , vect , tempSk , vv);

					/* Partial EM on Splitted clusters s & s2 */

					idxcluster_partial[0] = s;

					idxcluster_partial[1] = s2;

					em_mvgm(Z , Mnew , Snew , Pnew , epsi_partialem , lambda , idxcluster_partial , 2 , maxit_partial , 
						&likely_partial , &T_partial ,
						d , Knew , N , covtype , 
						R , G , 
						index , dZ , invSkdZ , 
						px , primalR , 
						invSk , tempSk , vect , vv , indx);


					/* Full EM */ 


					for (i = 0 ; i < Knew ; i++)
					{

						idxcluster_partial[i] = i;

					}

					em_mvgm(Z , Mnew , Snew , Pnew , epsi_fullem , lambda , idxcluster_partial , Knew , maxit_full , 
						&likely_new , &T_new ,
						d , Knew , N , covtype , 
						R , G , 
						index , dZ , invSkdZ  , 
						px , primalR , 
						invSk , tempSk , vect , vv , indx);

					if(likely_new > likely_current)
					{

						break;

					}		

				}

				if(likely_new > likely_current)
				{

					copy_GMM(M , S , P , Mnew , Snew , Pnew , d , Knew);

					Kcurrent       = Knew;

					likely_current = likely_new;

					fail           = 0;

					logl[0 + tlog] = likely_new;

					logl[1 + tlog] = likely_partial;

					logl[2 + tlog] = (double) T_partial;

					logl[3 + tlog] = likely_new;

					logl[4 + tlog] = (double) T_new;

					logl[5 + tlog] = 1.0;

					if(Kcurrent == Kmax)
					{

						next       = 1;
					}

				}

				else

				{

					fail          += 1;

					next           = 1;

					logl[0 + tlog] = likely_current;

					logl[1 + tlog] = likely_partial;

					logl[2 + tlog] = (double) T_partial;

					logl[3 + tlog] = likely_new;

					logl[4 + tlog] = (double) T_new;

					logl[5 + tlog] = 0.5;
				}			

			}


			else if((fail > 0) && (Kcurrent == Kmin) || (Kcurrent == Kmax))
			{

				fail = fail_exit;

			}

			if(fail >= fail_exit)
			{

				(*T_final)  = (t + 1);

				(*K_final)  = Kcurrent;

				break;


			}

			tlog    += 6;

		}

	}


	mxFree(R);

	mxFree(G);

	mxFree(px);

	mxFree(primalR);

	mxFree(dZ);

	mxFree(invSkdZ);

	mxFree(invSk);

	mxFree(tempSk);

	mxFree(vect);

	mxFree(vv);

	mxFree(indx);

	mxFree(index);

	mxFree(idxcluster_partial);

	mxFree(indexsort);

	mxFree(tempJsplit);

	mxFree(tempJmerge);

	mxFree(tempJmergeindex);

	mxFree(Jsplit);

	mxFree(Jmerge);

	mxFree(Mnew);

	mxFree(Snew);	

	mxFree(Pnew);

}



/* ---------------------------------------------------------------------------------------------------------------------------- */


void   copy_GMM(double *Mdest , double *Sdest , double *Pdest , 
				double *Msrc  , double *Ssrc  , double *Psrc , 
				int d , int K)

{

	int i , l , ld  , d2 = d*d;


	for (l = 0 ; l < K ; l++)
	{

		ld              = l*d;

		for(i = 0 + ld ; i < d + ld ; i++)
		{

			Mdest[i]    = Msrc[i];

		}

		ld              = l*d2;

		for(i = 0 + ld ; i < d2 + ld; i++)
		{

			Sdest[i]    = Ssrc[i];

		}

		Pdest[l]        = Psrc[l];		


	}

}


/* ---------------------------------------------------------------------------------------------------------------------------- */

void responsibility_GMM(double *Z , double *M , double *S , double *P , 
						double *R , double *G , 
						int d , int N , int *index , int K , int covtype , 
						double *px , 
						double *invSk , double *tempSk , double *vect , double *vv , int *indx)

{


	int i , j , n , k , d2 = d*d , d1 = d+1, nd , kk , kd , kd2 , kN , jd , id;

	double invdet , cteg , g , pg;

	register double sum , sumg , temp , prod;

	for (n = 0 ; n < N ; n++)
	{	

		px[n]      = tiny;

	}

	if(covtype == 0)   /* Full Cov */
	{

		for (k = 0 ; k < K ; k++)
		{


			kk          = index[k];

			kd          = kk*d;

			kd2         = kk*d2;

			kN          = kk*N;


			/* inv(S(: , : , k)) & det(S(: , : , k)) */

			invdet      = inv(S , invSk , 
				tempSk , vect , vv , indx , d , kd2);

			cteg        = ng*sqrt(fabs(invdet));


			/* dZ = Z - mean(Z)  (d x N) */

			nd             = 0;

			for (n = 0 ; n < N ; n++)
			{

				for (i = 0 ; i < d ; i++)
				{

					vect[i]  = (Z[i + nd] - M[i + kd]);  /* dZ */

				}

				sumg  = 0.0;

				jd    = 0;

				for(j = 0 ; j < d ; j++)
				{

					sum  = 0.0;

					for(i = 0 ; i < d ; i++)

					{
						sum    += invSk[i + jd]*vect[i];    /* Since invSk is symetric too ... */

					}


					sumg  += sum*vect[j];    /* inv(Sk)*dZ */

					jd    += d;

				}

				g         = cteg*exp(-0.5*sumg);

				G[n + kN] = g;

				pg        = P[kk]*g;

				R[n + kN] = pg;

				px[n]    += pg;

				nd       += d;

			}

		}
	}
	if(covtype == 1)  /* Elliptical Cov */
	{
		for (k = 0 ; k < K ; k++)
		{


			kk          = index[k];

			kd          = kk*d;

			kd2         = kk*d2;

			kN          = kk*N;


			/* inv(S(: , : , k)) & det(S(: , : , k)) */

			prod        = 1.0;

			for(i = 0 ; i < d ; i++)
			{

				temp          = 1.0/S[i*d1 + kd2];

				invSk[i*d1]   = temp;

				prod          *= temp;

			}

			cteg        = ng*sqrt(fabs(prod));


			nd             = 0;

			for (n = 0 ; n < N ; n++)
			{

				sumg       = 0.0;

				id         = 0;

				for (i = 0 ; i < d ; i++)
				{

					temp     = (Z[i + nd] - M[i + kd]);  /* dZ = Z - mean(Z)  (d x 1)  */

					sumg    += invSk[id]*temp*temp;

					id      += d1;

				}


				g         = cteg*exp(-0.5*sumg);

				G[n + kN] = g;

				pg        = P[kk]*g;

				R[n + kN] = pg;

				px[n]    += pg;

				nd       += d; 

			}

		}

	}

	if(covtype == 2)  /* Spherical Cov */
	{
		for (k = 0 ; k < K ; k++)
		{


			kk          = index[k];

			kd          = kk*d;

			kd2         = kk*d2;

			kN          = kk*N;


			/* inv(S(: , : , k)) & det(S(: , : , k)) */


			prod        = 1.0/S[kd2];

			cteg        = ng*sqrt(fabs(prod));

			nd          = 0;

			for (n = 0 ; n < N ; n++)
			{

				sumg       = 0.0;

				for (i = 0 ; i < d ; i++)
				{

					temp     = (Z[i + nd] - M[i + kd]); 	/* dZ = Z - mean(Z) */

					sumg    += (temp*temp);

				}


				g         = cteg*exp(-0.5*sumg*prod);

				G[n + kN] = g;

				pg        = P[kk]*g;

				R[n + kN] = pg;

				px[n]    += pg;

				nd       += d;
			}

		}

	}



	for (k = 0 ; k < K ; k++)
	{

		kN          = index[k]*N;

		for (n = 0 ; n < N ; n++)
		{

			R[n + kN] /= px[n];

		}

	}

}



/* ---------------------------------------------------------------------------------------------------------------------------- */

void Jsplit_GMM(double *R , double *G ,
				double *Jsplit ,
				int N , int K , int Kmax , 
				double *tempJsplit , int *indexsort)


{

	/*	
	Maximum Local Kullback Divergence Criterion
	Jsplit (Kmax x 2), K <= Kmax

	*/


	int k , n , kN;

	double sumRk , temp , f , sumJ;


	for (k = 0 ; k < K ; k++)
	{

		kN           = k*N;

		indexsort[k] = k;

		sumRk        = fudge;

		for(n = 0 + kN ; n < N + kN; n++)
		{

			sumRk   += R[n];

		}

		temp         = 1.0/sumRk;

		sumJ         = 0.0;

		for(n = 0 + kN ; n < N + kN ; n++)
		{

			f      = R[n]*temp;

			if(f > tiny)

			{

				sumJ     += f*log(f/(G[n] + tiny));   /* Maximum Local Kullback Divergence Criterion */

			}

		}

		tempJsplit[k]   = -sumJ;


	}



	qsindex( tempJsplit , indexsort , 0 , K - 1 );

	for(k = 0 ; k < K ; k++)
	{

		Jsplit[k]        = -tempJsplit[k];

		Jsplit[k + Kmax] = indexsort[k]; 

	}


}

/* ---------------------------------------------------------------------------------------------------------------------------- */

void Jmerge_GMM(double *R ,
				double *Jmerge ,
				int N , int K , int Kpair , 
				double *tempJmerge , int *tempJmergeindex , int *indexsort)


{
	/* 
	Approximation of Maximum Distribution Overlap Criterion
	Jmerge (Kpair x 2), Kpair = Kmax*(Kmax+1)/2, 

	*/


	int idx = 0 , k , l , n , kN , lN , Klimit = K*(K-1)/2;

	register double sumR;

	for (k = 0 ; k < K ; k++)
	{
		kN    = k*N;

		for(l = k+1 ; l < K ; l++)
		{
			lN   = l*N;

			sumR = 0.0;

			for( n = 0 ; n < N ; n++)
			{

				sumR   += R[n + kN]*R[n + lN];     /* Approximation of Maximum Distribution Overlap Criterion */

			}


			tempJmerge[idx]              = -sumR;

			tempJmergeindex[idx]         = k;

			tempJmergeindex[idx + Kpair] = l;

			idx++;

		}
	}


	for (l = 0 ; l < Klimit ; l++)  /* l < idx */
	{

		indexsort[l] = l;

	}

	qsindex( tempJmerge , indexsort , 0 , Klimit - 1 );

	for (l = 0 ; l < Klimit ; l++)  /* l < idx */
	{

		Jmerge[l]              =  -tempJmerge[l];

		Jmerge[l + Kpair]      =  tempJmergeindex[indexsort[l]];

		Jmerge[l + 2*Kpair]    =  tempJmergeindex[indexsort[l] + Kpair];

	}



}


/* ---------------------------------------------------------------------------------------------------------------------------- */



void mergeinit_GMM(double *M , double *S , double *P , int m1 , int m2  , int d )

{

	int i , d2 = d*d ;

	double Pm1 = P[m1] , Pm2 = P[m2] , temp = Pm1 + Pm2 , invtemp = 1.0/temp;


	for(i = 0 ; i < d ; i++)
	{

		M[i + m1*d]   = (Pm1*M[i + m1*d] + Pm2*M[i + m2*d])*invtemp;

	}

	for(i = 0 ; i < d2 ; i++)
	{

		S[i + m1*d2]  = (Pm1*S[i + m1*d2] + Pm2*S[i + m2*d2])*invtemp;

	}


	P[m1]      = temp;

	P[m2]      = 0.0;	

}


/* ---------------------------------------------------------------------------------------------------------------------------- */


void truncate_GMM(double *M , double *S , double *P , int m , int K , int d)

{

	int i , d2 = d*d , idxK = K-1;


	if(m < idxK)
	{

		for(i = 0 ; i < d ; i++)
		{

			M[i + m*d] = M[i + idxK*d];


		}
		for(i = 0 ; i < d2 ; i++)
		{

			S[i + m*d2] = S[i + idxK*d2];

		}

		P[m]           = P[idxK];

	}

}



/* ---------------------------------------------------------------------------------------------------------------------------- */


void extend_GMM(double *M , double *S , double *P , int d , int K)

{

	int i , d1 = (d+1) , d2 = d*d;



	for(i = 0 ; i < d2 ; i++)
	{

		S[i + K*d2]    = 0.0;

	}


	for(i = 0 ; i < d ; i++)
	{

		M[i + K*d]     = 0.0;

		S[i*d1 + K*d2] = 1.0;

	}


	P[K]           = 0.0;

}


/* ---------------------------------------------------------------------------------------------------------------------------- */

void splitinit_GMM(double *M , double *S , double *P , int k , int l , double split_init_epsi , 
				   int d , double *eigenvect , double *Sk , double *vv)

{

	int i , d1 = d+1 , d2 = d*d;

	double eigenval , stdeigen , temp , tmp;

	P[k]    *= 0.5;

	P[l]     = P[k];

	for(i = 0 ; i < d2 ; i++)
	{

		Sk[i] = S[i + k*d2];

	}

	eigenval  = maxeigen(Sk , d , eigenvect , vv);	

	stdeigen  = split_init_epsi*sqrt(eigenval);



	for(i = 0 ; i < d ; i++)
	{

		temp        = M[i + k*d];

		tmp         = eigenvect[i]*stdeigen;

		M[i + k*d]  = temp + tmp*randn();  

		M[i + l*d]  = temp + tmp*randn(); 



	}

	/*	
	for(i = 0 ; i < d2 ; i++)
	{

	S[i + k*d2]    *= 0.5;


	}

	*/


	for(i = 0 ; i < d2 ; i++)
	{

		S[i + k*d2]    = 0.0;


	}


	for(i = 0 ; i < d ;i++)
	{

		S[i*d1 + k*d2] = eigenval;


	}



	for(i = 0 ; i < d2 ; i++)
	{

		S[i + l*d2]    = S[i + k*d2];
	}

}


/* ---------------------------------------------------------------------------------------------------------------------------- */

double maxeigen(double *S , int d , double *eigenvect , double *v)

{
	/* Compute the First Eigenvector & Eigenvalue of the matrix S(d x d) */

	int i , j , t;

	int id;

	double r , temp  = 0.0 , newnorm , oldnorm = huge;

	for(i = 0 ; i < d ; i++)
	{


		v[i]  = rand();

		r     = fabs(v[i]);

		temp += (r*r);

	}

	newnorm   = sqrt(temp);

	temp      = 1.0/newnorm;

	for(i = 0 ; i < d ; i++)
	{

		v[i] *= temp;

	}


	for (t = 0 ; t < maxite_eigen ; t++)
	{


		for (i = 0 ; i < d ; i++)
		{

			temp = 0.0;

			id   = i*d;

			for(j = 0 ; j < d  ; j++)
			{


				temp += S[j + id]*v[j];   /* Since S = S^T */

			}

			eigenvect[i] = temp;

		}

		temp     = 0.0;

		for(i = 0 ; i < d ; i++)
		{


			v[i]  = eigenvect[i];

			r     = fabs(v[i]);

			temp += (r*r);

		}

		newnorm   = sqrt(temp);



		if(fabs(newnorm - oldnorm) < tol_eigen)
		{

			break;
		}

		else
		{

			temp      = 1.0/newnorm;

			for(i = 0 ; i < d ; i++)
			{

				v[i] *= temp;

			}

			oldnorm = newnorm;

		}


	}

	temp      = 1.0/newnorm;

	for(i = 0 ; i < d ; i++)
	{

		v[i]        *= temp;

		eigenvect[i] = v[i];

	}

	return newnorm;	

}



/* ---------------------------------------------------------------------------------------------------------------------------- */



void       em_mvgm(double *Z , double *M , double *S , double *P , double epsi , double lambda , int *idxcluster_partial , int nbcluster_partial , int maxit , 
				   double *likely_final , int *T_final ,
				   int d , int K , int N , int covtype , 
				   double *R , double *G , 
				   int *index , double *dZ , double *invSkdZ , 
				   double *px , double *primalR , 
				   double *invSk , double *tempSk , double *vect , double *vv , int *indx)

{


	int c , i , n , t  , ind , indR , nd , indM , indS , d1 = d+1 , d2 = d*d;

	double sumR , likely_old = likely_ini , likely , temp , invsumR , cteN = 1.0/N , cte_MDL = 0.0;



	if(nbcluster_partial < K)

	{

		for (c = 0 ; c < K ; c++)
		{

			index[c] = c;

		}

		responsibility_GMM(Z , M , S , P , 
			R , G , 
			d , N , index , K , covtype ,
			px , 
			invSk , tempSk , vect , vv , indx);


		for (n = 0 ; n < N ; n++)
		{

			sumR = 0.0;

			for (c = 0 ; c < nbcluster_partial ; c++)
			{

				sumR  += R[n + idxcluster_partial[c]*N];

			}

			primalR[n]  = sumR;

		}

	}

	else

	{

		for (n = 0 ; n < N ; n++)
		{


			primalR[n] = 1.0;

		}


		cte_MDL         = -0.5*log(N)*K*(d1 + d*(d1)/2);

	}


	for (t = 0 ; t < maxit ; t++)
	{


		/* E Step */

		responsibility_GMM(Z , M , S , P , 
			R , G , 
			d , N , idxcluster_partial , nbcluster_partial , covtype , 
			px , 
			invSk , tempSk , vect , vv , indx);


		likely = 0.0;

		for (n = 0 ; n < N ; n++)
		{

			for (c = 0 ; c < nbcluster_partial ; c++)
			{

				R[n + idxcluster_partial[c]*N]  *=primalR[n];

			}

			likely  += log(px[n]);

		}


		if (fabs((likely - likely_old)/likely) < epsi)

		{

			break;
		}

		likely_old = likely;




		/* M Step */



		for (c = 0 ; c < nbcluster_partial ; c++)
		{

			ind   = idxcluster_partial[c];

			indR  = ind*N;

			indM  = ind*d;

			indS  = ind*d2;

			sumR  = 0.0;

			for (n = 0 + indR; n < N + indR ; n++)
			{

				sumR  += R[n];

			}

			P[ind]    = sumR*cteN;                /* Weight's update */


			invsumR   = 1.0/(sumR + fudge);

			for (i = 0 ; i < d ; i++)
			{

				temp     = 0.0;

				nd       = i;

				for (n = 0 ; n < N ; n++)
				{

					temp  += (Z[nd] * R[n + indR] );

					nd    += d;

				}

				temp       *= invsumR;

				vect[i]     = temp;     /* mu(:,k) = (X*R(:,k)) ./ (sumR+fudge);, mean's update */

				M[i + indM] = temp; 


			}

			nd      = 0;

			for (n = 0 ; n < N ; n++)
			{

				for (i = 0 ; i < d ; i++)
				{

					dZ[i + nd]  = (Z[i + nd] - vect[i]);  /* (d x N) */

				}


				/* RdZ = dZ.*(repmat(R(: , k)' , d , 1) = (d x N) */


				temp  = R[n + indR];

				for (i = 0 + nd ; i < d  + nd ; i++)
				{

					invSkdZ[i]  = temp*dZ[i] ;  /* (d x N) */

				}

				nd    += d;

			}

			/* Sk = RdZ*Z' (d x N)*(N x d), Sk symetric */


			RZt(invSkdZ , dZ , tempSk , d , N , covtype);

			for(i = 0 ; i < d ; i++)
			{

				tempSk[i*d1] +=lambda;

			}


			invsumR   = 1.0/(sumR + lambda);

			for (i = 0 ; i < d2 ; i++)
			{

				tempSk[i]  *= invsumR;

				S[i + indS] = tempSk[i];

			}


		}

	}

	(*T_final)        = t;

	(*likely_final)   = likely + cte_MDL;


}


/* ---------------------------------------------------------------------------------------------------------------------------- */

void   RZt(double *R , double *Z , 
		   double *S , 
		   int d , int N , int covtype)


		   /* S = R*Z', where S is Symetric. Perform only (d*(d+1)/2) terms for full cov */

{

	register int i , j , n , nd , id , jd , d1=d+1 , d2 = d*d;

	register double sum;

	if(covtype == 0)  /* Full Cov */
	{

		for (j = 0 ; j < d ; j++)
		{

			id          = j;

			for(i = 0 ; i <= j ; i++)
			{

				sum     = 0.0;

				nd      = 0;

				for(n = 0 ; n < N ; n++)
				{

					sum  += R[i + nd]*Z[j + nd];

					nd   += d;
				}

				S[id]      = sum;

				id        += d;
			}
		}

		/* Copy terms  Lower(S)=>Upper(S) */

		jd       = 0;

		for (j = 0 ; j < d ; j++)
		{

			id   = j + jd + d;

			for(i = j+1 ; i < d ; i++)
			{

				S[id]      = S[i + jd];

				id        += d;

			}

			jd   += d;


		}
	}
	if(covtype == 1)   /* Elliptical Cov */
	{


		for(i = 0 ; i < d2 ; i++)
		{

			S[i]     = 0.0;

		}

		id          = 0;

		for(i = 0 ; i < d ; i++)
		{

			sum     = 0.0;

			nd      = i;

			for(n = 0  ; n < N; n++)
			{				
				sum  += R[nd]*Z[nd];	

				nd   += d;
			}

			S[id]   = sum;

			id      += d1;
		}

	}

	if(covtype == 2)   /* Spherical Cov */
	{

		sum     = 0.0;

		nd      = 0;

		for(n = 0 ; n < N ; n++)
		{			

			for (i = 0 + nd ; i < d + nd ; i++)
			{

				sum  += R[i]*Z[i];			

			}

			nd   += d;

		}


		for(i = 0 ; i < d2 ; i++)
		{

			S[i]     = 0.0;

		}


		id           = 0;

		sum          = sum/d;

		for(i = 0 ; i < d ; i++)
		{

			S[id]    = sum;

			id      += d1;

		}

	}	
}

/* ---------------------------------------------------------------------------------------------------------------------------- */


void qsindex (double  *a, int *index , int lo, int hi)
{
	/* 
	lo is the lower index, hi is the upper index
	of the region of array a that is to be sorted

	*/
	int i=lo, j=hi , ind;
	double x=a[(lo+hi)/2] , h;

	/*  partition */
	do
	{    
		while (a[i]<x) i++; 
		while (a[j]>x) j--;
		if (i<=j)
		{
			h        = a[i]; 
			a[i]     = a[j]; 
			a[j]     = h;
			ind      = index[i];
			index[i] = index[j];
			index[j] = ind;
			i++; 
			j--;
		}
	}
	while (i<=j);

	/*  recursion */
	if (lo<j) qsindex(a , index , lo , j);
	if (i<hi) qsindex(a , index , i , hi);
}



/* ---------------------------------------------------------------------------------------------------------------------------- */

double inv(double *S , double *invS , double *temp  , double *vect , double *vv , int *indx , int d , int offset)

{
	int i , j , jd;

	double dd , det = 1.0;

	for(i = 0 ; i < d*d ; i++)
	{

		temp[i] = S[i + offset];

	}


	if(ludcmp(temp , d , indx , &dd , vv ))
	{

		for(i = 0 ; i < d ; i++)

		{

			det *= temp[i + i*d];

		}

		for(j = 0; j < d; j++)
		{            
			for(i = 0; i < d; i++) 

			{
				vect[i] = 0.0;
			}

			jd      = j*d;

			vect[j] = 1.0;

			lubksb(temp , d , indx , vect);

			for(i = 0 ; i < d ; i++) 

			{

				invS[jd + i] = vect[i];

			}
		}

	}

	return (1.0/det);

}



/* ---------------------------------------------------------------------------------------------------------------------------- */

void lubksb(double *m, int n, int *indx, double *b)
{
	int i, ii = -1, ip, j , in;

	double sum;

	for(i = 0 ; i < n; i++)

	{
		ip        = indx[i];

		sum       = b[ip];

		b[ip]     = b[i];

		if(ii > -1)
		{
			for(j = ii; j <= i - 1; j++)
			{

				sum -= m[i + j*n] * b[j];

			}
		} 
		else if(sum)
		{

			ii = i;

		}

		b[i]     = sum;

	}

	for(i = n - 1 ; i >= 0 ; i--)

	{
		sum = b[i];

		in  = i*n;

		for(j = i + 1 ; j < n; j++)
		{
			sum -= m[i + j*n] * b[j];
		}

		b[i] = sum / m[i + in];
	}
}


/* ---------------------------------------------------------------------------------------------------------------------------- */


int ludcmp(double *m, int n, int *indx, double *d , double *vv)
{
	int i, imax, j, k , jn , kn , n1 = n - 1;

	double big, dum, sum , temp;


	d[0] = 1.0;

	for(i = 0; i < n; i++)
	{
		big = 0.0;

		for(j = 0; j < n; j++)
		{
			if((temp = fabs(m[i + j*n])) > big)

			{
				big = temp;
			}

		}
		if(big == 0.0)
		{

			return 0;
		}

		vv[i] = 1.0 / big;
	}

	for(j = 0; j < n; j++)
	{
		jn  = j*n;

		for(i = 0; i < j; i++)

		{
			sum = m[i + jn];

			for(k = 0 ; k < i; k++)

			{
				sum -= m[i + k*n ] * m[k + jn];
			}

			m[i + jn] = sum;
		}

		big = 0.0;

		for(i = j; i < n; i++)

		{
			sum = m[i + jn];

			for(k = 0; k < j; k++)

			{

				sum -= m[i + k*n] * m[k + jn];

			}

			m[i + jn] = sum;

			if((dum = vv[i] * fabs(sum)) >= big)

			{
				big  = dum;

				imax = i;
			}
		}

		if(j != imax)

		{
			for(k = 0; k < n; k++)

			{

				kn            = k*n;

				dum           = m[imax + kn];

				m[imax + kn]  = m[j + kn];

				m[j + kn]     = dum;

			}

			d[0]       = -d[0];

			vv[imax]   = vv[j];
		}

		indx[j] = imax;

		if(m[j + jn] == 0.0)

		{

			m[j + jn] = NUMERICS_FLOAT_MIN;

		}

		if(j != n1)

		{
			dum = 1.0 / (m[j + jn]);

			for(i = j + 1; i < n; i++)

			{

				m[i + jn] *= dum;

			}
		}
	}

	return 1;
}



/* ---------------------------------------------------------------------------------------------------------------------------- */



void randini(void)

{

	/* SHR3 Seed initialization */

	jsrseed  = (UL) time( NULL );

	jsr     ^= jsrseed;


	/* KISS Seed initialization */

#ifdef ranKISS

	z        = (UL) time( NULL );

	w        = (UL) time( NULL );

	jcong    = (UL) time( NULL );

	mix(z , w , jcong);

#endif


}


/* ---------------------------------------------------------------------------------------------------------------------------- */

void randnini(void)
{
	register const double m1 = 2147483648.0;

	register double  invm1;

	register double dn = 3.442619855899 , tn = dn , vn = 9.91256303526217e-3 , q;

	int i;


	/* Ziggurat tables for randn */

	invm1             = 1.0/m1;

	q                 = vn/exp(-0.5*dn*dn);

	kn[0]             = (dn/q)*m1;

	kn[1]             = 0;

	wn[0]             = q*invm1;

	wn[zigstep - 1 ]  = dn*invm1;

	fn[0]             = 1.0;

	fn[zigstep - 1]   = exp(-0.5*dn*dn);

	for(i = (zigstep - 2) ; i >= 1 ; i--)
	{
		dn              = sqrt(-2.*log(vn/dn + exp(-0.5*dn*dn)));

		kn[i+1]         = (dn/tn)*m1;

		tn              = dn;

		fn[i]           = exp(-0.5*dn*dn);

		wn[i]           = dn*invm1;
	}

}


/* ---------------------------------------------------------------------------------------------------------------------------- */


double nfix(void)
{
	const double r = 3.442620; 	/* The starting of the right tail */

	static double x, y;

	for(;;)

	{

		x = hz*wn[iz];

		if(iz == 0)

		{	/* iz==0, handle the base strip */

			do
			{
				x = -log(rand())*0.2904764;  /* .2904764 is 1/r */

				y = -log(rand());
			}

			while( (y + y) < (x*x));

			return (hz > 0) ? (r + x) : (-r - x);
		}

		if( (fn[iz] + rand()*(fn[iz-1] - fn[iz])) < ( exp(-0.5*x*x) ) )

		{

			return x;

		}


		hz = randint;

		iz = (hz & (zigstep - 1));

		if(abs(hz) < kn[iz])

		{
			return (hz*wn[iz]);

		}


	}

}


/* ---------------------------------------------------------------------------------------------------------------------------- */


double randn(void)

{

	hz = randint;

	iz = (hz & (zigstep - 1));

	return (abs(hz) < kn[iz]) ? (hz*wn[iz]) : ( nfix() );

}


/* ---------------------------------------------------------------------------------------------------------------------------- */
