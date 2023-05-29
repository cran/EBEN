//version2: use blas; user controlled memory usage		Dec152012

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <R_ext/Lapack.h>
#include <stdlib.h>
#ifndef FCONE
# define FCONE
#endif


//B: Binomial; F: full; Neg normal + exp + gamma
void LinearSolverBfNeg(double * a, double *logout, int N,int M,double *output);
void fEBInitializationBfNeg(double *Alpha, double * PHI2, int *Used, int *Unused, double *Mu2,
				double *BASIS, double *Targets, float *Scales, int * initial, int n, int *m, int kdim);



void fEBSigmoidBfNeg(double * y, double * PHI_Mu,int N);
double fEBDataErrorBfNeg(double dataError,double *y,double *PHI_Mu,double *Targets,int N);
void MatrixInverseBfNeg(double * a,int N);
void fEBCatPostModeBfNeg(double * Mu2, double *beta,double *SIGMA2, double * H2, double *PHI2,
								double *Targets, double *Alpha,int N, int M);
void fEBCatFullStatBfNeg(double * beta, double * SIGMA2, double * H2, float *S_in, float * Q_in,
				float * S_out, float * Q_out,  double *BASIS,float * Scales, double *PHI2,
				double * Targets, int * Used, double *Alpha, double * Mu2, double * BasisCache,
				int *n, int *m, int *kdim);

void fEBDeltaMLBfNeg(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
int *Used, int * Unused, float * S_out, float * Q_out, double *Alpha,
double *a_gamma, double *b_gamma, int m, int mBar);

void fEBDeltaMLgroup(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
				int *Used, int * Unused, float * S_out, float * Q_out, double *Alpha,
				double *a_gamma, double *b_gamma, int N_used, int N_unused,int k);

void fEBDeltaMLgroupMemorySave(double *deltaLogMarginal, int *selectedAction, double *newAlpha, int *nu,int *anyToDelete,
				int *Used, int * Unused, float * S_out, float * Q_out, double *Alpha,
				double *a_gamma, double *b_gamma, int N_used, int N_unused,int k);

void fEBBinaryMexBfNeg(int *Used, double *Mu2, double *SIGMA2, double *H2, double *Alpha, double *PHI2,
				double *BASIS, double * Targets, float *Scales, double *a_gamma, double *b_gamma,
				int *iteration, int *n, int *kdim, int *m,double * LOGlikelihood,int basisMax, int verbose, int group);

void fEBBinaryFull(double *BASIS, double * Targets, double *a_gamma, double * b_gamma,
				double * logLIKELIHOOD, double * Beta, double *wald,double *intercept, int *n, int *kdim,
				int *VB, int *bMax,int *gp)
{
	int N					= *n;
	int K					= *kdim;
	int M_full				= (K+1)*K/2;
	int verbose = VB[0];
	int group 	= gp[0];
	const int iter_max		= 50;
	const double err_max	= 1e-8;
	// set a limit for number of basis

	if(verbose >1) Rprintf("start EBLasso-NEG with a: %f, \tb: %f\n",a_gamma[0], b_gamma[0]);
//------------------------------------------------------------
	int basisMax			= *bMax;
	int basisMax_io 		= *bMax;
	if(K>=100) basisMax 	= 10*K;
//------------------------------------------------------------

	if(verbose >1) Rprintf("M_full: %d; basisMax: %d\n",M_full,basisMax);

	double vk				= 1e-30;
	double vk0				= 1e-30;
	double temp				= 0;
	int i,j,l,kk;
	float *Scales			= (float * ) Calloc(M_full, float);
	double scal 			= 1.0;
//		xixi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
	int inci =1;
	int incj =1;
	int inc0 				= 0;
	double *readPtr1;
	double zero_blas 		= 0;
	//assume Beta: Kx 4;
	for (i				=0;i<K;i++)
	{
		//Beta[i]			= -1;
		//Beta[K + i]		= -1;
		//Beta[K*2 + i]	= 0;
		//Beta[K*3 + i]	= 0;
		readPtr1 = &BASIS[i*N];
		temp = F77_CALL(ddot)(&N, readPtr1, &inci,readPtr1, &incj);
		if(temp ==0) temp	= 1;
		Scales[i]			=(float)sqrt(temp);
	}
	//PartII kk
	kk						= K;					//index starts at 0;
	for (i					=0;i<(K-1);i++)
	{
		for (j				=(i+1);j<K;j++)
		{
			temp = 0;
			for(l=0;l<N;l++) 	temp		= temp + BASIS[i*N + l]*BASIS[i*N + l]*BASIS[j*N + l]*BASIS[j*N + l];
			if (temp == 0)		temp		= 1;
			Scales[kk]			=(float)sqrt(temp);
			kk++;
		}
	}
	readPtr1 				= &Beta[basisMax_io*2];
	F77_CALL(dcopy)(&basisMax_io,&zero_blas,&inc0,readPtr1,&inci);  //dcopy(n, x, incx, y, incy) ---> y = x
	readPtr1 				= &Beta[basisMax_io*3];
	F77_CALL(dcopy)(&basisMax_io,&zero_blas,&inc0,readPtr1,&inci);  //dcopy(n, x, incx, y, incy) ---> y = x

	//
	int iter				= 0;
	double err				= 1000;
	double *Mu2, *SIGMA2, *H2, *Alpha, *PHI2;
	int * Used,*iteration, *m;
	//Rprintf("max basis: %d\n",basisMax); 		//R_alloc

	Used					= (int* ) Calloc(basisMax, int);
	Mu2						= (double * ) Calloc(basisMax, double);
	SIGMA2					= (double * ) Calloc(basisMax*basisMax, double);
	H2						= (double * ) Calloc(basisMax*basisMax, double);
	Alpha					= (double * ) Calloc(basisMax, double);
	PHI2					= (double * ) Calloc(N*basisMax, double);
	iteration				= (int* ) Calloc(1, int);
	m						= (int* ) Calloc(1, int);
	//Rprintf("outer loop starts");
	if(verbose >1) Rprintf("outer loop starts\n");
	m[0]			= 2;
	while (iter<iter_max && err>err_max)
	{
		iter				= iter + 1;

		vk0					= vk;
		iteration[0]		= iter;
		fEBBinaryMexBfNeg(Used, Mu2, SIGMA2, H2, Alpha,PHI2,	BASIS, Targets,Scales, a_gamma, b_gamma,
						iteration, n, kdim, m,logLIKELIHOOD,basisMax,verbose,group);

		vk = F77_CALL(dasum)(m, Alpha, &inci);

		err					= fabs(vk - vk0)/m[0];
//Rprintf("Iteration number: %d, err: %f\n",iter,err);
		if(verbose >2) Rprintf("Iteration number: %d, err: %f\n",iter,err);
	}

	// wald score
	int M					= m[0];
	double *tempW			= (double * ) Calloc(M,double);

	wald[0]					= 0;
	int index = 0;
	if(verbose >1) Rprintf("EBLASSO-NEG Finished, number of basis: %d\n",M);
	//Rprintf("number of basis: %d\n",(M-1));
	for(i=0;i<M;i++)
    {
		readPtr1 = &H2[i*M];
		tempW[i] = F77_CALL(ddot)(&M, readPtr1, &inci,Mu2, &incj);
	}
	 wald[0] = F77_CALL(ddot)(&M, tempW, &inci,Mu2, &incj);
	kk = 0;
	int locus1 = 0;
	int locus2 = 0;
	//int Meffect = M;

	//in Binary Model: M is actually length(effect) + 1;
	if(M>basisMax)
	{
		M = basisMax;
	}

		for(i=1;i<M;i++) 						//M is different from number of basis!!
		{
		// blup collection
		//Rprintf("test Used: %d\n",Used[i-1]);
			index				= Used[i-1] - 1; // 0 ~ 960-481-1
			//COMPUTE i,j
			if(index<K)
			{
				Beta[kk] 	= index + 1;
				Beta[basisMax_io+kk] 	= index + 1;
			}else //compute i,j
			{
				l = index;
				for(j=0;j<K;j++)
				{
					l = l -(K-j);
					if(l<0) //only work on negatives
					{
						locus1 = j-1;
						locus2 = l + K;
						break; //get out for loop
					}
				}
				Beta[kk] = locus1 +1;
				Beta[basisMax_io+kk] = locus2  +1;

			}
			scal = (double)Scales[index];
			Beta[basisMax_io*2 + kk]	= Mu2[i]/scal;
			Beta[basisMax_io*3 + kk]  = SIGMA2[i*M + i]/(scal*scal);
			kk = kk + 1;
		}

	//

	intercept[0]	= Mu2[0];
	intercept[1]	= SIGMA2[0];
	//Rprintf("fEB computation compelete!\n");



	Free(Scales);
	Free(Used);
	Free(Mu2);
	Free(SIGMA2);
	Free(H2);
	Free(Alpha);
	Free(PHI2);
	Free(iteration);
	Free(m);
	Free(tempW);
}




/************** outputs are passed by COPY in R, cann't dynamic realloc memory **************************/
/************** Not a problem in C */
// function [Used,Mu2,SIGMA2,H2,Alpha,PHI2]=fEBBinaryMexBfNeg(BASIS,Targets,PHI2,Used,Alpha,Scales,a,b,Mu2,iter)
void fEBBinaryMexBfNeg(int *Used, double *Mu2, double *SIGMA2, double *H2, double *Alpha, double *PHI2,
				double *BASIS, double * Targets, float *Scales, double *a_gamma, double *b_gamma,
				int *iteration, int *n, int *kdim, int *m,double * LOGlikelihood,int basisMax, int verbose, int group)
{
    //basis dimension
   int N,K,M_full,N_used,N_unused,M,i,j,L,kk,iter;
   	N					= *n;			// row number
    K					= *kdim;		// column number
    M_full				= K*(K+1)/2;
    kk					= K;


    double *beta	= (double *)Calloc(N,double);
	int *Unused = (int *) Calloc(M_full,int);
    iter				= *iteration;
	//Rprintf("Iteration number: %d\n",iter);
    const int	ACTION_REESTIMATE       = 0;
	const int	ACTION_ADD          	= 1;
	const int 	ACTION_DELETE        	= -1;
    const int   ACTION_TERMINATE        = 10;

    //[Alpha,PHI2,Used,Unused,Mu2]=InitialCategory(BASIS,Targets,Scales,PHI2,Used,Alpha,Mu2,IniLogic)
    int *IniLogic;
	IniLogic			= (int*) Calloc(1,int);
    if (iter<=1)
    {
        IniLogic[0]     = 0;
        m[0]            = 2;
		M				= m[0];
		N_used			= 1;
		N_unused		= M_full -1;

    }else
    {
		IniLogic[0]    = 1;
        M				= *m;          //Used + 1
		N_used			= M -1;
		N_unused		= M_full - N_used;
    }
    //
	//Rprintf("N_used is: %d; N_unused:%d, M: %d,sample size: %d \n",N_used,N_unused,M,N);

	fEBInitializationBfNeg(Alpha, PHI2, Used, Unused, Mu2, BASIS, Targets, Scales, IniLogic, N, m, K);
	if(verbose >3) Rprintf("\t Initialized basis %d, Alpha: %f, \n", Used[0],Alpha[0]);
	//Rprintf("\t Initialized basis %d, Alpha: %f, Mu: %f \n", Used[0],Alpha[0],Mu2[1]);
	//for(i=0;i<10;i++) Rprintf("PHI2: %f \t  %f; BASIS: %f\n",PHI2[i],PHI2[N+i],BASIS[181*N+i]/Scales[181]);
    double *basisCache;
	basisCache         = (double *) Calloc(N*K,double);
	int inci = 1;
	int incj = 1;
	double *readPtr1,*readPtr2;
	double scal;

    for(i=0;i<K;i++)
    {
        for(j=0;j<N;j++)            basisCache[i*N+j] = BASIS[i*N+j]*BASIS[i*N+j];
    }

	float *S_in, *Q_in, *S_out, *Q_out;
	S_in				= (float *) Calloc(M_full,float);
	Q_in				= (float *) Calloc(M_full,float);
	S_out				= (float *) Calloc(M_full,float);
	Q_out				= (float *) Calloc(M_full,float);
    //[beta,SIGMA2,Mu2,S_in,Q_in,S_out,Q_out,Intercept] ...
    //                   	= FullstatCategory(BASIS,Scales,PHI2,Targets,Used,Alpha,Mu2,BASIS_CACHE)
	fEBCatFullStatBfNeg(beta, SIGMA2, H2, S_in, Q_in, S_out, Q_out,  BASIS,Scales, PHI2,
				Targets, Used, Alpha, Mu2, basisCache,n, m, kdim);
//Rprintf("\t SIGMA2: %f, %f \t H2: %f, %f\n\t %f, %f, \t %f, %f\n",SIGMA2[0],SIGMA2[1], H2[0], H2[1],SIGMA2[2],SIGMA2[3],H2[2],H2[3]);
//Rprintf("\t 182th S_out: %f, Q_out: %f\n",S_out[181],Q_out[181]);
//for(i=0;i<10;i++) Rprintf("S_out: %f \t Q_out: %f\n",S_out[i],Q_out[i]);
    //              For: [DeltaML,Action,AlphaRoot,anyToDelete]     = fEBDeltaMLBfNeg(Used,Unused,S_out,Q_out,Alpha,a,b);
    double *DeltaML, *AlphaRoot, deltaLogMarginal,*phi,newAlpha,oldAlpha;
    double deltaInv,kappa,Mujj,s_ii,*tmp,*tmpp,mu_i;
    //
	int *Action, *anyToDelete,selectedAction;
	anyToDelete			= (int*) Calloc(1,int);
	DeltaML				=	(double *) Calloc(M_full,double);
	AlphaRoot			=	(double *) Calloc(M_full,double);
	Action				= (int *) Calloc(M_full,int);
  	phi					= (double *) Calloc(N,double);

  	tmp					= (double *) Calloc(basisMax,double);
  	tmpp				= (double *) Calloc(basisMax,double);

    int nu,jj,index;
    jj					= -1;
    int anyWorthwhileAction,UPDATE_REQUIRED;
    // mexPrintf("cols:%d \n",M);
  	//
    int i_iter;
    i_iter              = 0;
    int LAST_ITERATION  = 0;
	double logLikelihood,dL,logL0;
	logLikelihood		= 1e-30;
	dL					= 1e-30;
	double *PHI_Mu;
	PHI_Mu				= (double*) Calloc(N,double);
	//Rprintf("check point 3: before loop \n");
	if(verbose >3) Rprintf("check point 3: before loop \n");
    while(LAST_ITERATION!=1)
    {
        i_iter						= i_iter + 1;
		if(verbose >4) Rprintf("\t inner loop %d \n",i_iter);
		//Rprintf("\t inner loop %d \n",i_iter);
        //M							= N_used + 1;
     	//N_unused					= M_full - N_used;
		logL0						= logLikelihood;
        //[DeltaML,Action,AlphaRoot,anyToDelete]     = fEBDeltaMLBfNeg(Used,Unused,S_out,Q_out,Alpha,a,b);

		if(group ==0)
		{
		fEBDeltaMLBfNeg(DeltaML, Action, AlphaRoot,anyToDelete,Used, Unused, S_out, Q_out, Alpha,
				a_gamma,b_gamma, N_used, N_unused);
		}else		//group eblasso
		{
		fEBDeltaMLgroup(DeltaML, Action, AlphaRoot,anyToDelete,Used, Unused, S_out, Q_out, Alpha,
				a_gamma,b_gamma, N_used, N_unused,K);
		}
		deltaLogMarginal			= 0;
        nu							= -1;
		//fEBDeltaMLgroupMemorySave(&deltaLogMarginal, &selectedAction, &newAlpha,&nu,anyToDelete,Used, Unused, S_out, Q_out, Alpha,
		//		a_gamma,b_gamma, N_used, N_unused,K);


        for(i=0;i<M_full;i++)
        {
            if(DeltaML[i]>deltaLogMarginal)
            {
                deltaLogMarginal    = DeltaML[i];
                nu                  = i;
            }
        }
        //
        if(nu==-1)
		{
			anyWorthwhileAction     = 0;
			selectedAction          = -10;
		}else
		{
			anyWorthwhileAction	= 1;
		    selectedAction              = Action[nu];
			newAlpha                    = AlphaRoot[nu];
		}
		//Rprintf("\t ActionOn nu= %d, deltaML: %f, selectedAction: %d\n",nu+1, DeltaML[nu],selectedAction);
        if(selectedAction==ACTION_REESTIMATE || selectedAction==ACTION_DELETE)
        {
            index                   = nu + 1;
            for(i=0;i<N_used;i++)
            {
                if (Used[i]==index)	jj  = i;
            }
        }
        kk                          = K;
        for(i=0;i<K;i++)
        {
            if (i==nu)
            {
				readPtr1 = &BASIS[i*N];
				scal = 1.0/(double)Scales[i];
				F77_CALL(dcopy)(&N,readPtr1,&inci,phi,&incj);
				F77_CALL(dscal)(&N,&scal,phi,&inci);

            }else if(i<(K-1))
            {
                for(j=i+1;j<K;j++)
                {
                    if(kk==nu)
                    {
                        for(L=0;L<N;L++)    phi[L] =BASIS[i*N+L]*BASIS[j*N+L]/(double)Scales[kk];
                    }
                    kk              = kk + 1;
                }
            }
        }
        //newAlpha                    = AlphaRoot[nu];
        if(anyWorthwhileAction==0)  selectedAction = ACTION_TERMINATE;
        if(selectedAction==ACTION_REESTIMATE)
        {
            if (fabs(log(newAlpha)-log(Alpha[jj]))<=1e-3 && anyToDelete[0] ==0)
            {
    			//cout<<"Terminated by small deltaAlpha on basis:"<<nu<<endl;
                selectedAction		= ACTION_TERMINATE;
            }
        }
	//Rprintf("\t selectedAction: %d\n",selectedAction);
        //
        UPDATE_REQUIRED				= 0;
        if(selectedAction==ACTION_REESTIMATE)
        {
			if(verbose >4) Rprintf("\t\t Action: Reestimate : %d \t deltaML: %f\n",nu + 1, deltaLogMarginal);
			//Rprintf("\t\t Action: Reestimate : %d \t deltaML: %f\n",nu + 1, deltaLogMarginal);
            oldAlpha				= Alpha[jj];
            Alpha[jj]				= newAlpha;
            index					= jj + 1;

            deltaInv				= 1.0/(newAlpha-oldAlpha);
            kappa					= 1.0/(SIGMA2[index*M+index] + deltaInv);
            Mujj					= Mu2[jj+1];
			scal = - Mujj *kappa;
			readPtr1 = &SIGMA2[index*M];
			F77_CALL(daxpy)(&M,&scal,readPtr1,&inci,Mu2,&incj);

            UPDATE_REQUIRED			= 1;
        }
        /////////////////////////////////////////////////////////////////////////////////
        else if(selectedAction==ACTION_ADD)
        {
			if(verbose >4) Rprintf("\t\t Action:add : %d \t deltaML: %f\n",nu + 1,deltaLogMarginal);
			//Rprintf("\t\t Action:add : %d \t deltaML: %f\n",nu + 1,deltaLogMarginal);
			//Rprintf("\t\t newAlpha: %f\n",newAlpha);
            // B_phi*PHI2*SIGMA2        tmp = B_phi*PHI2
            index					= N_used + 1;
			if(index > (basisMax -10)&& (N*K) > 1e7) {
				Rprintf("bases: %d, warning: out of Memory, alloc more to Neffect!\n",index);
			}
            for(i=0;i<index;i++)
            {
                tmp[i]				= 0;
                for(j=0;j<N;j++) tmp[i] = tmp[i] + beta[j]*phi[j]*PHI2[i*N+j]; //M+1   x 1
            }
            for(i=0;i<index;i++)
            {
				readPtr1 = &SIGMA2[i*index];
				tmpp[i] = F77_CALL(ddot)(&index, tmp, &inci,readPtr1, &incj);

            }//tmpp is tmp in matlab.

            //Alpha: M -> M + 1
			Alpha[N_used]			= newAlpha;                     //new element

            //
			readPtr1 = &PHI2[index*N];
			F77_CALL(dcopy)(&N,phi,&inci,readPtr1,&incj);
            //

            s_ii					= 1.0/(newAlpha + (double)S_in[nu]);
            mu_i					= -s_ii*(double)Q_in[nu];

			F77_CALL(daxpy)(&index,&mu_i,tmpp,&inci,Mu2,&incj);

            Mu2[index]				= mu_i;							//new element
            //

            Used[N_used]			= nu + 1;						//new element
			N_used					= N_used + 1;

            //
			N_unused				= N_unused - 1;
            for(i=0;i<N_unused;i++)
            {
                if(Unused[i]== (nu + 1))		Unused[i] =Unused[N_unused];
            }
			m[0]					= N_used + 1;
			M						= m[0];
            UPDATE_REQUIRED			= 1;
			//for(i=1;i<M;i++) Rprintf(" \t\t basis: %d :new weight: %f \n",Used[i-1],Mu2[i]);

        }
        //////////////////////////////////////////////////////////////////////////////////
        else if(selectedAction==ACTION_DELETE)
        {
			if(verbose >4) Rprintf("\t\t Action: delete : %d deltaML: %f \n",nu + 1,deltaLogMarginal);
			//Rprintf("\t\t Action: delete : %d deltaML: %f \n",nu + 1,deltaLogMarginal);
            index					= N_used - 1;
            //
			Alpha[jj]				= Alpha[index];           //Alpha: M -> M - 1

            //
			readPtr1 = &PHI2[(jj+1)*N];
			readPtr2 = &PHI2[N_used*N];
			F77_CALL(dcopy)(&N,readPtr2,&inci,readPtr1,&incj);


            Mujj					= Mu2[jj+1];
			scal = - Mujj/SIGMA2[(jj+1)*M + jj+1];
			readPtr1 = &SIGMA2[(jj+1)*M];
			F77_CALL(daxpy)(&M,&scal,readPtr1,&inci,Mu2,&incj);
//for(i=0;i<M;i++)  Mu2[i] = Mu2[i] - Mujj*SIGMA2[(jj+1)*M +i]/SIGMA2[(jj+1)*M + jj+1];

            //
			Mu2[jj+1]				= Mu2[N_used];

            //Used; Unused;
            Used[jj]				= Used[index];
			N_used					= index;

            //
			N_unused				= N_unused + 1;
			Unused[N_unused -1]		= nu + 1;

			m[0]					= N_used + 1;
			M						= m[0];
            UPDATE_REQUIRED			= 1;
        //for(i=1;i<M;i++) Rprintf(" \t\t basis: %d :new weight: %f \n",Used[i-1],Mu2[i]);
		}
        //Rprintf("\t\t Update_required: %d\n",UPDATE_REQUIRED);
        //
        if(UPDATE_REQUIRED==1)
        {
			fEBCatFullStatBfNeg(beta, SIGMA2, H2, S_in, Q_in, S_out, Q_out,  BASIS,Scales, PHI2,
			Targets, Used, Alpha, Mu2, basisCache,n, m, kdim);

        }
		//Rprintf("\t\t selected Action: %d\n",selectedAction);
		//
        if(selectedAction==ACTION_TERMINATE) LAST_ITERATION =1;
        if(i_iter==1000)   LAST_ITERATION = 1;
		logLikelihood = 0;
		for(i = 0;i<N;i++)
		{
			readPtr1 = &PHI2[i];
			PHI_Mu[i] = F77_CALL(ddot)(&M, readPtr1, &N,Mu2, &incj);
//for(j = 0;j<M;j++)		PHI_Mu[i]	= PHI_Mu[i] + PHI2[j*N+i]*Mu2[j];
			logLikelihood			= logLikelihood + Targets[i]*log(exp(PHI_Mu[i])/(1+exp(PHI_Mu[i]))) +
										(1-Targets[i])*log(1/(1+exp(PHI_Mu[i])));
		}
		dL							= fabs((logLikelihood - logL0)/logL0);
		//Rprintf("\t\t likelihoodChange: %f\n",dL);
		if(dL <1e-3) LAST_ITERATION = 1;
		//Rprintf("\t\t Last_iteration value: %d\n",LAST_ITERATION);

    }

	LOGlikelihood[0] = logLikelihood;


	Free(beta);
	Free(Unused);
	Free(IniLogic);
	Free(basisCache);
	Free(S_in);
	Free(Q_in);
	Free(S_out);
	Free(Q_out);
	Free(anyToDelete);
	Free(DeltaML);
	Free(AlphaRoot);
	Free(Action);
	Free(phi);
	Free(tmp);
	Free(tmpp);
	Free(PHI_Mu);
}

/****************************************************************************/

// [Alpha,PHI2,Used,Unused,Mu2]=InitialCategory(BASIS,Targets,Scales,PHI2,Used,Alpha,Mu2,IniLogic)  //IniLogic: whether input is empty or not

void fEBInitializationBfNeg(double *Alpha, double * PHI2, int *Used, int *Unused, double *Mu2,
				double *BASIS, double *Targets, float *Scales, int * initial, int N, int *m, int K)
{
    //basis dimension
    int M,M_full,N_used,i,j,k,kk,index;

    M_full					= K*(K+1)/2;
	int IniLogic			= *initial;
    //INPUT
    if(IniLogic==0)            // is empty
    {
		m[0]				= 2;
		M					= m[0];
        N_used				= 1;
    }else                   // not empty
    {
       	N_used              = m[0]-1;
		M					= m[0];
    }
    //output
    const double init_alpha_max     = 1e3;
    const double init_alpha_min     = 1e-3;

	if(IniLogic==0)            // is empty
    {
		//Rprintf("\t Inside Initialization, M: %d, K: %d\n",M, K);
        double *TargetPseudo,proj_ini,proj;
        int loc1 = 0;
		int loc2 = 0;
		int inci =1;
		int incj = 1;
		double *readPtr1;
		double scal;

        TargetPseudo		= (double *) Calloc(N,double);
		scal = -1;
		F77_CALL(dcopy)(&N,&scal,&loc1,TargetPseudo,&inci);
		scal = 2;
		F77_CALL(daxpy)(&N, &scal,Targets, &inci,TargetPseudo, &incj);
//for(i=0;i<N;i++)            TargetPseudo[i]     = 2*Targets[i] -1;
        kk					= K;
        proj_ini			= 0;
		Used[0]			= 1;
        for(i=0;i<K;i++)
        {
			readPtr1 = &BASIS[i*N];
			proj = F77_CALL(ddot)(&N, readPtr1, &inci,TargetPseudo, &incj);
//for(j=0;j<N;j++)                proj    = proj + BASIS[i*N+j]*TargetPseudo[j];
            proj			= proj/(double)Scales[i];
            if(fabs(proj) > fabs(proj_ini))
            {
                proj_ini    = proj;
                loc1        = i;
                loc2        = i;
                Used[0]		= i + 1;
            }
        }
        for(i=0;i< (K - 1);i++)
        {
            for(j= (i+ 1); j<K; j++)
            {
                proj			= 0;
                for(k=0;k<N;k++)            proj    = proj + BASIS[i*N+k]*BASIS[j*N+k]*TargetPseudo[k];
                proj            = proj/(double)Scales[kk];
                if(fabs(proj) > fabs(proj_ini))
                {
                    proj_ini    = proj;
                    loc1        = i;
                    loc2        = j;
                    Used[0]		= kk  + 1;
                }
                kk              = kk + 1;
            }
        }
        //PHI2, duplicate for linear solver
		double *PHIqr;
		PHIqr					= (double *) Calloc(N*M,double);
		scal = 1;
		i = 0;
		F77_CALL(dcopy)(&N,&scal,&i,PHI2,&inci);
		F77_CALL(dcopy)(&N,&scal,&i,PHIqr,&inci);

        //
        double * PHI;
        PHI						= (double *) Calloc(N,double);
        if(loc1==loc2)
        {
			readPtr1 = &BASIS[loc1*N];
			F77_CALL(dcopy)(&N,readPtr1,&inci,PHI,&incj);
			scal = 1/(double)Scales[loc1];
			F77_CALL(dscal)(&N,&scal, PHI,&inci);
			readPtr1 = &PHI2[N];
			F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
			readPtr1 = &PHIqr[N];
			F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
        }else
        {
            index				= Used[0] -1;
           	for(i=0;i<N;i++)
            {
                PHI[i]			= BASIS[loc1*N+i]*BASIS[loc2*N+i]/(double)Scales[index];
			}
			readPtr1 = &PHI2[N];
			F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
			readPtr1 = &PHIqr[N];
			F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
        }
		double *logout;
        logout                      = (double *) Calloc(N,double);
        for (i=0;i<N;i++)            logout[i]               = log(((TargetPseudo[i] * 0.9 + 1)/2)/(1-(TargetPseudo[i] * 0.9 + 1)/2));
		// Call function
		LinearSolverBfNeg(PHIqr, logout, N, M, Mu2);

        if(Mu2[1] == 0) Alpha[0] 	= 1;
        else            Alpha[0]    = 1/(Mu2[1]*Mu2[1]);
        if(Alpha[0]< init_alpha_min) Alpha[0]				= init_alpha_min;
        if(Alpha[0]> init_alpha_max) Alpha[0]				= init_alpha_max;


		Free(TargetPseudo);
		Free(PHIqr);
		Free(PHI);
		Free(logout);

	}

	int IsUsed			= 0;
    kk                  = 0;
    for(i=0;i<M_full;i++)
    {
        IsUsed          = 0;
        for(j=0;j<N_used;j++)
        {
            //index   = Used[j];
            if ((i+1)==Used[j])     IsUsed  = 1;
        }
        if(IsUsed==0)
        {
            Unused[kk]  = (i+ 1);
            kk          = kk + 1;
        }
    }
}
void LinearSolverBfNeg(double * a, double *logout, int N, int M,double *output)
{
	const int nrhs		= 1;
	const double Rcond	= 1e-5;
	int rank			= M;
	int *jpvt;
	jpvt				= (int * ) Calloc(M,int);
	const int lwork	= M*N + 4*N;
	double * work;
	work				= (double *) Calloc(lwork,double);

	int info			= 0;
	// *************************Call LAPACK library ************************
	F77_CALL(dgelsy)(&N, &M, &nrhs, a, &N, logout, &N,jpvt,&Rcond, &rank, work, &lwork,&info);
	if(info!=0)
	{
		Rprintf("Call linear solver degls error!\n");
		return;
	}

	//Rprintf("Matrix inversed!\n");
	//output				= (double *) realloc(NULL,M*sizeof(double));
	int inci = 1;
	int incj = 1;
	F77_CALL(dcopy)(&M,logout,&inci,output,&incj);
//for (i=0;i<M;i++) output[i] = logout[i];


	Free(jpvt);
	Free(work);
}

 /// ***********************************************************************************************


//[Mu2 beta SIGMA2] 	= fEBCatPostModeBfNeg(PHI2,Targets,Alpha,Mu2);


/// *********[beta,SIGMA2,Mu2,S_in,Q_in,S_out,Q_out,BASIS_B_PHI,Intercept] ...
///                       	= FullstatCategory(BASIS,Scales,PHI2,Targets,Used,Alpha,Mu2,BASIS_CACHE) ************
    //Mu2 is the same size of M in PHI2; one dimension more than Alpha
    //Function fEBCatPostModeBfNeg:
    //    [Mu2 U beta]     	= fEBCatPostModeBfNeg(PHI2,Targets,Alpha,Mu2);
    //Targets: nx1

void fEBCatFullStatBfNeg(double * beta, double * SIGMA2, double * H2, float *S_in, float * Q_in,
				float * S_out, float * Q_out,  double *BASIS,float * Scales, double *PHI2,
				double * Targets, int * Used, double *Alpha, double * Mu2, double * BasisCache,
				int *n, int *m, int* kdim)
{
    //basis dimension
    int N,K,M,i,j,L,p,kk;
   	N					= *n;			// row number
    K					= *kdim;		// column number
    //M_full				= K*(K+1)/2;
    M					= *m;
    kk					= K;

    //output mxArrays; BASIS_B_PHI is not counte

    //[Mu2 Ui beta SIGMA2]     	= fEBCatPostModeBfNeg(PHI2,Targets,Alpha,Mu2);
	fEBCatPostModeBfNeg(Mu2, beta,SIGMA2, H2, PHI2,Targets, Alpha,N, M);

    //	y				= fEBSigmoidBfNeg(PHI2 * Mu2);
    double *PHI_Mu,*y;
    PHI_Mu        		= (double *) Calloc(N,double);
	y					= (double *) Calloc(N,double);
	int inci =1;
	int incj =1;
	double *readPtr1, *readPtr2;
	double scal;
    for(i=0;i<N;i++)
    {
		readPtr1 = &PHI2[i];
		PHI_Mu[i] = F77_CALL(ddot)(&M, readPtr1, &N,Mu2, &incj);
//for(j=0;j<M;j++)            PHI_Mu[i]           = PHI_Mu[i] + PHI2[j*N+i]*Mu2[j];
    }
	fEBSigmoidBfNeg(y, PHI_Mu,N);

    //e=Targets-y
    double *e;
    e					= (double *) Calloc(N,double);
	F77_CALL(dcopy)(&N,Targets,&inci,e,&incj);
	scal = -1;
	F77_CALL(daxpy)(&N, &scal,y, &inci,e, &incj);
//for(i=0;i<N;i++)    e[i]        = Targets[i] - y[i];

    //Main loop
        //temp parameters: BPvector
    double *BPvector,*temp,*temp1,*temp2;
    BPvector			= (double *) Calloc(M,double);
    temp				= (double *) Calloc(M,double);
    temp1				= (double *) Calloc(M*N,double);
    temp2				= (double *) Calloc(N,double);
    double tempSum,BBsquare,tempZE;
    //Cache BASIS.^2	outside the Inner loop of the program: save computation
	//double Print;
	//Rprintf("K is: %d\n",K);
    for(i=0;i<K;i++)
    {
		scal = (double)Scales[i];
        for(p=0;p<M;p++)
        {
            BPvector[p]     = 0;
            for(j=0;j<N;j++)
            {
                temp1[p*N+j]= BASIS[i*N+j]*PHI2[p*N+j]*beta[j];
                BPvector[p] = BPvector[p] + temp1[p*N+j];
            }
            BPvector[p]     = BPvector[p]/scal;
        }
        //temp             	= (BPvector*Ui).^2;
		//temp				= BPvector*SIGMA*BPvector'

        tempSum             = 0;
		for(p=0;p<M;p++)
        {
			readPtr1  = &SIGMA2[p*M];
			temp[p]      	= F77_CALL(ddot)(&M, readPtr1, &inci,BPvector, &incj);
		}
//for(j=0;j<M;j++)                temp[p]     = temp[p] + BPvector[j]*SIGMA2[p*M+j];

		//temp[p]         = temp[p]*BPvector[p];
        //tempSum         = tempSum + temp[p];
        tempSum = F77_CALL(ddot)(&M, temp, &inci,BPvector, &incj);
		for(p=0;p<N;p++)
        {
            temp2[p]           = BASIS[i*N+p]*e[p];
        }

        //S_in(i)           = beta'*BASIS2(:,i)/Scales(i)^2 -sum(temp);
		readPtr1 = &BasisCache[i*N];
		BBsquare = F77_CALL(ddot)(&N, beta, &inci,readPtr1, &incj);

		readPtr1 = &BASIS[i*N];
		tempZE = F77_CALL(ddot)(&N, readPtr1, &inci,e, &incj);

        S_in[i]             = (float)(BBsquare/(scal*scal)-tempSum);
        Q_in[i]             = (float)(tempZE/scal);
		//Print				= S_in[i];
		//Rprintf("S_in: %f\n",Print);
        S_out[i]            = S_in[i];
        Q_out[i]            = Q_in[i];
        //Interactions
        if(i!=(K-1))
        {
            for(L=(i+1);L<K;L++)
            {
                for(p=0;p<M;p++)
                {
					readPtr1 = &temp1[p*N];
					readPtr2 = &BASIS[L*N];
                    BPvector[p]     = F77_CALL(ddot)(&N, readPtr1, &inci,readPtr2, &incj);
//for(j=0;j<N;j++)	BPvector[p] = BPvector[p] + temp1[p*N+j]*BASIS[L*N+j];
                }
				scal = 1/(double)Scales[kk];
				F77_CALL(dscal)(&M,&scal, BPvector,&inci);
                //temp             	= (BPvector*Ui).^2;
                for(p=0;p<M;p++)
                {
					readPtr1 = &SIGMA2[p*M];
					temp[p] = F77_CALL(ddot)(&M, readPtr1, &inci,BPvector, &incj);
//for(j=0;j<M;j++)	temp[p]     = temp[p] + BPvector[j]*SIGMA2[p*M+j];
//                    temp[p]         = temp[p]*BPvector[p];
//                    tempSum         = tempSum + temp[p];
                }
				tempSum = F77_CALL(ddot)(&M, temp, &inci,BPvector, &incj);
                //S_in(kk)       	= beta'*(BASIS2(:,i).*BASIS2(:,L))/Scales(kk)^2 -sum(temp);
                BBsquare            = 0;
                for(p=0;p<N;p++)
                {
                    BBsquare        = BBsquare + beta[p]*BasisCache[i*N+p]*BasisCache[L*N+p];
				}
				readPtr1 = &BASIS[L*N];
				tempZE = F77_CALL(ddot)(&N, temp2, &inci,readPtr1, &incj);
//tempZE        	= tempZE + temp2[p]*BASIS[L*N+p];
                scal = (double)Scales[kk];
                S_in[kk]            = (float)(BBsquare/(scal*scal)-tempSum);
                Q_in[kk]            = (float)(tempZE/scal);

                S_out[kk]         	= S_in[kk];
                Q_out[kk]          	= Q_in[kk];
                kk                  = kk +1;
            }
        }
    }// Main loop ends

    //S_out(Used),Q_out(Used)
    //S_out(Used)				= (Alpha .* S_in(Used)) ./ (Alpha - S_in(Used));
    //_out(Used)       			= (Alpha .* Q_in(Used)) ./ (Alpha - S_in(Used));
    int N_used,index;
    N_used						= M-1;
	float tempAlpha;
    for(i=0;i<N_used;i++)
    {
        index					= Used[i] -1;
        // mexPrintf("index: \t %d\n",index);
		tempAlpha = (float)Alpha[i];
        S_out[index]			= tempAlpha*S_in[index]/(tempAlpha-S_in[index]);
        Q_out[index]			= tempAlpha*Q_in[index]/(tempAlpha-S_in[index]);
    	/*Print				= S_out[index];
		Rprintf("S_out: %f\n",Print);
			Print				= Alpha[i];
			Rprintf("Alpha: %f\n",Print);*/
	}

	Free(PHI_Mu);
	Free(y);
	Free(e);
	Free(BPvector);
	Free(temp);
	Free(temp1);
	Free(temp2);

}


// #*******************************************************************
//[Mu2 beta SIGMA2] 			= fEBCatPostModeBfNeg(PHI2,Targets,Alpha,Mu2);
void fEBCatPostModeBfNeg(double * Mu2, double *beta,double *SIGMA2, double * H2, double *PHI2,
								double *Targets, double *Alpha,int N, int M)
{
	int i,j,k,L;
    //control parameters
    const double GRADIENT_MIN  	= 1e-6;
    double temp                 = pow(2.0,8);
    const double STEP_MIN       = 1/temp;
    const int   itsMax    		= 25;

    //Function fEBDataErrorBfNeg [error,y]   =fEBDataErrorBfNeg(PHI_Mu,Targets)
    //Calculate PHI_Mu: size: Nx1
    double *PHI_Mu, *y;
	PHI_Mu						= (double *) Calloc(N,double);
	int inci =1;
	int incj =1;
	int inc0 = 0;
	double *readPtr1;
	double scal;


    for(i = 0;i<N;i++)
    {
		readPtr1 = &PHI2[i];
		PHI_Mu[i] = F77_CALL(ddot)(&M, readPtr1, &N,Mu2, &incj);
//for(j = 0;j<M;j++)		PHI_Mu[i]	= PHI_Mu[i] + PHI2[j*N+i]*Mu2[j];
    }
	double dataError			= 0;
	y							= (double *) Calloc(N,double);
	dataError					= fEBDataErrorBfNeg(dataError,y,PHI_Mu,Targets,N);
    //
    double regulariser			= 0;
    for(i = 1;i<M;i++)			regulariser		= regulariser + Alpha[i-1]*Mu2[i]*Mu2[i]/2;
    double newTotalError,g0,h0;
    newTotalError				= regulariser + dataError;
    double * errorLog,*e,*g2;
    errorLog					= (double *) Calloc(itsMax,double);
	e							= (double *) Calloc(N,double);
	g2							= (double *) Calloc(M,double);
    //
    g0							= 0;
    h0							= 0;
    //setup H
	int countGrad;
    double *DeltaMu, *Mu_new2;
    DeltaMu						= (double *) Calloc(M,double);
    Mu_new2						= (double *) Calloc(M,double);
    double step					= 1.0;

    //*****************************************************************************
    // Loops start here
    for(i = 0;i<itsMax;i++)
    {
        errorLog[i]				= newTotalError;    //log error value
        //Gradient

		F77_CALL(dcopy)(&N,Targets,&inci,e,&incj);
		scal = -1;
		F77_CALL(daxpy)(&N, &scal,y, &inci,e, &incj);
		scal = 1;
		g0 = F77_CALL(ddot)(&N, &scal, &inc0,e, &incj);
		for(j = 0;j<N;j++)
        {
            beta[j]				= y[j]*(1-y[j]);
        }
		h0 = F77_CALL(ddot)(&N, &scal, &inc0,beta, &incj);

        g2[0]					= g0;
        H2[0]					= h0;
        // first row and first column of H
        for(j=1;j<M;j++)
        {
			readPtr1 = &PHI2[j*N];
            g2[j]				= F77_CALL(ddot)(&N, e, &inci,readPtr1, &incj);
            H2[j]				= F77_CALL(ddot)(&N, beta, &inci,readPtr1, &incj);

            g2[j]				= g2[j]-Alpha[j-1]*Mu2[j];
            H2[j*M]				= H2[j];
        }
        //Hessian
        //h01
        //H;
        for(j=1;j<M;j++)
        {
            for(k = 1;k<M;k++)
            {
                H2[k*M+j]       = 0;
                for(L = 0;L<N;L++)		H2[k*M+j]   = H2[k*M+j] + PHI2[j*N+L]*beta[L]*PHI2[k*N+L];
                if(j==k)                H2[k*M+j]   = H2[k*M+j] + Alpha[k-1];
            }
        }
		//save a copy of H
		j = M*M;
		F77_CALL(dcopy)(&j,H2,&inci,SIGMA2,&incj);

		MatrixInverseBfNeg(SIGMA2,M);				//inverse of H2 is needed for wald score
		//Rprintf("INSIDE fEBPostMode: iter %d \n", i);
		//Rprintf("\t SIGMA2: %f, %f \t H2: %f, %f\n\t %f, %f, \t %f, %f\n",SIGMA2[0],SIGMA2[1], H2[0], H2[1],SIGMA2[2],SIGMA2[3],H2[2],H2[3]);
        countGrad               = 0;
        for(j = 0;j<M;j++)
		{
            if(fabs(g2[j])<GRADIENT_MIN)	countGrad++;
		}
        if(countGrad==M)						//end loop
        {
            for (k = 0;k<M;k++)
            {
                readPtr1 =&SIGMA2[k];
				DeltaMu[k]      = F77_CALL(ddot)(&M, g2, &inci,readPtr1, &M);
//for(L = 0;L<M;L++)	DeltaMu[k]  = DeltaMu[k] + SIGMA2[L*M+k]*g2[L];
            }

            break;
        }
        //Comput Full Newton step H^(-1)*g;
         //sigma*g
		for (k = 0;k<M;k++)
		{
				readPtr1 =&SIGMA2[k];
				DeltaMu[k]      = F77_CALL(ddot)(&M, g2, &inci,readPtr1, &M);
//for(L=0;L<M;L++)	DeltaMu[k]  = DeltaMu[k] + SIGMA2[L*M+k]*g2[L];
		}
        //
        step					= 1;
        while (step>STEP_MIN)
        {
			F77_CALL(dcopy)(&M,Mu2,&inci,Mu_new2,&incj);
			F77_CALL(daxpy)(&M,&step,DeltaMu,&inci,Mu_new2,&inci);
//for(j=0;j<M;j++)	Mu_new2[j]      = Mu2[j] + step*DeltaMu[j];
            //
            for(j=0;j<N;j++)
            {
				readPtr1 = &PHI2[j];
				PHI_Mu[j] = F77_CALL(ddot)(&M, readPtr1, &N,Mu_new2, &incj);
//for(k=0;k<M;k++)	PHI_Mu[j]   = PHI_Mu[j] + PHI2[k*N+j]*Mu_new2[k];
            }
            //Function 2:   fEBDataErrorBfNeg
            dataError			= fEBDataErrorBfNeg(dataError,y,PHI_Mu,Targets,N);
            //regulariser update
            regulariser			= 0;
            for(j=1;j<M;j++)	regulariser  	= regulariser + Alpha[j-1]*Mu_new2[j]*Mu_new2[j]/2;
            newTotalError       = dataError + regulariser;
            if(newTotalError>=errorLog[i])      step		= step/2;
            else
            {
				F77_CALL(dcopy)(&M,Mu_new2,&inci,Mu2,&incj);
//for(j=0;j<M;j++)				Mu2[j]		= Mu_new2[j];
                step            =0;
            }
        }//end of while
        if(step==1) break;
    }

	Free(PHI_Mu);
	Free(y);
	Free(errorLog);
	Free(e);
	Free(g2);
	Free(DeltaMu);
	Free(Mu_new2);

}



 //************************************************************************************************
double fEBDataErrorBfNeg(double dataError,double *y,double *PHI_Mu,double *Targets,int N)
{
	int i;
    dataError					= 0;		// scalar
    // Call function fEBSigmoidBfNeg;
	fEBSigmoidBfNeg(y,PHI_Mu,N);
    for (i = 0;i<N;i++)
    {
        if(y[i]!=0)		dataError   = dataError - Targets[i]*log(y[i]);
        if (y[i]!=1)	dataError   = dataError - (1-Targets[i])*log(1-y[i]);
    }
	return dataError;
}

void fEBSigmoidBfNeg(double * y, double * PHI_Mu,int N) // from left to right
{
	int i;
    for (i=0;i<N;i++)	y [i]	= 1/(1+exp(-PHI_Mu[i]));

}

//****************************************************************
//** Matrix Inverse by call lapack library of cholesky decomposition and linear solver *
void MatrixInverseBfNeg(double * a,int N)
{
	const char uplo			= 'U';
	int info				= 0;
	//*************************Call LAPACK library ************************
	F77_CALL(dpotrf)(&uplo, &N,a, &N,&info FCONE);
	if(info!=0)
	{
		Rprintf("Call 1st function. dpotrf error, Ill-conditioned Hessian!\n");
		return;
	}
	F77_CALL(dpotri)(&uplo,&N,a,&N,&info FCONE);
		if(info!=0)
	{
		Rprintf("Call 2nd function dpotri error!\n");
		return;
	}
	//a is upper triangular
	int i,j;
	for (i=1;i<N;i++)
	{
		for(j=0;j<i;j++)	a[j*N+i]	= a[i*N+j];
	}
}

//***************************************************************************************************
void fEBDeltaMLBfNeg(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
				int *Used, int * Unused, float * S_out, float * Q_out, double *Alpha,
				double *a_gamma, double *b_gamma, int N_used, int N_unused)
{
    //basis dimension
    int M_full,i,index;
	anyToDelete[0] 						= 0;
    M_full								= N_used + N_unused;
    //
    //Parameters setup      //Action is return as int
    const int	ACTION_REESTIMATE		= 0;
	const int 	ACTION_ADD				= 1;
	const int	ACTION_DELETE			= -1;
    //const double    ACTION_OUT		= 0;      //for Unused index and stay unused. not needed: mxArray initialized as zeros
    const double    logL_inf			= 0.0;
    double alpha, beta, gamma,delta,root1,root2,logL,oldAlpha,a,b;
	a									= *a_gamma;
	b									= *b_gamma;
    //int   anyToDelete					= 0;
    int   anyToAdd						= 0;
    const int CNPriorityAddition		= 0;
    const int CNPriorityDeletion		= 1;
	//Rprintf("Inside fEBDeltaMLBfNeg: a: %f, b %f \n",a, b);
	double Schange,Qchange;
    for(i=0;i<N_used;i++)
    {
        //anyToDelete				= false;
        index						= Used[i] -1;
		DeltaML[index]				= 0;
        //    mexPrintf("N_used: \t %d,\t%d\n",N_used,index);
		Schange = (double)S_out[index];
		Qchange = (double)Q_out[index];
        alpha						= 2*a +2 + b*Schange - b*Qchange*Qchange;
        beta						= (4*a + 5)*Schange + b*Schange*Schange -Qchange*Qchange;
        gamma						= (2*a + 3)*Schange*Schange;
        delta						= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha>0 && delta>0 && beta<0)
        {
            root1					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            if (logL > logL_inf)
            {
                AlphaRoot[index]    = root1;
                Action[index]       = ACTION_REESTIMATE;
                //
                oldAlpha            = Alpha[i];
                DeltaML[index]  	= 0.5*(log(root1*(oldAlpha+Schange)/(oldAlpha*(root1 + Schange)))+
                                    Qchange*Qchange*(1/(root1+Schange)-1/(oldAlpha+Schange)))+
                                    (a+1)*log((root1*(1+b*oldAlpha))/(oldAlpha*(1+b*root1)));
            }
        }
        //case 2
        else if(alpha==0 && beta<0)
        {
            root1               = -gamma/beta;
            AlphaRoot[index]    = root1;
            Action[index]       = ACTION_REESTIMATE;
            //
            oldAlpha            = Alpha[i];
            DeltaML[index]  	= 0.5*(log(root1*(oldAlpha+Schange)/(oldAlpha*(root1 + Schange)))+
                                    Qchange*Qchange*(1/(root1+Schange)-1/(oldAlpha+Schange)))+
                                    (a+1)*log((root1*(1+b*oldAlpha))/(oldAlpha*(1+b*root1)));
        }
        //case 3
        else if(alpha <0)
        {
            root2               = (- beta-sqrt(delta))/(2*alpha);
            AlphaRoot[index]    = root2;
            Action[index]       = ACTION_REESTIMATE;
            //
            oldAlpha            = Alpha[i];
            DeltaML[index]  	= 0.5*(log(root2*(oldAlpha+Schange)/(oldAlpha*(root2 + Schange)))+
                                    Qchange*Qchange*(1/(root2 + Schange)-1/(oldAlpha + Schange)))+
                                    (a+1)*log((root2*(1 + b*oldAlpha))/(oldAlpha*(1+b*root2)));
        }
        //DELETE
        else if (N_used>1)
        {
            anyToDelete[0]      = 1;
            Action[index]       = ACTION_DELETE;
            oldAlpha            = Alpha[i];
            logL                = (log(oldAlpha/(oldAlpha+Schange))+pow(Qchange,2)/(oldAlpha+Schange))*0.5 + (a+1)*log(b*oldAlpha/(1+b*oldAlpha));
            DeltaML[index]      = - logL;
        }
    }
    //ADDITION
    for(i=0;i<N_unused;i++)
    {
        index					= Unused[i] -1;
		DeltaML[index]			= 0;
		Schange = (double)S_out[index];
		Qchange = (double)Q_out[index];

        alpha					= 2*a + 2 + b*Schange - b*Qchange*Qchange;
        beta					= (4*a + 5)*Schange + b*Schange*Schange -Qchange*Qchange;
        gamma					= (2*a + 3)*Schange*Schange;
        delta					= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha>0 && delta>0 && beta<0)
        {
            root1					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            if (logL > logL_inf)
            {
                AlphaRoot[index]    = root1;
                Action[index]       = ACTION_ADD;
                anyToAdd            = 1;
                //
                DeltaML[index]  	= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            }
        }
        //case 2
        else if(alpha==0 && beta<0)
        {
            root1               = -gamma/beta;
            AlphaRoot[index]    = root1;
            Action[index]       = ACTION_ADD;
            anyToAdd            = 1;
            //
            DeltaML[index]  		= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
        }
        //case 3
        else if(alpha <0)
        {
            root2               = (- beta-sqrt(delta))/(2*alpha);
            AlphaRoot[index]    = root2;
            Action[index]       = ACTION_ADD;
            anyToAdd            = 1;
            //
            DeltaML[index]  	= (log(root2/(root2+Schange))+pow(Qchange,2)/(root2+Schange))*0.5 + (a+1)*log(b*root2/(1+b*root2));
        }
    }

	//
    if((anyToAdd==1 && CNPriorityAddition==1) || (anyToDelete[0]==1 && CNPriorityDeletion==1))
    {
        for(i=0;i<M_full;i++)
        {
            if (Action[i] == ACTION_REESTIMATE)											DeltaML[i]     = 0;
			else if (Action[i] == ACTION_DELETE)
            {
                    if(anyToAdd==1 && CNPriorityAddition==1 && CNPriorityDeletion!=1)    DeltaML[i]     = 0;
            }else if (Action[i] == ACTION_ADD)
            {
                    if(anyToDelete[0] ==1 && CNPriorityDeletion==1 && CNPriorityAddition!=1) DeltaML[i] = 0;
            }
        }
    }
}


//----------------------
//group
//***************************************************************************************************
void fEBDeltaMLgroup(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
				int *Used, int * Unused, float * S_out, float * Q_out, double *Alpha,
				double *a_gamma, double *b_gamma, int N_used, int N_unused,int k)
{
    //basis dimension
    int M_full,i,index,Ninteraction;
	anyToDelete[0] 						= 0;
    M_full								= N_used + N_unused;
	Ninteraction = M_full - k;
	double weight = 1.0;

	//
    //Parameters setup      //Action is return as int
    const int	ACTION_REESTIMATE		= 0;
	const int 	ACTION_ADD				= 1;
	const int	ACTION_DELETE			= -1;
    //const double    ACTION_OUT		= 0;      //for Unused index and stay unused. not needed: mxArray initialized as zeros
    const double    logL_inf			= 0.0;
    double alpha, beta, gamma,delta,root1,root2,logL,oldAlpha,a,b;
	a									= *a_gamma;
	b									= *b_gamma;
	if(Ninteraction != 0)
	{
		weight 		= sqrt(Ninteraction);
		//Rprintf("Inside fEBDeltaML: weight %f \n",weight);
		//b 	= *b_gamma/weight;
	}
    //int   anyToDelete					= 0;
    int   anyToAdd						= 0;
    const int CNPriorityAddition		= 0;
    const int CNPriorityDeletion		= 1;
	//Rprintf("Inside fEBDeltaML: a: %f, b %f \n",a, b);
	double Schange,Qchange;
    for(i=0;i<N_used;i++)
    {
		b									= *b_gamma;
        //anyToDelete				= false;
        index						= Used[i] -1;
		if(index>=k) b = *b_gamma/weight;
		DeltaML[index]				= 0;
		Schange = (double)S_out[index];
		Qchange = (double)Q_out[index];
        //    mexPrintf("N_used: \t %d,\t%d\n",N_used,index);
        alpha						= 2*a +2 + b*Schange - b*Qchange*Qchange;
        beta						= (4*a + 5)*Schange + b*Schange*Schange -Qchange*Qchange;
        gamma						= (2*a + 3)*Schange*Schange;
        delta						= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha>0 && delta>0 && beta<0)
        {
            root1					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            if (logL > logL_inf)
            {
                AlphaRoot[index]    = root1;
                Action[index]       = ACTION_REESTIMATE;
                //
                oldAlpha            = Alpha[i];
                DeltaML[index]  	= 0.5*(log(root1*(oldAlpha+Schange)/(oldAlpha*(root1 + Schange)))+
                                    Qchange*Qchange*(1/(root1+Schange)-1/(oldAlpha+Schange)))+
                                    (a+1)*log((root1*(1+b*oldAlpha))/(oldAlpha*(1+b*root1)));
            }
        }
        //case 2
        else if(alpha==0 && beta<0)
        {
            root1               = -gamma/beta;
            AlphaRoot[index]    = root1;
            Action[index]       = ACTION_REESTIMATE;
            //
            oldAlpha            = Alpha[i];
            DeltaML[index]  	= 0.5*(log(root1*(oldAlpha+Schange)/(oldAlpha*(root1 + Schange)))+
                                    Qchange*Qchange*(1/(root1+Schange)-1/(oldAlpha+Schange)))+
                                    (a+1)*log((root1*(1+b*oldAlpha))/(oldAlpha*(1+b*root1)));
        }
        //case 3
        else if(alpha <0)
        {
            root2               = (- beta-sqrt(delta))/(2*alpha);
            AlphaRoot[index]    = root2;
            Action[index]       = ACTION_REESTIMATE;
            //
            oldAlpha            = Alpha[i];
            DeltaML[index]  	= 0.5*(log(root2*(oldAlpha+Schange)/(oldAlpha*(root2 + Schange)))+
                                    Qchange*Qchange*(1/(root2 + Schange)-1/(oldAlpha + Schange)))+
                                    (a+1)*log((root2*(1 + b*oldAlpha))/(oldAlpha*(1+b*root2)));
        }
        //DELETE
        else if (N_used>1)
        {
            anyToDelete[0]      = 1;
            Action[index]       = ACTION_DELETE;
            oldAlpha            = Alpha[i];
            logL                = (log(oldAlpha/(oldAlpha+Schange))+pow(Qchange,2)/(oldAlpha+Schange))*0.5 + (a+1)*log(b*oldAlpha/(1+b*oldAlpha));
            DeltaML[index]      = - logL;
        }
    }
    //ADDITION
    for(i=0;i<N_unused;i++)
    {
        index					= Unused[i] -1;
		b									= *b_gamma;
		if(index>=k) b = *b_gamma/weight;
		DeltaML[index]			= 0;
		Schange = (double)S_out[index];
		Qchange = (double)Q_out[index];
        alpha					= 2*a + 2 + b*Schange - b*Qchange*Qchange;
        beta					= (4*a + 5)*Schange + b*Schange*Schange -Qchange*Qchange;
        gamma					= (2*a + 3)*Schange*Schange;
        delta					= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha>0 && delta>0 && beta<0)
        {
            root1					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            if (logL > logL_inf)
            {
                AlphaRoot[index]    = root1;
                Action[index]       = ACTION_ADD;
                anyToAdd            = 1;
                //
                DeltaML[index]  	= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            }
        }
        //case 2
        else if(alpha==0 && beta<0)
        {
            root1               = -gamma/beta;
            AlphaRoot[index]    = root1;
            Action[index]       = ACTION_ADD;
            anyToAdd            = 1;
            //
            DeltaML[index]  		= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
        }
        //case 3
        else if(alpha <0)
        {
            root2               = (- beta-sqrt(delta))/(2*alpha);
            AlphaRoot[index]    = root2;
            Action[index]       = ACTION_ADD;
            anyToAdd            = 1;
            //
            DeltaML[index]  	= (log(root2/(root2+Schange))+pow(Qchange,2)/(root2+Schange))*0.5 + (a+1)*log(b*root2/(1+b*root2));
        }
    }

	//
    if((anyToAdd==1 && CNPriorityAddition==1) || (anyToDelete[0]==1 && CNPriorityDeletion==1))
    {
        for(i=0;i<M_full;i++)
        {
            if (Action[i] == ACTION_REESTIMATE)											DeltaML[i]     = 0;
			else if (Action[i] == ACTION_DELETE)
            {
                    if(anyToAdd==1 && CNPriorityAddition==1 && CNPriorityDeletion!=1)    DeltaML[i]     = 0;
            }else if (Action[i] == ACTION_ADD)
            {
                    if(anyToDelete[0] ==1 && CNPriorityDeletion==1 && CNPriorityAddition!=1) DeltaML[i] = 0;
            }
        }
    }
}


//***************************************************************************************************
void fEBDeltaMLgroupMemorySave(double *deltaLogMarginal, int *selectedAction, double *newAlpha, int *nu,int *anyToDelete,
				int *Used, int * Unused, float * S_out, float * Q_out, double *Alpha,
				double *a_gamma, double *b_gamma, int N_used, int N_unused,int k)
{
    //basis dimension
    int M_full,i,index,Ninteraction;
	anyToDelete[0] 						= 0;
    M_full								= N_used + N_unused;
	Ninteraction = M_full - k;
	double weight = 1.0;

	float *DeltaML,*AlphaRoot;
	int *Action;
	DeltaML				=	(float *) Calloc(M_full,float);
	AlphaRoot			=	(float *) Calloc(M_full,float);
	Action				= (int *) Calloc(M_full,int);


	//
    //Parameters setup      //Action is return as int
    const int	ACTION_REESTIMATE		= 0;
	const int 	ACTION_ADD				= 1;
	const int	ACTION_DELETE			= -1;
    //const double    ACTION_OUT		= 0;      //for Unused index and stay unused. not needed: mxArray initialized as zeros
    const double    logL_inf			= 0.0;
    double alpha, beta, gamma,delta,root1,root2,logL,oldAlpha,a,b;
	a									= *a_gamma;
	b									= *b_gamma;
	if(Ninteraction != 0)
	{
		weight 		= sqrt(Ninteraction);
		//Rprintf("Inside fEBDeltaML: weight %f \n",weight);
		//b 	= *b_gamma/weight;
	}
    //int   anyToDelete					= 0;
    int   anyToAdd						= 0;
    const int CNPriorityAddition		= 0;
    const int CNPriorityDeletion		= 1;
	//Rprintf("Inside fEBDeltaML: a: %f, b %f \n",a, b);
	double Schange,Qchange;
    for(i=0;i<N_used;i++)
    {
		b									= *b_gamma;
        //anyToDelete				= false;
        index						= Used[i] -1;
		if(index>=k) b = *b_gamma/weight;
		DeltaML[index]				= 0;
		Schange = (double)S_out[index];
		Qchange = (double)Q_out[index];
        //    mexPrintf("N_used: \t %d,\t%d\n",N_used,index);
        alpha						= 2*a +2 + b*Schange - b*Qchange*Qchange;
        beta						= (4*a + 5)*Schange + b*Schange*Schange -Qchange*Qchange;
        gamma						= (2*a + 3)*Schange*Schange;
        delta						= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha>0 && delta>0 && beta<0)
        {
            root1					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            if (logL > logL_inf)
            {
                AlphaRoot[index]    = (float)( root1);
                Action[index]       = ACTION_REESTIMATE;
                //
                oldAlpha            = Alpha[i];
                DeltaML[index]  	=(float)( 0.5*(log(root1*(oldAlpha+Schange)/(oldAlpha*(root1 + Schange)))+
                                    Qchange*Qchange*(1/(root1+Schange)-1/(oldAlpha+Schange)))+
                                    (a+1)*log((root1*(1+b*oldAlpha))/(oldAlpha*(1+b*root1))));
            }
        }
        //case 2
        else if(alpha==0 && beta<0)
        {
            root1               = -gamma/beta;
            AlphaRoot[index]    = (float)( root1);
            Action[index]       = ACTION_REESTIMATE;
            //
            oldAlpha            = Alpha[i];
            DeltaML[index]  	=(float)( 0.5*(log(root1*(oldAlpha+Schange)/(oldAlpha*(root1 + Schange)))+
                                    Qchange*Qchange*(1/(root1+Schange)-1/(oldAlpha+Schange)))+
                                    (a+1)*log((root1*(1+b*oldAlpha))/(oldAlpha*(1+b*root1))));
        }
        //case 3
        else if(alpha <0)
        {
            root2               = (- beta-sqrt(delta))/(2*alpha);
            AlphaRoot[index]    = (float)(root2);
            Action[index]       = ACTION_REESTIMATE;
            //
            oldAlpha            = Alpha[i];
            DeltaML[index]  	= (float)( 0.5*(log(root2*(oldAlpha+Schange)/(oldAlpha*(root2 + Schange)))+
                                    Qchange*Qchange*(1/(root2 + Schange)-1/(oldAlpha + Schange)))+
                                    (a+1)*log((root2*(1 + b*oldAlpha))/(oldAlpha*(1+b*root2))));
        }
        //DELETE
        else if (N_used>1)
        {
            anyToDelete[0]      = 1;
            Action[index]       = ACTION_DELETE;
            oldAlpha            = Alpha[i];
            logL                = (log(oldAlpha/(oldAlpha+Schange))+pow(Qchange,2)/(oldAlpha+Schange))*0.5 + (a+1)*log(b*oldAlpha/(1+b*oldAlpha));
            DeltaML[index]      = (float)(- logL);
        }
    }
    //ADDITION
    for(i=0;i<N_unused;i++)
    {
        index					= Unused[i] -1;
		b									= *b_gamma;
		if(index>=k) b = *b_gamma/weight;
		DeltaML[index]			= 0;
		Schange = (double)S_out[index];
		Qchange = (double)Q_out[index];
        alpha					= 2*a + 2 + b*Schange - b*Qchange*Qchange;
        beta					= (4*a + 5)*Schange + b*Schange*Schange -Qchange*Qchange;
        gamma					= (2*a + 3)*Schange*Schange;
        delta					= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha>0 && delta>0 && beta<0)
        {
            root1					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1));
            if (logL > logL_inf)
            {
                AlphaRoot[index]    =(float)(root1);
                Action[index]       = ACTION_ADD;
                anyToAdd            = 1;
                //
                DeltaML[index]  	=(float)( (log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1)));
            }
        }
        //case 2
        else if(alpha==0 && beta<0)
        {
            root1               = -gamma/beta;
            AlphaRoot[index]    = (float)( root1);
            Action[index]       = ACTION_ADD;
            anyToAdd            = 1;
            //
            DeltaML[index]  		= (float)((log(root1/(root1+Schange))+pow(Qchange,2)/(root1+Schange))*0.5 + (a+1)*log(b*root1/(1+b*root1)));
        }
        //case 3
        else if(alpha <0)
        {
            root2               = (- beta-sqrt(delta))/(2*alpha);
            AlphaRoot[index]    = (float)( root2);
            Action[index]       = ACTION_ADD;
            anyToAdd            = 1;
            //
            DeltaML[index]  	= (float)((log(root2/(root2+Schange))+pow(Qchange,2)/(root2+Schange))*0.5 + (a+1)*log(b*root2/(1+b*root2)));
        }
    }

	//
    if((anyToAdd==1 && CNPriorityAddition==1) || (anyToDelete[0]==1 && CNPriorityDeletion==1))
    {
        for(i=0;i<M_full;i++)
        {
            if (Action[i] == ACTION_REESTIMATE)											DeltaML[i]     = 0;
			else if (Action[i] == ACTION_DELETE)
            {
                    if(anyToAdd==1 && CNPriorityAddition==1 && CNPriorityDeletion!=1)    DeltaML[i]     = 0;
            }else if (Action[i] == ACTION_ADD)
            {
                    if(anyToDelete[0] ==1 && CNPriorityDeletion==1 && CNPriorityAddition!=1) DeltaML[i] = 0;
            }
        }
    }

        deltaLogMarginal[0]			= 0;
        nu[0]							= -1;
        for(i=0;i<M_full;i++)
        {
            if((double)DeltaML[i]>deltaLogMarginal[0])
            {
                deltaLogMarginal[0]    = (double)DeltaML[i];
                nu[0]                  = i;
            }
        }
		index = nu[0];
		selectedAction[0] = Action[index];
		newAlpha[0]          = (double)AlphaRoot[index];

		//Free
    Free(DeltaML);
	Free(AlphaRoot);
	Free(Action);
}


