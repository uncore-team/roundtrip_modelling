/*=================================================================
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2006 The MathWorks, Inc.
 *
 *	Evaluation of fitting functions during the fitting of a 3-loglogistic
 *
 *=================================================================*/
/* $Revision: 1.10.6.4 $ */
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define NUMINPUTS 3
/* X (3x1), S (nx1), N */ 

/* Output Arguments */

#define NUMOUTPUTS 2
/* F (3x1), J (3x3) - optional */

#define WARNING_PREFIX "[Loglogistic fitting warning] "


#define SMALLVAL 1e-9

//#define DEBUG


static int calculatefunctions(double a, double b, double c, const double *samples, unsigned n,
							  double *res1, double *res2, double *res3)
/* Evaluate the three functions that are being optimized.
   1-> ok, -1->NaN */
{unsigned f;
 double minx,invc,sum1,sum2,sum3,xma,aux,b2invc,xma2invc,logxma,logb,bc,nn;

	/* The following checks should not execute, since they distort the input parameters */

	if (a<=0.0) 
	{
        #ifdef DEBUG
    		mexPrintf("a");
        #endif
		a=SMALLVAL;
	}
    if (b<=0.0) 
	{
        #ifdef DEBUG
    		mexPrintf("b");
        #endif
		b=SMALLVAL;
	}
    if (c<=0.0) /* Care: a value too close to 0 produces overflows later on */
	{
        #ifdef DEBUG
    		mexPrintf("c");
        #endif
		c=SMALLVAL;
	}

	minx=samples[0];
	for (f=1; f<n; f++)
		if (samples[f]<minx)
			minx=samples[f];
    if (a>=minx) 
	{
        #ifdef DEBUG
    		mexPrintf("m");
        #endif
		a=minx-SMALLVAL;
	}
    
	invc=1.0/c;
	b2invc=pow(b,invc);
	logb=log(b);
	bc=b*c;
	nn=(double)n;

	sum1=0.0;
	sum2=0.0;
	sum3=0.0;
	for (f=0; f<n; f++)
	{
		xma=samples[f]-a;	
		xma2invc=pow(xma,invc);
		logxma=log(xma);
		aux=xma2invc+b2invc;

		sum1+=( (1.0+invc)-(2.0*b2invc/c)/aux )/xma;

		sum2+=1.0/aux;

		sum3+=logxma-2.0*(  log(xma/b)   /*(logxma-logb)*/    /( pow(xma/b,b2invc)   /*xma2invc/b2invc*/+1.0));
	}
	sum2=( nn-2.0*b2invc*sum2 )/(bc);
	sum3=( -nn*(logb+c)+sum3 )/(c*c);

	*res1=sum1;
	*res2=sum2;
	*res3=sum3;
    
    #ifdef DEBUG
            mexPrintf("a=%f, b=%f, c=%f, n=%d, f1=%f, f2=%f, f3=%f\n",a,b,c,n,sum1,sum2,sum3);
    #endif

	return(1);
}

static int calculatejacobian(double a, double b, double c, const double *samples, unsigned n,
							 double *jacob)
/* Evaluate the jacobian matrix of the functions that are being optimized */
{double f1a,f1b,f1c,f2a,f2b,f2c,f3a,f3b,f3c;
 double minx,invc,b2invc,bc,nn,logb,c2,c3;
 double xma,xma2invc,logxma,xma2,xmab,xmab2invc,logxmab,aux,aux2,aux1,aux12;
 unsigned f;
	
	/* The following checks should not execute, since they distort the input parameters */

	if (a<=0.0) 
	{
        #ifdef DEBUG
    		mexPrintf("A");
        #endif
		a=SMALLVAL;
	}
    if (b<=0.0) 
	{
        #ifdef DEBUG
    		mexPrintf("B");
        #endif
		b=SMALLVAL;
	}
    if (c<=0.0) 
	{
        #ifdef DEBUG
    		mexPrintf("C");
        #endif
		c=SMALLVAL;
	}

	minx=samples[0];
	for (f=1; f<n; f++)
		if (samples[f]<minx)
			minx=samples[f];
    if (a>=minx) 
	{
        #ifdef DEBUG
    		mexPrintf("M");
        #endif
		a=minx-SMALLVAL;
	}

	invc=1.0/c;
	b2invc=pow(b,invc);
	logb=log(b);
	bc=b*c;
	nn=(double)n;
	c2=c*c;
	c3=c2*c;

	f1a=0.0;
	f1b=0.0;
	f1c=0.0;
	f2a=0.0;
	f2b=0.0;
	f2c=0.0;
	f3a=0.0;
	f3b=0.0;
	f3c=0.0;
	for (f=0; f<n; f++)
	{
		xma=samples[f]-a;	
		xma2invc=pow(xma,invc);
		logxma=log(xma);
		aux=xma2invc+b2invc;

		xma=samples[f]-a;
		xma2invc=pow(xma,invc);
		logxma=log(xma);
		aux=xma2invc+b2invc;

        xma2=xma*xma;
		xmab=xma/b;
		xmab2invc=pow(xmab,invc);
		logxmab=log(xmab);

		aux2=aux*aux;
		aux1=xmab2invc+1.0;
		aux12=aux1*aux1;

		f1a=f1a+(1.0+invc)/xma2-2.0*b2invc/c*(c*aux+xma2invc)/(c*xma2*aux2);

		f1b=f1b-2.0/(c*xma)*b2invc*xma2invc/(bc*aux2);

		f1c=f1c-1.0/(c2*xma)+2.0*b2invc*(c*aux+logb*xma2invc-xma2invc*logxma)/(c3*xma*aux2);

		f2a=f2a+xma2invc/(xma*aux2);

		f2b=f2b+(c*aux-xma2invc)/aux2;

		f2c=f2c+(c*aux+logb*xma2invc-xma2invc*logxma)/aux2;

		f3a=f3a-1.0/xma+2.0*(c*xmab2invc-xmab2invc*logxmab+c)/(c*xma*aux12);

		f3b=f3b+(xmab2invc*logxmab-c*aux1)/aux12;

		f3c=f3c+logxma+(logxmab*(xmab2invc*logxmab-2.0*c*aux1))/(c*aux12);
	}
	f2a=-2.0*b2invc/(bc*c)*f2a;
	f2b=-nn/(bc*b)+2.0*b2invc/(bc*bc)*f2b;
	f2c=-nn/(bc*c)+2.0*b2invc/(bc*c*c)*f2c;
	f3a=invc*invc*f3a;
	f3b=-nn/(bc*c)-2.0/(bc*c*c)*f3b;
	f3c=2.0*nn*logb/(c*c*c)+nn/(c*c)-2.0/(c*c*c)*f3c;

	jacob[0]=f1a;
	jacob[1]=f2a;
	jacob[2]=f3a;

	jacob[3]=f1b;
	jacob[4]=f2b;
	jacob[5]=f3b;

	jacob[6]=f1c;
	jacob[7]=f2c;
	jacob[8]=f3c;
	
	
	
/*	for (f=0; f<9; f++) jacob[f]=(double)f;*/
	return(1);
}

static void printusage(void)
{
	mexPrintf("USAGE: \n"
			  "  F = fittingfunctionsLL3( abccolumnvector, nonorderedsamplescolumnvector, lengthofsecondvector )\n"
			  "\t\tF -> (a;b;c)\n"	
			  " [F,J] = fittingfunctionsLL3( ... )\n");
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{double *parms[NUMINPUTS];
 double *output,*jacob,res;
 int resu;
 unsigned f,n;
 unsigned char withjacobian;

    if (nrhs != NUMINPUTS) 
	{
		printusage();
		mexErrMsgTxt("Invalid number of input arguments."); 
	}
    if (nlhs > NUMOUTPUTS)
	{
		printusage();
		mexErrMsgTxt("Invalid number of output arguments."); 
	}

	for (f=0; f<NUMINPUTS; f++)
		parms[f]=mxGetPr(prhs[f]);
	plhs[0]=mxCreateDoubleMatrix(3,1,mxREAL);
	output=mxGetPr(plhs[0]);
	if (nlhs==2)
	{
		withjacobian=1;
		plhs[1]=mxCreateDoubleMatrix(3,3,mxREAL);
		jacob=mxGetPr(plhs[1]);
	}
	else withjacobian=0;
    
    /* Do the actual computations in a subroutine */
	resu=calculatefunctions((parms[0])[0],(parms[0])[1],(parms[0])[2],
					 		parms[1],(unsigned)(*(parms[2])),
							&(output[0]),&(output[1]),&(output[2]));
	if (resu==-1)
		mexWarnMsgTxt(WARNING_PREFIX"Nan in function calculations");
	if (withjacobian)
	{
		resu=calculatejacobian((parms[0])[0],(parms[0])[1],(parms[0])[2],
					 			parms[1],(unsigned)(*(parms[2])),
							 	jacob);
		if (resu==-1)
			mexWarnMsgTxt(WARNING_PREFIX"Nan in jacobian calculations");
	
	}

}


