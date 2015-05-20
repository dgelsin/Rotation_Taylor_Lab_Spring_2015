/** 
   \file SumLogs.c 
   Implementation of functions for sum-of-logs.
   
   A sum-of-logs function is a function of the following type:
     F(x) = SUM_{i=0..n-1} (A_i*ln(1-x*B_i)) + A_n*ln(x)
   The function has one main argument (x), and auxilliary arguments: 
   - A[] (n+1 array) also called coefficients
   - B[] (n array)   also called factors
   
   This file contains functions for evaluating sum-of-logs and its derivatives, as well as finding
   maximum points
*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include "NumericOpt.h"
#include "SumLogs.h"
#include "Utils.h"

// for debugging numerical maximization procedure
#define DEBUG_MAX_NOT
#define DEBUG_BOUNDARIES_NOT
#define BOUND_PERCISION 1e-12


/******************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATION                             ******/
/******************************************************************************************************/



/** compute F(x) */
double sumLog_wrapper(double x, void* args);

/** compute F'(x)  */
double sumLogDeriv_wrapper(double x, void* args);

/** simultaneously compute F(x) and F'(x) */
void sumLogAndDeriv_wrapper(double x, void* args, double* derivs_out, int dummy);

/** compute F''(x)  */
double sumLogDeriv2_wrapper(double x, void* args);

/** simultaneously compute F'(x) and F''(x) */
void sumLogFirst2Derivs_wrapper(double x, void* args, double* deriv1, double* deriv2);

/** simultaneously compute (F'(x))^2 and its derivative ((F'(x))^2)' = 2F'(X)F''(X) */
void   sumLogDerivSquared_derivs(double x, void* args, double* derivs_out, int dummy);



/******************************************************************************************************/
/******                                     INTERFACE FUNCTIONS                                  ******/
/******************************************************************************************************/


/** evalSumLogs
	evaluates a sum-of-logs function (using log1p for fast computations)
    @param x        main function argument
    @param args     auxilliary arguments 
    @return F_args(x)
*/
double evalSumLogs(double x, SumLogsArguments* args) {
	int term;
	double val = 0.0;

	if(args->coefficients[args->numTerms] > 0.0) {
		if(1.0-x < 1.0) {
			val += args->coefficients[args->numTerms]*log1p(x-1);
		} else {
			val += args->coefficients[args->numTerms]*log(x);
		}
	}
	
	for(term=0; term<args->numTerms; term++) {
		if(args->coefficients[term] > 0.0) {
			val += args->coefficients[term]*log1p(-args->factors[term]*x);
		}
	}// end of for(block)

return val;
}
/** end of evalSumLogs **/



/** sumLogDerivatives
	evaluates derivatives of a sum-of-logs function
    @param x          main function argument
    @param args       auxilliary arguments 
    @param numDerivs  number of derivatives to compute
    @param derivs_out pre-allocated array of output values for derivatives
*/
void sumLogDerivatives(double x, SumLogsArguments* args, int numDerivs, double* derivs_out) {
	int ind, der;
	double res, factor;
	
	/*--- initialize all derivatives with contribution of n'th term ---*/
	res = -args->coefficients[args->numTerms];
	if(x != 0.0) {
		factor = -1.0/x;
	} else if(res == 0.0){
		factor = 1.0;
	} else {
		fprintf(stderr, "Fatal error: sum-log derivatives undefined with non-zero coefficient and log of zero.\n");
		for(der=0; der<numDerivs; der++) {
			derivs_out[der] = 0.0;
		}// end of for(der)
		return;
	}
	
	for(der=0; der<numDerivs; der++) {
		res *= factor;
		if(der > 0) {
			res *= der;
		}
		derivs_out[der] = res;
	}// end of for(der)

	/*--- compute contribution of each of the first 0..n-1 terms ---*/
	for(ind=0; ind<args->numTerms; ind++) {
		res = args->coefficients[ind];
		if(res > 0.0) {
			factor = args->factors[ind]/(1.0-args->factors[ind]*x);
			for(der=0; der<numDerivs; der++) {
				res *= factor;
				if(der > 0) {
					res *= der;
				}
				derivs_out[der] -= res;
			}// end of for(der)
		}
	}// end of for(ind)
	
	return;
}
/** end of sumLogDerivatives **/


/** sumLogsBoundaries
	given auxilliary arguments of a sum-of-logs function, estimates bounds (lower and upper) on its maximum point
	uses pre-computed minimum and maximum factors, if available (computes those, if not)
    @param args       auxilliary arguments 
    @param bounds_out a pre-allocated array with two places for lower and upper bound (placed in that order)
*/
void sumLogsBoundaries(SumLogsArguments* args, double* bounds_out) {
	int term;
	double sumOfProducts;
	double lastCoefficient = args->coefficients[args->numTerms];
	
	/*--- default boundaries in case last coefficient is zero ---*/
	if(lastCoefficient == 0.0) {
		bounds_out[0] = 0.0;
		bounds_out[1] = 0.0;
		return;
	}
	
	/*--- compute minimum and maximum factors, if not pre-computed ---*/
//	if(args->minFactor < 0.0 || args->maxFactor < 0.0) {
/*--- DEBUG		printf(" !!! computing minimum and maximum factors for boundaries !!!\n"); ---*/
		args->minFactor = args->maxFactor = -1.0;
		for(term=0; term<args->numTerms; term++) {
			if(args->coefficients[term] > 0.0) {
				if(args->minFactor == -1.0 || args->factors[term] < args->minFactor) {
					args->minFactor = args->factors[term];
				}
				if(args->maxFactor == -1.0 || args->factors[term] > args->maxFactor) {
					args->maxFactor = args->factors[term];
				}
			}
		}// end of for(term)
		/*--- set inverse of max/min factors  ---*/ 
		args->maxFactorInverse = 1/args->maxFactor;
		args->minFactorInverse = 1/args->minFactor;
//	}
	
	/*--- compute sumOfProducts = SUM_i(A_i*B_i) ---*/
	sumOfProducts = 0.0;
	for(term=0; term<args->numTerms; term++) {
		sumOfProducts += args->coefficients[term]*args->factors[term];
	}// end of for(block)

	/*--- compute lower and upper bounds ---*/
	bounds_out[0] = lastCoefficient / (lastCoefficient*args->maxFactor + sumOfProducts);
	bounds_out[1] = lastCoefficient / (lastCoefficient*args->minFactor + sumOfProducts);
	
	if(bounds_out[1] <=args->maxFactorInverse- BOUND_PERCISION && bounds_out[0] >= BOUND_PERCISION) {
		return;
	}

/*--- DEBUG		printf(" c[%d]=%g, f[%d]=%g,  ", term+1,args->coefficients[term],term+1,args->factors[term]); ---*/
/*--- DEBUG		printf(" c[%d]=%g, | ", term+1,args->coefficients[term]); ---*/
/*--- DEBUG		printf(" sum of products is %g\n", sumOfProducts); ---*/
/*--- DEBUG		printf("\nfactor bounds (%d terms) are [%g,%g], and optimization bounds are [%g,%g]\n", args->numTerms, args->minFactor, args->maxFactor, bounds_out[0], bounds_out[1]); ---*/
	

	/*--- take care of end cases where bounds are near bounds of regions where F(x) is defined ---*/

	if(args->maxFactorInverse < 2*BOUND_PERCISION) {
		bounds_out[0] = bounds_out[1] = args->maxFactorInverse/2.0;
	} else {
		if(bounds_out[0] < BOUND_PERCISION) {
			bounds_out[0] = BOUND_PERCISION;
			if(bounds_out[1] < bounds_out[0]) {
				bounds_out[1] = bounds_out[0];
				return;
			}
		}
		if(bounds_out[1] > args->maxFactorInverse- BOUND_PERCISION) {
/*--- DEBUG					printf("upper bound (%g) exceeds 1/maxFactor (%g).\n", bounds_out[1], args->maxFactorInverse- BOUND_PERCISION); ---*/
			bounds_out[1] = args->maxFactorInverse- BOUND_PERCISION;
			if(bounds_out[1] < bounds_out[0]) {
				bounds_out[0] = bounds_out[1];
				return;
			}
		}
	}
	return;
}
/** end of sumLogsBoundaries **/



/** maximizeSumLogs
    computes maximum of a sum-of-logs function using a numerical optimization method.
    @param args       auxilliary arguments 
    @param val        starting point for optimixation. on return holds maximum
    @param funVal     function evaluation at starting point (0 if not computed). on return holds function at maximum
    @param optMethod  indicator for which optimization method to use
    @return 0
*/
int maximizeSumLogs(SumLogsArguments* args, double* val, double* funVal, OptMethod optMethod) {
	
	double boundaries[2];
	double x, initVal = *val;
	double fx;
	
	/*--- compute boundaries for maximum ---*/
	sumLogsBoundaries(args, boundaries);

#ifdef DEBUG_MAX
	printf("\n -- maximizing in region [%g,%g], with start %g, val %g\n", boundaries[0], boundaries[1],initVal, *funVal);
/*	printf("Computing F(x) and F'(x) in a range of values.\n");
	printf("%10s %10s %10s\n", " x ", " f(x) ", " f'(x) ");
	for(x=boundaries[1]*0.90; x<=boundaries[1]*1.1; x+=boundaries[1]*0.01) {
		printf("%10g %10g %10g\n",x,sumLog_wrapper(x, (void*)args),sumLogDeriv_wrapper(x, (void*)args));
	}
*/
#endif //DEBUG_MAX

	if(boundaries[1] - boundaries[0] < BOUND_PERCISION) {
		/*--- boundaries are identical - don't employ optimization ---*/
		*val = boundaries[0];
		if(funVal != NULL) {
			fx = evalSumLogs(boundaries[0],args);
			if(fx < *funVal) {
				*val = initVal;
				return 1;
			}
			*funVal = fx;
		}
		return 0;
	}
	
	x = initVal;
	
	if(x<boundaries[0] || x>boundaries[1]) {
		x = boundaries[1];
		fx = evalSumLogs(x,args);
#ifdef DEBUG_MAX
		printf(" setting to top boundary %g", x);
		if(funVal != NULL) {
			printf(" val %g", *funVal);
		}
		printf("\n");			
#endif //DEBUG_MAX

if(funVal != NULL && fx > *funVal) {
			*funVal = fx;
			*val = x;
			return 0;
		}
		if(funVal == NULL && fx > evalSumLogs(*val,args)) {
			*val = x;
			return 0;
		}
		x = (boundaries[0] + boundaries[1]) / 2;
#ifdef DEBUG_MAX
		printf(" resetting x to %g, ", x);
#endif //DEBUG_MAX
	}

#ifdef DEBUG_MAX
		printf(" numerically maximizing ");
#endif //DEBUG_MAX
	switch(optMethod) {
		case BFGS_DIRECT:
			/*--- find maximum of F(x) directly using BFGS ---*/
			x = optimizeBFGS(x, boundaries[0], boundaries[1], (void*)args, 1/*maximize*/, sumLogAndDeriv_wrapper);
			break;
		case BFGS_DERIV_SQUARED:
			/*--- find minimum of (F'(x))^2 using BFGS ---*/
			x = optimizeBFGS(x, boundaries[0], boundaries[1], (void*)args, 0/*minimize*/, sumLogDerivSquared_derivs); 
			break;
		case GSL_FINDROOT:
			/*--- use Newton's method to find root of F'(x) [ GSL LIBRARY ] ---*/
//			x = findRootGSL_FDF(x, (void*)args, sumLogDeriv_wrapper, sumLogDeriv2_wrapper, sumLogFirst2Derivs_wrapper); 
			x = findRootGSL(boundaries[0], boundaries[1], (void*)args, sumLogDeriv_wrapper); 
			break;
	}
	

	if(funVal != NULL) {
		fx = evalSumLogs(x,args);
		if(fx < *funVal) {
			*val = initVal;
			return 1;
		}
		*funVal = fx;
	}
	*val    = x;

#ifdef DEBUG_MAX
	printf(" --> found %g, val %g.\n",x, *funVal);
/**/	printf("Computing F(x) and F'(x) in a range of values.\n");
	printf("%10s %10s %10s\n", " x ", " f(x) ", " f'(x) ");
	for(x=*val*0.99; x<=*val*1.01; x+=*val*0.001) {
		printf("%12g %12g %12g\n",x,sumLog_wrapper(x, (void*)args),sumLogDeriv_wrapper(x, (void*)args));
	}
/**/
#endif //DEBUG_MAX
	
	
	return 0; 
}
/** end of maximizeSumLogs **/ 



/******************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATION                             ******/
/******************************************************************************************************/



/** sumLog_wrapper
	a wrapper for a function that computes F(x) for a sum-of-logs function
    @param x          main function argument
    @param args       auxilliary arguments 
    @return F(x)
*/
double sumLog_wrapper(double x, void* args) {
	return evalSumLogs(x, (SumLogsArguments*)args);
}
/** end of sumLog_wrapper **/



/** sumLogDeriv_wrapper
	a wrapper for a function that computes F'(x) for a sum-of-logs function
    @param x          main function argument
    @param args       auxilliary arguments 
    @return F'(x)
*/
double sumLogDeriv_wrapper(double x, void* args) {
	double deriv[1];
	sumLogDerivatives(x, (SumLogsArguments*)args, 1 /*numDerivs*/, deriv);
	return deriv[0];
}
/** end of sumLogDeriv_wrapper **/


/** sumLogAndDeriv_wrapper
	a wrapper for a sum-of-logs function (and derivatives), with args masked out as void*
    @param x          main function argument
    @param args       auxilliary arguments 
    @return F_args(X)
*/
void sumLogAndDeriv_wrapper(double x, void* args, double* derivs_out, int dummy) {
	/*--- expects dummy to be 2 ---*/

	derivs_out[0] =  evalSumLogs(x, (SumLogsArguments*)args);
	sumLogDerivatives(x, (SumLogsArguments*)args, 1 /*numDerivs*/, derivs_out+1);
	
	return;
}
/** end of sumLogAndDeriv_wrapper **/



/** sumLogDeriv2_wrapper
	a wrapper for a function that computes F''(x) for a sum-of-logs function
    @param x          main function argument
    @param args       auxilliary arguments 
    @return F''(x)
*/
double sumLogDeriv2_wrapper(double x, void* args) {
	double derivs[2];
	sumLogDerivatives(x, (SumLogsArguments*)args, 2 /*numDerivs*/, derivs);
	return derivs[1];
}
/** end of sumLogDeriv2_wrapper **/



/** sumLogFirst2Derivs_wrapper
	a wrapper for a function that computes F'(x) and F''(x) for a sum-of-logs function
    @param x          main function argument
    @param args       auxilliary arguments
    @param deriv1     pointer to space where first derivative is to be outputted
    @param deriv2     pointer to space where second derivative is to be outputted
*/
void sumLogFirst2Derivs_wrapper(double x, void* args, double* deriv1, double* deriv2) {
	double derivs[2];
	sumLogDerivatives(x, (SumLogsArguments*)args, 2 /*numDerivs*/, derivs);
	deriv1[0] = derivs[0];
	deriv2[0] = derivs[1];
	return;
}
/** end of sumLogFirst2Derivs_wrapper **/



/** sumLogDerivSquared_derivs
	a wrapper for a function that computes (F'(x))^2 and its first derivative for a sum-of-logs function
    @param x          main function argument
    @param args       auxilliary arguments 
    @return F_args(X)
*/
void sumLogDerivSquared_derivs(double x, void* args, double* derivs_out, int dummy) {
	/*--- expects dummy to be 2 ---*/

	sumLogDerivatives(x, (SumLogsArguments*)args, 2 /*numDerivs*/, derivs_out);
	derivs_out[1] *= 2*derivs_out[0];		// ((F'(x))^2)' = 2F'(x)F"(x)
	derivs_out[0] *=   derivs_out[0];		//  (F'(x))^2
	
	return;
}
/** end of sumLogDerivSquared_derivs **/



/******************************************************************************************************/
/******                                          END OF FILE                                     ******/
/******************************************************************************************************/
