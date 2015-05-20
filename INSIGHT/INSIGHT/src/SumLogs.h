#ifndef SUM_LOGS_H
#define SUM_LOGS_H
/** 
   \file SumLogs.c 
   Header file of functions for sum-of-logs.
   
   A sum-of-logs function is a function of the following type:
     F(x) = SUM_{i=0..n-1} (A_i*ln(1-x*B_i)) + A_n*ln(x)
   The function has one main argument (x), and auxilliary arguments: 
   - A[] (n+1 array) also called coefficients
   - B[] (n array)   also called factors
   
   This file contains functions for evaluating sum-of-logs and its derivatives, as well as finding
   maximum points
*/


#include "NumericOpt.h"



/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/** Auxilliary arguments for sum-of-logs function */
typedef struct {
	int numTerms;			// number of terms (n)
	double* coefficients;	// As in sumLogs function (array of length numTerms+1)
	double* factors;		// Bs in sumlogs function (array of length numTerms)
	double  minFactor;		// saves minimum factor (-1.0 if uninitialized)
	double  maxFactor;		// saves maximum factor (-1.0 if uninitialized)
	double  minFactorInverse;	// 1/minFactor (-1.0 if uninitialized)
	double  maxFactorInverse;	// 1/maxFactor (-1.0 if uninitialized)
	
}SumLogsArguments;



/******************************************************************************************************/
/******                                     INTERFACE FUNCTIONS                                  ******/
/******************************************************************************************************/


/** evalSumLogs
	evaluates a sum-of-logs function (using log1p for fast computations)
    @param x        main function argument
    @param args     auxilliary arguments 
    @return F_args(x)
*/
double evalSumLogs(double x, SumLogsArguments* args);



/** sumLogDerivatives
	evaluates derivatives of a sum-of-logs function
    @param x          main function argument
    @param args       auxilliary arguments 
    @param numDerivs  number of derivatives to compute
    @param derivs_out pre-allocated array of output values for derivatives
    @return F(x)
*/
void sumLogDerivatives(double x, SumLogsArguments* args, int numDerivs, double* derivs_out);


/** sumLogsBoundaries
	given auxilliary arguments of a sum-of-logs function, estimates bounds (lower and upper) on its maximum point
	uses pre-computed minimum and maximum factors, if available (computes those, if not)
    @param args       auxilliary arguments 
    @param bounds_out a pre-allocated array with two places for lower and upper bound (placed in that order)
*/
void sumLogsBoundaries(SumLogsArguments* args, double* bounds_out);



/** maximizeSumLogs
    computes maximum of a sum-of-logs function using a numerical optimization method.
    @param args       auxilliary arguments 
    @param val        starting point for optimixation. on return holds maximum
    @param funVal     function evaluation at starting point (0 if not computed). on return holds function at maximum
    @param optMethod  indicator for which optimization method to use
    @return 0
*/
int maximizeSumLogs(SumLogsArguments* args, double* val, double* funVal, OptMethod optMethod);


#endif
