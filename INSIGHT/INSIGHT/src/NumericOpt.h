#ifndef NUMERIC_OPT_H
#define NUMERIC_OPT_H
/** 
   \file NumOptimizers.h 
   Header file for various optimization procedures of one-dimensional functions
*/



/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/** Different possible optimization methods */
typedef enum {
	BFGS_DIRECT,			// uses the BFGS algorithm to maximize F(x)
	BFGS_DERIV_SQUARED,		// uses the BFGS method to minimize (F'(x))^2
	GSL_FINDROOT			// uses Newton's method for finding roots, as implemented in GSL
}OptMethod;




/******************************************************************************************************/
/******                                     INTERFACE FUNCTIONS                                  ******/
/******************************************************************************************************/



/*----------------------------------------------------------*/
/*------------------------   GSL   -------------------------*/
/*----------------------------------------------------------*/

/*---------------------------------
  constants for GSL optimization.
  GSL_MAX_ITER        - maximum number of optimization iterations
  GSL_CONVERGENCE_EPS - double indicating convergence criterion (smaller --> tighter)
------------------------*/
#define GSL_MAX_ITER 20
#define GSL_CONVERGENCE_EPS 1e-4


/** findRootGSL
    Uses GSL procedures to find root of function using calls to function, using a bounding interval
    @param lowerBound lower bound for root search interval
    @param upperBound upper bound for root search interval
    @param funArgs arguments of generic function
    @param fun pointer to function
    @return optimum point found by procedure
*/
double findRootGSL(double lowerBound, double upperBound, void* funArgs, double (*fun)( double x, void* funArgs)  );



/** findRootGSL_FDF
    Uses GSL procedures to find root of function using calls to function and first derivative
    @param initVal initial value for root search
    @param funArgs arguments of generic function
    @param fun pointer to function
    @param deriv pointer to function calculating derivative
    @param funAndDeriv pointer to function calculating function and its derivative
    @return optimum point found by procedure
*/
double findRootGSL_FDF(double initVal, void* funArgs,
						double (*fun)( double x, void* funArgs),
						double (*deriv)( double x, void* funArgs),
						void (*funAndDeriv)( double x, void* funArgs, double* fx, double* dfx)  );
						

/*-----------------------------------------------------------*/
/*------------------------   BFGS   -------------------------*/
/*-----------------------------------------------------------*/

/*---------------------------------
  constants for BFGS optimization.
  BFGS_LOG determines how much will be logged to stdout (-1 --> no output)
  FACTR is the multiple of machine precision that result will be
  PGTOL is size of derivative on exit
  MVAL is not recommended to be higher than 20
  decreasing FACTR, pgtol will increase precision
------------------------*/
#define BFGS_LOG -1
#define FACTR 1.0e10
#define PGTOL 1.0e-4
#define MVAL 20
/** optimizeBFGS
    Uses BFGS procedure to find minimum/maximum of function
    @param initVal a file descriptor for element data
    @return optimum point found by procedure
*/
double optimizeBFGS(double initVal, double lowerBound, double upperBound, void* funArgs, int maxORmin,
					void (*fun)( double x, void* funArgs, double* results, int numDerivs));


#endif
