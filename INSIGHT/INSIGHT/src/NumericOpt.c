/** 
   \file NumOptimizers.c 
   Various implementation of one-dimensional optimization procedures
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericOpt.h"
#include "bfgs.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define GSL_FDF_OPT_METHOD gsl_root_fdfsolver_secant
#define GSL_F_OPT_METHOD gsl_root_fsolver_brent

// for debugging numerical optimization procedure
#define DEBUG_OPT_NOT

/******************************************************************************************************/
/******                                     INTERFACE FUNCTIONS                                  ******/
/******************************************************************************************************/



/** findRootGSL
    Uses GSL procedures to find root of function using calls to function, using a bounding interval
    @param lowerBound lower bound for root search interval
    @param upperBound upper bound for root search interval
    @param funArgs arguments of generic function
    @param fun pointer to function
    @return optimum point found by procedure
*/
double findRootGSL(double lowerBound, double upperBound, void* funArgs, double (*fun)( double x, void* funArgs)  ) {
	
	double x = upperBound;

#ifdef DEBUG_OPT
	double xPrev = x;
#endif // DEBUG_OPT
	int status;
	int iter;

	const gsl_root_fsolver_type *gslSolverType;
	gsl_root_fsolver *gslSolver;
	gsl_function gslFunc;
     
	gslFunc.function = fun;
	gslFunc.params = funArgs;
     
	gslSolverType = GSL_F_OPT_METHOD;
	gslSolver = gsl_root_fsolver_alloc (gslSolverType);
	iter = 0;

#ifdef DEBUG_OPT
	printf ("\n==> using %s method to find root \n", gsl_root_fsolver_name (gslSolver));
     
	printf ("%5s [%9s, %9s] %9s %10s %9s\n","iter", "lower","upper","root","err(est)","val");
	printf ("  %5d %10.7f %10.7f %10.7f %10.7f %10.7f\n", iter, lowerBound, upperBound, x, x - xPrev, fun(x,funArgs));
#endif // DEBUG_OPT
	gsl_root_fsolver_set (gslSolver, &gslFunc, lowerBound, upperBound);
     
	status = GSL_CONTINUE;
	while (status == GSL_CONTINUE && iter < GSL_MAX_ITER){
		iter++;
		status = gsl_root_fsolver_iterate (gslSolver);
#ifdef DEBUG_OPT
		xPrev = x;
#endif // DEBUG_OPT
		x = gsl_root_fsolver_root (gslSolver);
		lowerBound = gsl_root_fsolver_x_lower(gslSolver);
		upperBound = gsl_root_fsolver_x_upper(gslSolver);
		status = gsl_root_test_interval (lowerBound, upperBound, 0.0, GSL_CONVERGENCE_EPS);
     
#ifdef DEBUG_OPT
		if (status == GSL_SUCCESS) {
			printf ("  Converged:\n");
		}
		printf ("  %5d %10.7f %10.7f %10.7f %10.7f %10.7f\n", iter, lowerBound, upperBound, x, x - xPrev, fun(x,funArgs));
#endif // DEBUG_OPT
	}
     
	gsl_root_fsolver_free (gslSolver);

	return x;
}
/** end of findRootGSL **/



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
						void (*funAndDeriv)( double x, void* funArgs, double* fx, double* dfx)  ) {
	
	double x = initVal;

	double xPrev = x;
	int status;
	int iter;

	const gsl_root_fdfsolver_type *gslSolverType;
	gsl_root_fdfsolver *gslSolver;
	gsl_function_fdf gslFuncAndDeriv;
     
	gslFuncAndDeriv.f = fun;
	gslFuncAndDeriv.df = deriv;
	gslFuncAndDeriv.fdf = funAndDeriv;
	gslFuncAndDeriv.params = funArgs;
     
	gslSolverType = GSL_FDF_OPT_METHOD;
	gslSolver = gsl_root_fdfsolver_alloc (gslSolverType);
	gsl_root_fdfsolver_set (gslSolver, &gslFuncAndDeriv, x);
     
	status = GSL_CONTINUE;
	iter = 0;
	
#ifdef DEBUG_OPT
	printf ("\n==> using %s method to find root \n", gsl_root_fdfsolver_name (gslSolver));
     
	printf ("  %-5s %10s %10s %10s\n", "iter", "root", "err(est)", "val");
		printf ("  %5d %g %g %g\n", iter, x, x - xPrev, fun(x,funArgs));
#endif // DEBUG_OPT
	while (status == GSL_CONTINUE && iter < GSL_MAX_ITER){
		iter++;
		status = gsl_root_fdfsolver_iterate (gslSolver);
		xPrev = x;
		x = gsl_root_fdfsolver_root (gslSolver);
		status = gsl_root_test_delta (x, xPrev, 0.0, GSL_CONVERGENCE_EPS);
     
#ifdef DEBUG_OPT
		if (status == GSL_SUCCESS) {
			printf ("  Converged:\n");
		}
		printf ("  %5d %10.7f %10.7f %10.7f\n", iter, x, x - xPrev, fun(x,funArgs));
#endif // DEBUG_OPT
	}
     
	gsl_root_fdfsolver_free (gslSolver);

	return x;
}
/** end of findRootGSL_FDF **/



/** optimizeBFGS
    Uses BFGS procedure to find minimum/maximum of function
    @param initVal a file descriptor for element data
    @return optimum point found by procedure
*/
double optimizeBFGS(double initVal, double lowerBound, double upperBound, void* funArgs, int maxORmin,
						void (*fun)( double x, void* funArgs, double* results, int numDerivs)) {

	double x, funEvals[2];

	/*--- space allocated for setULbounds() ---*/
	double wa[2*MVAL+4 + 12*MVAL*MVAL + 12*MVAL];
	double dsave[29];
	int    isave[44];
	int    iwa[3];
	char   csave[60];
	int    lsave[4];
	char   task[60];

	/*--- constants used by setULbounds() ---*/
	int    m = MVAL;
	double factr = FACTR;
	double pgtol = PGTOL;
	int    noisy = BFGS_LOG;
	int    numPars = 1;
	int    nbd = 2;
	

	strcpy(task, "START                                                      ");
	task[59] = ' ';
	
	if (fun == NULL) {
		fprintf(stderr,"Error: cannot optimize a null function.\n");
		return lowerBound-1.0;
	}
  

	x = initVal;
	(*fun)(x, funArgs, funEvals,2);
	if(maxORmin) {
		funEvals[0] = -funEvals[0];
		funEvals[1] = -funEvals[1];
	}
	while (1) {
		setULbounds(&numPars, &m, &x, &lowerBound, &upperBound, &nbd, funEvals, funEvals+1, 
				&factr, &pgtol, wa, iwa, task, &noisy, csave, lsave, isave, dsave);
		/*--- DEBUG  printf(" -- interval is [%g,%g] --", lowerBound,upperBound); ---*/

		if (task[0]=='F' && task[1]=='G') {
			(*fun)(x, funArgs, funEvals,2);
			if(maxORmin) {
				funEvals[0] = -funEvals[0];
				funEvals[1] = -funEvals[1];
			}
			if (BFGS_LOG > 1)
				printf("funEval=%f\n", funEvals[0]);
		} else if (strncmp(task, "NEW_X", 5)==0) {
			continue;
		} else {
			break;
		}
	}// end of while(1)
	
	if (BFGS_LOG > 1)
		printf("\n");

	return x;
}
/** end of optimizeBFGS **/



/******************************************************************************************************/
/******                                          END OF FILE                                     ******/
/******************************************************************************************************/

/*******
	
	CLEANUP    PREVIOUS UNUSED VERSION OF optimizeBFGS
	
** optimizeBFGS
    Uses BFGS procedure to find minimum/maximum of function
    @param initVal a file descriptor for element data
    @return optimum point found by procedure
*
double optimizeBFGS_bak(double initVal, double lowerBound, double upperBound, void* funArgs, int maxORmin,
					double (*fun)( double x, void* funArgs),
					double (*funDeriv)( double x, void* funArgs)   ) {

	double x, funEval, derivEval;

	*--- space allocated for setULbounds() ---*
	double wa[2*MVAL+4 + 12*MVAL*MVAL + 12*MVAL];
	double dsave[29];
	int    isave[44];
	int    iwa[3];
	char   csave[60];
	int    lsave[4];
	char   task[60];

	*--- constants used by setULbounds() ---*
	int    m = MVAL;
	double factr = FACTR;
	double pgtol = PGTOL;
	int    noisy = BFGS_LOG;
	int    numPars = 1;
	int    nbd = 2;
	
	int sign = (maxORmin)?(-1):(1);

	strcpy(task, "START                                                      ");
	task[59] = ' ';
	
	if (fun == NULL || funDeriv==NULL) {
		return lowerBound-1.0;
	}
  

	x = initVal;
	funEval   = sign *     (*fun)(x, funArgs);
	derivEval = sign *(*funDeriv)(x, funArgs);

	while (1) {
		setULbounds(&numPars, &m, &x, &lowerBound, &upperBound, &nbd, &funEval, &derivEval, 
				&factr, &pgtol, wa, iwa, task, &noisy, csave, lsave, isave, dsave);

		if (task[0]=='F' && task[1]=='G') {
			funEval   = sign *     (*fun)(x, funArgs);
			derivEval = sign *(*funDeriv)(x, funArgs);
			if (BFGS_LOG > 1)
				printf("funEval=%f\n", funEval);
		} else if (strncmp(task, "NEW_X", 5)==0) {
			continue;
		} else {
			break;
		}
	}// end of while(1)
	
	if (BFGS_LOG > 1)
		printf("\n");

	return x;
}
** end of optimizeBFGS **
*******/
