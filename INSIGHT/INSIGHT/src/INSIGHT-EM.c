/** 
   \file INSIGHT-EM.c 
   Main file for INSIGHT EM algorithm and data processing.
   Contains model inference (through EM) and other related tools.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "Utils.h"
#include "SumLogs.h"

#define FIELD_TOKEN " \t"
#define ID_LEN  100

#define VERSION "1.1"
#define DATE    "February 2013"

/*--- method used for numerical optimization (eta and gamma updates) ---*/
/*--- see OptMethod data type in NumericOpt.h for options ---*/
#define OPT_METHOD GSL_FINDROOT

/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/** Different possible polymorphism classes for sites */
typedef enum {
	MONOMORPHIC,		// for monomorphic sites
	POLY_L,				// for polymorphic sites with low minor allele frequency
	POLY_H				// for polymorphic sites with high minor allele frequency
}PolyClass;



/** Data held for each site */
typedef struct {
	PolyClass polyClass;			// polymorphism class
	int numSamples;					// number of samples in population for which there is data
	double pZeqXmaj;				// prior probability that deep ancestral allele (Z) equals major allele
	double pZeqXmin;				// prior probability that deep ancestral allele (Z) equals minor allele
	char   siteId[ID_LEN];			// string id for site
}SiteData;



/** Auxilliary information held for each site for the selection EM  */
/** This auxilliary information is a function of data at that site and the neutral parameters  */
typedef struct {
	int    numSites;	// number of sites with the exact same pattern
	double neut;		// expressions for P[S='n']
	double noDiv;		// expressions for P[S='s' & Z = X]
	double div;			// expressions for P[S='s' & Z \neq X]
	double wattA;		// Watterson's a
}AuxStatPattern_SelectionEM;


/** Auxilliary information held for each site for the beta1 EM  */
/** This auxilliary information is a function of data at that site and the lambda parameters  */
typedef struct {
	double beta1;		// expressions for P[A=Xmaj]
	double beta3;		// expressions for P[A=Xmin]
}AuxStat_beta1EM;



/** Data held for a collection of sites in blocks */
typedef struct {
	int numBlocks;			// number of blocks in collection
	int numSites;			// total number of sites across all blocks
	int maxNumSamples;		// number of individuals sampled in population
	int* numSitesPerBlock;	// array holding number of sites per block
	char** blockNames;		// array of names for blocks (elements)
	SiteData** siteData;	// 2D array for site data indexed by block id and site id (in block)
}CollectionData;



/** Model Parameters */
typedef struct {
	int numBlocks;			// number of blocks modeled
	/*--- neutral parameters ---*/
	double* theta;			// array of theta neutral polymorphism rate parameters (one per block)
	double* lambda;			// array of lambda neutral divergence parameters (one per block)
	double  beta[3];		// array for beta parameters for the derived allele frequency spectrum (indexed -1)
	/*--- selection parameters ---*/
	double  rho;			// rho parameter   - selection mixture
	double  eta;			// eta parameter   - divergence scale for sites under selection compared to neutral sites
	double  gamma;			// gamma parameter - polymorphism scale for sites under selection compared to neutral sites
	/*--- auxilliary stats ---*/
	double meanLambda;	// weighted mean of lambda across blocks (weighted by num of element sites in each block)
	double meanTheta;		// weighted mean of lambda across blocks (weighted by num of element sites in each block)
	double meanLamThe;	// weighted mean of lambda*theta across blocks (weighted by num of element sites in each block)
}ModelParameters;



/** EM Setup arguments*/
typedef struct {
	/*--- i/o ---*/
	char   inFileName[FILENAME_LEN+1];		// name of main input file name
	char   logFileName[FILENAME_LEN+1];		// name of log file (or empty if no logging)
	char   postFileName[FILENAME_LEN+1];	// name of file to which to output posterior counts (or empty if no posterior count output)
	int    verbose;							// flag for more logging onto screen
	int    logInterval;						// number of iterations between logs to file

	/*--- EM halting ---*/
	int    maxNumIterations;				// maximum number of EM iterations
	double lnLd_epsilon;					// stopping condition for EM (delta in ln-likelihood)

	/*--- EM initialization ---*/
	double rhoInitial;						// initial value for rho   parameter
	double etaInitial;						// initial value for eta   parameter
	double gammaInitial;					// initial value for gamma parameter
	double beta1Initial;					// initial value for beta1 parameter (for beta1EM)

	/*--- EM fixed parameters ---*/
	unsigned short fixRho;					// boolean flag: 1 - fix rho   at initial value
	unsigned short fixEta;					// boolean flag: 1 - fix eta   at initial value
	unsigned short fixGamma;				// boolean flag: 1 - fix gamma at initial value

	unsigned short computeConfidence;		// boolean flag: 1 - compute confidence intervals
	
	unsigned short useSimplifiedModel;		// boolean flag: 1 - compute likelihood based on simplified version of the model
	unsigned short estimateBeta;			// boolean flag: 1 - run EM for estimating the ratio beta1/(beta1+beta3)
}EMsetup;



/** EM Status used for return value of EM */
typedef struct {
	int    numIterations;					// number of iterations for EM run
	double lnLikelihood;					// ln-likelihood of solution
	double lnLd_diff;						// difference in ln-likelihood at last iteration
	char   endState[20];					// string identifying state at end of EM:
											//  "converged" - EM converged
											//  "overshoot" - reached reduction in likelihood (lnLd_doff < 0)
											//  "timeout"   - reached maximum number of iterations
											//  "zero-likelihood" - no possible solution (when restrictions are applied to parameters)
											//  "error"     - general error in EM
}EMstatus;



/******************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATION                             ******/
/******************************************************************************************************/



/*--- main functions ---*/
int     printUsage(const char* progName);
int     processInputArgs(char* argv[], int argc, EMsetup* emSetup);
int     processElementSiteFile(const char* inFileName, EMsetup* emSetup, CollectionData* data, ModelParameters* params);
int     processCountFile(const char* countFileName, CollectionData* data, ModelParameters* params);
int     readData(FILE* dataFile, CollectionData* data);
int     selectionEM(CollectionData* elementData, ModelParameters* params, EMsetup* emSetup, EMstatus* emStatus);
int     beta1EM(CollectionData* data, double* lambda, double* beta1_p, EMsetup* emSetup, EMstatus* emStatus);
int     computePosteriors(CollectionData* elementData, ModelParameters* params, unsigned short useSimplifiedModel, const char* outFileName, double* totalCounts);
int     computeFisherConfidence(CollectionData* elementData, ModelParameters* params, EMsetup* emSetup, double* confidence_out);
int     computeFisherIM(CollectionData* elementData, ModelParameters* params, unsigned short useSimplifiedModel, double** fisherIM);
/*--- auxilliary functions ---*/
int     computeWattersons_a(int numValues, double* valueArray);
int     freeCollectionData(CollectionData* data);
int     freeModelParameters(ModelParameters* params);
void    print3X3matrix(FILE*, double** matrix3X3);
int     invertFisherIM(double** fisherIM, double** varianceMatrix);


/******************************************************************************************************/
/******                                          MAIN                                            ******/
/******************************************************************************************************/



/** printUsage
    Prints usage for program
    @param progName the name of executable as used in commandline
    @return 0
*/
int printUsage(const char* progName) {

	printf("+--------------------------------------------------------------------------------------\n");
	printf("| INSIGHT-EM (v%s)  - program for estimating selection parameters from poly/div patterns\n", VERSION);
	printf("+--------------------------------------------------------------------------------------\n");
	printf("| Usage: ' %s infile [optional flags] '\n",progName);
	printf("|  infile contains a summary of sequence information across a given set of genomic positions\n");
	printf("+--------------------------------------------------------------------------------------\n");
	printf("| optional flags:\n");
	printf("|\n");
	printf("| ~~~~ General ~~~~\n");
	printf("| -h %11s : show this usage message\n","--help");
//	printf("| -s %11s : reverts to simpler version of the model conditioning on known ancestral states\n","--simple");
	printf("| -b %11s : runs EM on neutral 'L' sites to estimate beta1/(beta1+beta3) [ optional initial value, default = 0.5 ]\n","--beta1-3");
	printf("|\n");
	printf("| ~~~~ I/O ~~~~\n");
	printf("| -v %11s : run with more messages outputed to screen\n","--verbose");
	printf("| -f %11s : log filename for EM       ( required if -l --log option is used )\n","--log-file");
	printf("| -p %11s : produce posterior counts of all site types into a specified file\n","--post-cnt");
	printf("| -l %11s : number of iterations between log printouts ( default = 100 )\n","--log-iter");
	printf("| -c %11s : do NOT compute confidence intervals for parameters ( computed by default )\n","--no-conf ");
	printf("|\n");
	printf("| ~~~~ EM halting conditions ~~~~\n");
	printf("| -i %11s : upper bound on number of EM iterations     ( default = 20,000   )\n","--max-iter");
	printf("| -d %11s : ln-likelihood difference at which EM stops ( default = 0.000001 )\n","--min-diff");
	printf("|\n");
	printf("| ~~~~ EM initialization ~~~~\n");
	printf("| -r %11s : initial value for rho   parameter          ( default = 0.6 )\n","--rho-init");
	printf("| -e %11s : initial value for eta   parameter          ( default = 1.0 )\n","--eta-init");
	printf("| -g %11s : initial value for gamma parameter          ( default = 0.5 )\n","--gam-init");
	printf("|\n");
	printf("| ~~~~ EM limit updates ~~~~\n");
	printf("| -fr %10s : do not update rho   parameter\n","--fix-rho");
	printf("| -fe %10s : do not update eta   parameter\n","--fix-eta");
	printf("| -fg %10s : do not update gamma parameter\n","--fix-gam");
	printf("+--------------------------------------------------------------------------------------\n");
	
	return 0;
}
/** end of printUsage **/



/** main
    @return 0 if all OK, and -1 otherwise
*/
int main(int argc, char** argv) {
	char timeString[100];
	int res, block;
	
	EMsetup emSetup;
	EMstatus emStatus;
	ModelParameters params;
	CollectionData data;
	double alpha, tau; 	// ALPHA-TAU
	double Dp, Pw, beta1;
	double* WattersonsArray;
	double totalPosteriorCounts[11];
	
	// confidence intervals for the five selection parameters (rho, eta, gamma, E[Dp], E[Pw], alpha, tau)
	double confidenceIntervals[7];
	
	startTime();
	
	/*--- initialize EM setup arguments and process input arguments ---*/
	res = processInputArgs(argv, argc, &emSetup);
	
	if(res < 0) {
		return -1;
	}
	
	if(emSetup.verbose) {
		printf("----------------------------------------------------------------------------------------------\n");
		printf("       INSIGHT-EM v%s, %s \n",VERSION,DATE);
		printf("----------------------------------------------------------------------------------------------\n");
	}
	
	
	/*--- apply EM for ratio beta1/(beta1+beta3) ---*/

	if(emSetup.estimateBeta) {
		if(emSetup.verbose) {
			printf("==> Estimating the ratio beta1/(beta1+beta3) using neutral sites with low MAF.\n");
			printf("----------------------------------------------------------------------------------------------\n");
			printf("==> Processing site data file '%s' and lambda parameters\n",emSetup.inFileName);
		}
		res = processElementSiteFile(emSetup.inFileName, &emSetup, &data, &params);
		if(res < 0) {
			return -1;
		}

		if(emSetup.verbose) {
			printf("==> Performing EM on ratio beta1/(beta1+beta3) [ initial value = %g ]\n", emSetup.beta1Initial);
			printf("    - EM stops after %d iterations or when log-likelihood increase is below %g\n", emSetup.maxNumIterations, emSetup.lnLd_epsilon);
			if(emSetup.logFileName[0] != '\0') {
				printf("    - progress logged every %d iterations as '.' and in file '%s'\n", emSetup.logInterval, emSetup.logFileName);
			} else {
				printf("    - '.' = %d iterations\n", emSetup.logInterval);
			}
			printf("----------------------------------------------------------------------------------------------\n");
		}

		beta1 = emSetup.beta1Initial;
		res = beta1EM(&data, params.lambda, &beta1, &emSetup, &emStatus);
		
		printTime(timeString);
		if(emSetup.verbose) {
			printf("Done. Running time %s.\n",timeString);
			//printf("----------------------------------------------------------------------------------------------\n");
		}
		printf("\n----------------------------------------------------------\n");
		printf("          %6s \n",             "beta1");
		printf("Estimates: %7lf \n",    beta1);
		printf("----------------------------------------------------------------------------------------------\n");
		printf("            %5s %8s %10s %12s\n",             "iter"  , "lnLd",   "diff"   , "status");
		printf("EM status:  %6d %10g %7g  %s\n",    emStatus.numIterations,   emStatus.lnLikelihood ,  emStatus.lnLd_diff , emStatus.endState);
		printf("----------------------------------------------------------------------------------------------\n");
	
	
		freeModelParameters(&params);
		freeCollectionData(&data);
		return 0;
	}
		
	
	/*--- apply EM for selection parameters ---*/

	/*--- process count file and set neutral parameters ---*/

	if(emSetup.useSimplifiedModel) {
		if(emSetup.verbose) {
			printf("==> Processing count file '%s' and estimating neutral model parameters\n",emSetup.inFileName);
		}
		res = processCountFile(emSetup.inFileName, &data, &params);
	} else {
		if(emSetup.verbose) {
			printf("==> Processing site data file '%s' and extracting neutral model parameters\n",emSetup.inFileName);
		}
		res = processElementSiteFile(emSetup.inFileName, &emSetup, &data, &params);
		if(emSetup.verbose) {
			printf("==> Data file consists of %d sites in %d genomic blocks.\n",data.numSites,data.numBlocks);
		}
	}
	
	
	if(res < 0) {
		return -1;
	}
	
	if(emSetup.verbose) {
		printf("==> Performing EM on selection parameters\n");
		printf("    - EM stops after %d iterations or when log-likelihood increase is below %g\n", emSetup.maxNumIterations, emSetup.lnLd_epsilon);
		if(emSetup.logFileName[0] != '\0') {
			printf("    - progress logged every %d iterations as '.' and in file '%s'\n", emSetup.logInterval, emSetup.logFileName);
		} else {
			printf("    - '.' = %d iterations\n", emSetup.logInterval);
		}
//		printf("    - starting with rho = %f, eta=%f, gamma=%f.\n", emSetup.rhoInitial, emSetup.etaInitial, emSetup.gammaInitial);
		if(emSetup.fixRho || emSetup.fixEta || emSetup.fixGamma) {
			printf("    - fixing the following parameters:");
			if(emSetup.fixRho)
				printf(" rho = %f", emSetup.rhoInitial);
			if(emSetup.fixEta)
				printf(" eta = %f", emSetup.etaInitial);
			if(emSetup.fixGamma)
				printf(" gamma = %f", emSetup.gammaInitial);		
			printf(".\n");
		}			
		if(!emSetup.fixRho || !emSetup.fixEta || !emSetup.fixGamma) {
			printf("    - estimating the following parameters:");
			if(!emSetup.fixRho)
				printf(" rho (init = %f)", emSetup.rhoInitial);
			if(!emSetup.fixEta)
				printf(" eta (init = %f)", emSetup.etaInitial);
			if(!emSetup.fixGamma)
				printf(" gamma (init = %f)", emSetup.gammaInitial);		
			printf(".\n");
		}
		if(emSetup.useSimplifiedModel) {
			printf("    - using simplified version of the model which assumes ancestral states Zi are known.\n");
		} else {
			printf("    - using complete version of the model integrating over assignments to the ancestral states Zi.\n");
		}
		printf("----------------------------------------------------------------------------------------------\n");
	}
	
	selectionEM(&data, &params, &emSetup, &emStatus);

	/*--- compute weighted mean for neutral parameters ---*/
	
		
	/*-----------------------------------------------------------------*/
	/*---------------------   INITIALIZE MATRIX   ---------------------*/
	/*-----------------------------------------------------------------*/
	
	/*--- pre-compute Watterson's a values for all possible numbers of samples ---*/
	params.meanLambda = 0.0;
	params.meanTheta  = 0.0;
	params.meanLamThe = 0.0;
	for(block=0; block<data.numBlocks; block++) {
		params.meanLambda += params.lambda[block]*data.numSitesPerBlock[block];
		params.meanTheta  += params.theta[block] *data.numSitesPerBlock[block];
		params.meanLamThe += params.lambda[block]*params.theta[block]*data.numSitesPerBlock[block];
	}// end of for(block)
	/*--- normalize by number of KBs in elements ---*/
	params.meanLambda /= (((double)data.numSites)/1000);
	params.meanTheta  /= (((double)data.numSites)/1000);
	params.meanLamThe /= (((double)data.numSites)/1000);

	
	WattersonsArray = (double*)malloc((data.maxNumSamples+1) * sizeof(double));
	if(WattersonsArray == NULL) {
		fprintf(stdout, "*** Error: OOM allocating space for Watterson's array in main.\n");
		return -1;
	}
	computeWattersons_a(data.maxNumSamples, WattersonsArray);
	params.meanTheta  *= WattersonsArray[data.maxNumSamples];
	params.meanLamThe *= WattersonsArray[data.maxNumSamples];
	free(WattersonsArray);

	if(params.rho == 0.0 || params.eta == 0) {
		alpha = 0.0;
	} else {
		alpha = 1/(1+(1-params.rho)/(params.rho*params.eta));	// ALPHA-TAU
	}
	
	if(params.rho == 0.0 || params.gamma == 0) {
		tau = 0.0;
	} else {
		tau   = 1/(1+(1-params.rho)/(params.rho*params.gamma));	// ALPHA-TAU
	}

	Dp = params.rho * params.eta * params.meanLambda;
	Pw = params.rho * params.gamma * (params.meanTheta - params.eta * params.meanLamThe);

	printTime(timeString);
	if(emSetup.verbose) {
		printf("Done. Running time %s.\n",timeString);
		// printf("----------------------------------------------------------------------------------------------\n");
	}
	
	if(emSetup.postFileName[0] != '\0') {
		computePosteriors(&data, &params, emSetup.useSimplifiedModel, emSetup.postFileName, totalPosteriorCounts);
		if(emSetup.verbose) {
			printf("----------------------------------------------------------------------------------------------\n");
			printf("==> Posterior counts output to file %s.\n",emSetup.postFileName);
		}
	}		
		
	
	printf("----------------------------------------------------------------------------------------------\n");
	printf("          %6s %8s %8s %8s %8s %8s %8s\n",             "rho", "eta", "gamma", "Dp", "Pw", "alpha", "tau");				// ALPHA-TAU
	printf("Estimates: %7lf %7lf %7lf %7lf %7lf %7lf %7lf\n",    params.rho ,  params.eta ,  params.gamma,   Dp ,  Pw, alpha, tau);	// ALPHA-TAU
	if(emSetup.computeConfidence && 0 != strcmp(emStatus.endState,"zero-likelihood")) {
		computeFisherConfidence(&data, &params, &emSetup, confidenceIntervals);
		printf("StndrdErr: %7lf %7lf %7lf %7lf %7lf %7lf %7lf\n",   confidenceIntervals[0], confidenceIntervals[1], confidenceIntervals[2],			// ALPHA-TAU
												confidenceIntervals[3], confidenceIntervals[4], confidenceIntervals[5], confidenceIntervals[6]);	// ALPHA-TAU
	}
	if(emSetup.postFileName[0] != '\0') {
		printf("Posterior: %7lf -------- -------- %7lf %7lf %7lf %7lf\n",		// ALPHA-TAU
	/*rho*/		   (totalPosteriorCounts[2]+totalPosteriorCounts[5]+totalPosteriorCounts[10])/(totalPosteriorCounts[0]+totalPosteriorCounts[3]+totalPosteriorCounts[6]),
	/*Dp*/	   	(totalPosteriorCounts[5])/(((double)data.numSites)/1000),
	/*Pw*/		   (totalPosteriorCounts[10])/(((double)data.numSites)/1000),
	/*alpha*/	   (totalPosteriorCounts[5])/(totalPosteriorCounts[3]),			// ALPHA-TAU
	/*tau*/		   (totalPosteriorCounts[10])/(totalPosteriorCounts[6]));		// ALPHA-TAU
// DEBUG POSTERIORS		printf("cSel = %f, cNeut = %f, ratio = %f.\n", totalPosteriorCounts[2]+totalPosteriorCounts[5]+totalPosteriorCounts[10],
// DEBUG POSTERIORS			   totalPosteriorCounts[0]+totalPosteriorCounts[3]+totalPosteriorCounts[6]-totalPosteriorCounts[2]-totalPosteriorCounts[5]-totalPosteriorCounts[10],
// DEBUG POSTERIORS			   (totalPosteriorCounts[2]+totalPosteriorCounts[5]+totalPosteriorCounts[10])/(totalPosteriorCounts[0]+totalPosteriorCounts[3]+totalPosteriorCounts[6]));
	}
	printf("----------------------------------------------------------------------------------------------\n");
	printf("            %5s %8s %10s %12s\n",             "iter"  , "lnLd",   "diff"   , "status");
	printf("EM status:  %6d %10g %7g  %s\n",    emStatus.numIterations,   emStatus.lnLikelihood ,  emStatus.lnLd_diff , emStatus.endState);
	printf("----------------------------------------------------------------------------------------------\n");
	
	
	freeModelParameters(&params);
	freeCollectionData(&data);
	

  return 0;
}
/** end of main **/



/******************************************************************************************************/
/******                               INTERNAL FUNCTION IMPLEMENTATION                           ******/
/******************************************************************************************************/



/** processInputArgs
    Processes program input arguments
    @param argv array of string arguments
    @param argc length of argv[]
    @return 0 if all OK, and -1 otherwise
*/
int processInputArgs(char* argv[], int argc, EMsetup* emSetup) {
	int ind;
	unsigned short res, isInFileName, doLog;
	double tmpDouble;
	
	if(argc < 2) {
		printUsage(argv[0]);
		printf("| Error: need to supply input file name.\n");
		printf("+--------------------------------------------------------------------------------------\n");
		printf("\n");
		return -1;
	}
	
	/*--- init EM setup parameters ---*/
	emSetup->inFileName[0]    = '\0';
	emSetup->logFileName[0]   = '\0';
	emSetup->postFileName[0]   = '\0';
	emSetup->verbose          = 0;
	emSetup->maxNumIterations = 20000;
	emSetup->lnLd_epsilon     = 0.000001;
	emSetup->logInterval      = 100;
	emSetup->rhoInitial       = 0.6;
	emSetup->etaInitial       = 1.0;
	emSetup->gammaInitial     = 0.5;
	emSetup->beta1Initial     = 0.5;
	emSetup->fixRho           = 0;
	emSetup->fixEta           = 0;
	emSetup->fixGamma         = 0;
	emSetup->computeConfidence = 1;
	emSetup->useSimplifiedModel = 0;
//	emSetup->useSimplifiedModel = 1;
	emSetup->estimateBeta = 0;
	
	/*--- parse input options ---*/
	isInFileName = 0;
	doLog = 0;
	for(ind=1; ind<argc; ind++) {
		/*--- options with NO additional arguments---*/
		if ( 0 == strcmp(argv[ind], "-h") || 0 == strcmp(argv[ind],     "--help") ) {
			printUsage(argv[0]);
			printf("\n");
			return -1;
		} else if( 0 == strcmp(argv[ind], "-v") || 0 == strcmp(argv[ind], "--verbose") ) {
			emSetup->verbose = 1;
			continue;
		} else if( 0 == strcmp(argv[ind], "-s!!")) {
			/*--- depracated option. Should be used only internally with old-style count files ---*/
			emSetup->useSimplifiedModel = 1;
			continue;
		} else if( 0 == strcmp(argv[ind], "-c") || 0 == strcmp(argv[ind], "--no-conf") ) {
			emSetup->computeConfidence = 0;
			continue;
		} else if( 0 == strcmp(argv[ind], "-fr") || 0 == strcmp(argv[ind], "--fix-rho") ) {
			emSetup->fixRho = 1;
			continue;
		} else if( 0 == strcmp(argv[ind], "-fe") || 0 == strcmp(argv[ind], "--fix-eta") ) {
			emSetup->fixEta = 1;
			continue;
		} else if( 0 == strcmp(argv[ind], "-fg") || 0 == strcmp(argv[ind], "--fix-gam") ) {
			emSetup->fixGamma = 1;
			continue;
		} else if( 0 == strcmp(argv[ind], "-b") || 0 == strcmp(argv[ind], "--beta1-3") ) {
			emSetup->estimateBeta = 1;
			if(ind<argc-1) {
				ind++;
				res = sscanf(argv[ind], "%lf", &tmpDouble);
				if(res<1)   ind--;
				else        emSetup->beta1Initial = tmpDouble;
			}
			if(emSetup->beta1Initial <= 0.0 || emSetup->beta1Initial>=1.0) {
				printUsage(argv[0]);
				printf("| Error: initial value for beta1 should be number in (0,1) (%g).\n", emSetup->beta1Initial);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
			continue;
		} else if(ind == argc-1) {
			if(emSetup->inFileName[0] == '\0') {
				isInFileName = 1;
			} else {
				printUsage(argv[0]);
				printf("| Error: last option '%s' invalid or requires an additional argument.\n", argv[ind]);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
		}
		
		/*--- options with additional arguments---*/
		if(isInFileName) {
			/*--- do noting ---*/
		} else if ( 0 == strcmp(argv[ind], "-i") || 0 == strcmp(argv[ind], "--max-iter") ) {
			ind++;
			res = sscanf(argv[ind], "%d", &emSetup->maxNumIterations);
			if(res<1 || emSetup->maxNumIterations < 1) {
				printUsage(argv[0]);
				printf("| Error: --max-iter option expects positive integer number (%s).\n", argv[ind]);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
		} else if( 0 == strcmp(argv[ind], "-d") || 0 == strcmp(argv[ind], "--min-diff") ) {
			ind++;
			res = sscanf(argv[ind], "%lf", &emSetup->lnLd_epsilon);
			if(res<1 || emSetup->lnLd_epsilon < 0) {
				printUsage(argv[0]);
				printf("| Error: --min-diff option expects positive real number (%s).\n", argv[ind]);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
		} else if( 0 == strcmp(argv[ind], "-f") || 0 == strcmp(argv[ind], "--log-file") ) {
			ind++;
			strncpy(emSetup->logFileName, argv[ind], FILENAME_LEN);
			if(strlen(argv[ind]) > FILENAME_LEN) {
				printf("Warning: logfile name %s is truncated to %d characters.\n", argv[ind], FILENAME_LEN);
			}
		} else if( 0 == strcmp(argv[ind], "-p") || 0 == strcmp(argv[ind], "--post-cnt") ) {
			ind++;
			strncpy(emSetup->postFileName, argv[ind], FILENAME_LEN);
			if(strlen(argv[ind]) > FILENAME_LEN) {
				printf("Warning: file name for posterior counts %s is truncated to %d characters.\n", argv[ind], FILENAME_LEN);
			}
		} else if( 0 == strcmp(argv[ind], "-l") || 0 == strcmp(argv[ind], "--log-iter") ) {
			ind++;
			doLog = 1;
			res = sscanf(argv[ind], "%d", &emSetup->logInterval);
			if(res<1 || emSetup->logInterval < 1) {
				printUsage(argv[0]);
				printf("| Error: --log-iter option expects positive integer number (%s).\n", argv[ind]);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
		} else if( 0 == strcmp(argv[ind], "-r") || 0 == strcmp(argv[ind], "--rho-init") ) {
			ind++;
			res = sscanf(argv[ind], "%lf", &emSetup->rhoInitial);
			if(res<1 || emSetup->rhoInitial < 0.0 || emSetup->rhoInitial > 1.0) {
				printUsage(argv[0]);
				printf("| Error: --rho-init option expects real number in [0,1] (%s).\n", argv[ind]);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
		} else if( 0 == strcmp(argv[ind], "-e") || 0 == strcmp(argv[ind], "--eta-init") ) {
			ind++;
			res = sscanf(argv[ind], "%lf", &emSetup->etaInitial);
			if(res<1 || emSetup->etaInitial < 0) {
				printUsage(argv[0]);
				printf("| Error: --eta-init option expects positive real number (%s).\n", argv[ind]);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
		} else if( 0 == strcmp(argv[ind], "-g") || 0 == strcmp(argv[ind], "--gam-init") ) {
			ind++;
			res = sscanf(argv[ind], "%lf", &emSetup->gammaInitial);
			if(res<1 || emSetup->gammaInitial < 0) {
				printUsage(argv[0]);
				printf("| Error: --gam-init option expects positive real number (%s).\n", argv[ind]);
				printf("+--------------------------------------------------------------------------------------\n");
				printf("\n");
				return -1;
			}
		} else if(emSetup->inFileName[0] == '\0' ) {
			isInFileName = 1;
		} else {
			printUsage(argv[0]);
			printf("| Error: unknown option '%s' ( '%s' assumed to be input filename ).\n", argv[ind], emSetup->inFileName);
			printf("+--------------------------------------------------------------------------------------\n");
			printf("\n");
			return -1;
		}
		
		/*--- if argument does not match, it is assumed to be input filename ---*/
		if(isInFileName) {
			isInFileName = 0;
			strncpy(emSetup->inFileName, argv[ind], FILENAME_LEN);
			if(strlen(argv[ind]) > FILENAME_LEN) {
				printf("Warning: infile name '%s' is truncated to first %d characters ('%s').\n", argv[ind], FILENAME_LEN, emSetup->inFileName);
			}
		}
	}// end of for(ind)

	/*--- check validity of parameter settings ---*/
	res = 1;
	if(doLog && emSetup->logFileName[0] == '\0') {
		if(res == 1) {
			printUsage(argv[0]);
		}
		res = 0;
		printf("| Error: --log-iter option was used without log file name. Use --log-file option\n");
	}
	if(!emSetup->fixRho && (emSetup->rhoInitial== 0.0 || emSetup->rhoInitial== 1.0)) {
		if(res == 1) {
			printUsage(argv[0]);
		}
		res = 0;
		printf("| Error: Should not initialize rho at boundary value %f without explicitly fixing it (--fix-rho).\n",emSetup->rhoInitial);
	}
	if(!emSetup->fixEta && emSetup->etaInitial== 0.0) {
		if(res == 1) {
			printUsage(argv[0]);
		}
		res = 0;
		printf("| Error: Should not initialize eta at boundary value %f without explicitly fixing it (--fix-eta).\n",emSetup->etaInitial);
	}
	if(!emSetup->fixGamma && emSetup->gammaInitial== 0.0) {
		if(res == 1) {
			printUsage(argv[0]);
		}
		res = 0;
		printf("| Error: Should not initialize gamma at boundary value %f without explicitly fixing it (--fix-eta).\n",emSetup->gammaInitial);
	}
	
	if(res == 0) {		
		printf("+--------------------------------------------------------------------------------------\n");
		printf("\n");
		return -1;
	}
	
	
	return 0;
}
/** end of processInputArgs **/



/** processElementSiteFile
    Reads data from element site file
    Initializes neutral parameters according to "block" statements in file
    Summarizes element data in CollectionData data structure
    @param inFileName a file name for element data
    @param emSetup EM setup structure
    @param data Empty structure to be initialized
    @param params Empty structure to be initialized
    @return 0 if all OK, and -1 otherwise
*/
int processElementSiteFile(const char* inFileName, EMsetup* emSetup, CollectionData* data, ModelParameters* params) {
	char line[101];
	char* token;
	int block,numBlocks,site,i, res;
	
	/*--- constant used as default minimal value for theta and lambda in case they are estimated to be 0 ---*/
	double minThetaLambda = 0.00001;
	
	int numSites, numSamples;
	
	// unsigned short ignoreBlock;
	
	if(data == NULL || params == NULL) {
		fprintf(stdout, "*** Error: cannot initialize null CollectionData or ModelParameters in processElementSiteFile().\n");
		return -1;
	}
	
	FILE* inFile = fopen (inFileName, "r");
	if(inFile == NULL) {
		fprintf(stdout, "*** Error: cannot open input file '%s' in processElementSiteFile().\n",inFileName);
		return -1;
	}
	
	/*--- first pass on file to get number of blocks, set beta parameters ---*/
	/*--- and check general structure of input file  ---*/
	res = 0;
	numSites = 0;
	numBlocks = 0;
	numSamples = -1;
	params->beta[0] = params->beta[1] = params->beta[2] = -1.0;
	while ( fgets(line, 100, inFile) != NULL) {
		token = strtok(line, FIELD_TOKEN);
		if(0 == strcmp(token, "samples")) {
			token = strtok(NULL, FIELD_TOKEN);
			res =  sscanf(token, "%d",&numSamples);
			if(res < 1 || numSamples < 0) {
				fprintf(stdout, "*** Error: expecting positive number of chromosome samples in 'samples' line, not %s, in processElementSiteFile().\n",
						token);
				res = -1;
				break;
			}
		} else if(0 == strcmp(token, "beta")) {
			/*--- set beta parameters ---*/
			if(params->beta[0] == -1.0) {
				for(i=0; i<3; i++) {
					token = strtok(NULL, FIELD_TOKEN);
					res =  sscanf(token, "%lf",&params->beta[i]);
					if(res < 1 || params->beta[i] < 0.0 || params->beta[i] > 1.0) {
						fprintf(stdout, "*** Error: expecting number in [0,1] for beta_%d, not %s, in processElementSiteFile().\n",
								i+1, token);
						res = -1;
						break;
					}
				}// end of for(i)
				if(res < 0.0) {
					break;
				}
				if(fabs(params->beta[0] + params->beta[1] + params->beta[2] - 1.0) > 0.00001) {
					fprintf(stdout, "*** Error: sum of beta parameters should be 1.0, not %g+%g+%g=%g, in processElementSiteFile().\n",
								params->beta[0], params->beta[1], params->beta[2], params->beta[0] + params->beta[1] + params->beta[2]);
					res = -1;
					break;
				}
			} else {
				fprintf(stdout, "*** Error: expecting only one 'beta' line in element site file, in processElementSiteFile().\n");
				res = -1;
				break;
			}
		} else if(0 == strcmp(token, "block")) {
			numBlocks++;
		} else if(0 == strcmp(token, "site")) {
			numSites++;
			if(numBlocks == 0) {
				fprintf(stdout, "*** Error: found 'site' line before any block line, in processElementSiteFile().\n");
				res = -1;
				break;
			}
		} else if(token != NULL && token[0] != '#') {
			/*--- allow for empty lines and comment lines ---*/
			fprintf(stdout, "*** Error: illegal first token '%s' found in line of element site file, in processElementSiteFile().\n", token);
			res = -1;
			break;
		}
	}// end of while(inFile)
	
	if(res == 0 && !emSetup->estimateBeta && params->beta[0] == -1.0) {
		fprintf(stdout, "*** Error: expecting 'beta' line in element site file to set beta parameters, in processElementSiteFile().\n");
		res = -1;
	}
	
	if(res < 0) {
		fclose(inFile);
		return -1;
	}
	
	/*--- allocate memory for data and parameters ---*/
	
	data->numSitesPerBlock = NULL;
	data->blockNames = NULL;
	data->siteData = NULL;
	params->theta  = NULL;
	params->lambda = NULL;
	
	data->numSitesPerBlock = (int*)malloc(numBlocks *  sizeof(int));
	data->blockNames = (char**)malloc(numBlocks *  sizeof(char*));
	if(data->blockNames != NULL) {
		data->blockNames[0] = NULL;
		data->blockNames[0] = (char*)malloc(numBlocks * ID_LEN * sizeof(char));
	}
	data->siteData = (SiteData**)malloc(numBlocks *  sizeof(SiteData*));
	params->theta  = (double*) malloc(numBlocks *  sizeof(double));
	params->lambda = (double*) malloc(numBlocks *  sizeof(double));
	
	if(data->siteData == NULL || data->numSitesPerBlock == NULL || data->blockNames == NULL || data->blockNames[0] == NULL || params->theta == NULL || params->lambda == NULL) {
		fprintf(stdout, "*** Error: OOM allocating arrays for CollectionData and ModelParameters in processElementSiteFile().\n");
		freeModelParameters(params);
		freeCollectionData(data);
		fclose(inFile);		
		return -1;
	}
	
	/*--- second pass on file to get lambda and beta parameters ---*/
	/*--- and number of sites per block ---*/
	res = 0;
	block = -1;
	rewind(inFile);
	while ( fgets(line, 100, inFile) != NULL) {
		token = strtok(line, FIELD_TOKEN);
		if(0 == strcmp(token, "block")) {
			block++;
			if(block >= numBlocks) {
				fprintf(stdout, "*** Error: Unexpected error in second parse of site element file in processElementSiteFile().\n");
				res = -1;
				break;
			}
			data->numSitesPerBlock[block] = 0;
			data->siteData[block] = NULL;

			/*--- read through block id and then theta and lambda parameters ---*/
			token = strtok(NULL, FIELD_TOKEN);
			data->blockNames[block] = data->blockNames[0] + block*ID_LEN;
			strncpy(data->blockNames[block], token, ID_LEN-1);
			for(i=0; i<2; i++) {
				token = strtok(NULL, FIELD_TOKEN);
				if(0 == strcmp(token, "lambda")) {
					token = strtok(NULL, FIELD_TOKEN);
					res = sscanf(token, "%lf", &params->lambda[block]);
					if(res < 1) {
						fprintf(stdout, "*** Error: expecting to read in value for lambda parameter for block %d, got '%s', in processElementSiteFile().\n",
								block+1, token);
						res = -1;
						break;
					}
				} else if(0 == strcmp(token, "theta")) {
					token = strtok(NULL, FIELD_TOKEN);
					res = sscanf(token, "%lf", &params->theta[block]);
					if(res < 1) {
						fprintf(stdout, "*** Error: expecting to read in value for theta parameter for block %d, got '%s', in processElementSiteFile().\n",
								block+1, token);
						res = -1;
						break;
					}
				} else {
					fprintf(stdout, "*** Error: expecting only 'theta' and 'lambda' assignments in 'block' line of block %d, got '%s', in processElementSiteFile().\n",
							block+1, token);
					res = -1;
					break;
				}					
			}// end of for(i)
			if(res < 0) {
				break;
			}
			if(params->theta[block] < minThetaLambda) {
				params->theta[block] = minThetaLambda;
			}
			if(params->lambda[block] < minThetaLambda) {
				params->lambda[block] = minThetaLambda;
			}
		} else if(0 == strcmp(token, "site")) {
			if(block < 0) {
				fprintf(stdout, "*** Error: Unexpected error in second parse of site element file in processElementSiteFile().\n");
				res = -1;
				break;
			}
			data->numSitesPerBlock[block]++;
		}
	}// end of while(inFile)
	
	if(block != numBlocks-1) {
		fprintf(stdout, "*** Error: Unexpected error after second parse of site element file in processElementSiteFile().\n");
		res = -1;
	}

	if(res < 0.0) {
		freeModelParameters(params);
		freeCollectionData(data);
		fclose(inFile);		
		return -1;
	}
	
	/*--- third pass on file to get actual site data ---*/
	/*--- ignoring blocks with no lambda/theta (negative values) ---*/
	res = 0;
	block = -1;
	// ignoreBlock = 1;
	data->numSites = 0;
	rewind(inFile);
	site = 0;
	data->maxNumSamples = numSamples;
	while ( fgets(line, 100, inFile) != NULL) {
		token = strtok(line, FIELD_TOKEN);
		if(0 == strcmp(token, "block")) {
			block++;
			// if(!ignoreBlock) {
			// 	block++;
			// }
			if(block >= numBlocks) {
				fprintf(stdout, "*** Error: Unexpected error in second parse of site element file in processElementSiteFile().\n");
				res = -1;
				break;
			}
			/*--- ignore blocks with weak information on neutral poly/div ---*/
			if(params->theta[block] < 0.0 || params->lambda[block] < 0.0) {
				token = strtok(NULL, FIELD_TOKEN);
				printf("- Ignoring block '%s' because of uncertain local neutral parameters lambda and/or theta.\n",token);
				// ignoreBlock = 1;
				numBlocks--;
			} else {
				// ignoreBlock = 0;
				data->numSites += data->numSitesPerBlock[block];
				site = 0;
				data->siteData[block] = (SiteData*)malloc(data->numSitesPerBlock[block] * sizeof(SiteData));
				if(data->siteData[block] == NULL) {
					fprintf(stdout, "*** Error: OOM allocating site data for block %d in processElementSiteFile().\n", block+1);
					res = -1;
					break;
				}
/* DEBUG 	printf("Block %d with lambda %.8lf and %d sites.\n",block+1,params->lambda[block], data->numSitesPerBlock[block]);	*/	
			}
		} else if(0 == strcmp(token, "site")) {
			if(block < 0) {
				fprintf(stdout, "*** Error: Unexpected error in second parse of site element file in processElementSiteFile().\n");
				res = -1;
				break;
			}
			// if(ignoreBlock) {
			// 	continue;
			// }
			/*--- read through site id, then number of samples (if needed) ---*/
			token = strtok(NULL, FIELD_TOKEN);
			data->siteData[block][site].siteId[0] = '\0';
			strncpy(data->siteData[block][site].siteId,token,ID_LEN-1);
			if(numSamples <= 0) {
				token = strtok(NULL, FIELD_TOKEN);
				res = sscanf(token,"%d",&data->siteData[block][site].numSamples);
				if(res < 1 || data->siteData[block][site].numSamples <= 0) {
					fprintf(stdout, "*** Error: expecting positive number of chromosome samples for site %d in block %d, got '%s', in processElementSiteFile().\n",
							site+1, block+1, token);
					res = -1;
					break;
				}
				if(data->siteData[block][site].numSamples > data->maxNumSamples) {
					data->maxNumSamples = data->siteData[block][site].numSamples;
				}
			} else {
				data->siteData[block][site].numSamples = numSamples;
			}
			
			/*--- read site type ---*/
			token = strtok(NULL, FIELD_TOKEN);
			if(0 == strcmp(token,"M")) {
				data->siteData[block][site].polyClass  = MONOMORPHIC;
			} else if(0 == strcmp(token,"L")) {
				data->siteData[block][site].polyClass  = POLY_L;
			} else if(0 == strcmp(token,"H")) {
				data->siteData[block][site].polyClass  = POLY_H;
			} else {
				fprintf(stdout, "*** Error: illegal site type designation '%s' for site %d in block %d. Expecting 'M', 'L', or 'H', in processElementSiteFile().\n",
						token, site+1, block+1);
				res = -1;
				break;
			}				
			
			/*--- read Pr[Z = Xmaj] ---*/
			token = strtok(NULL, FIELD_TOKEN);
			res = sscanf(token,"%lf",&data->siteData[block][site].pZeqXmaj);
			if(res < 1 || data->siteData[block][site].pZeqXmaj < 0.0 || data->siteData[block][site].pZeqXmaj > 1.0) {
				fprintf(stdout, "*** Error: expecting probability for ancestral allele in [0,1] for site %d in block %d, got '%s', in processElementSiteFile().\n",
						site+1, block+1, token);
				res = -1;
				break;
			}
			/*--- read Pr[Z = Xmaj] ---*/
			if(data->siteData[block][site].polyClass == POLY_L || data->siteData[block][site].polyClass == POLY_H) {
				token = strtok(NULL, FIELD_TOKEN);
				res = sscanf(token,"%lf",&data->siteData[block][site].pZeqXmin);
				if(res < 1 || data->siteData[block][site].pZeqXmin < 0.0 || data->siteData[block][site].pZeqXmin > 1.0) {
					fprintf(stdout, "*** Error: expecting probability for ancestral allele in [0,1] for site %d in block %d, got '%s', in processElementSiteFile().\n",
							site+1, block+1, token);
					res = -1;
					break;
				}
				/*--- check posterior sum and consider rounding errors  ---*/
				if(data->siteData[block][site].pZeqXmaj + data->siteData[block][site].pZeqXmin - 1.0 > 1e-6) {
					fprintf(stdout, "*** Error: expecting sum of major and minor probability for ancestral allele to be in [0,1] for site %d in block %d, got %g+%g=%g, in processElementSiteFile().\n",
							site+1, block+1, data->siteData[block][site].pZeqXmaj, data->siteData[block][site].pZeqXmin,
							data->siteData[block][site].pZeqXmaj + data->siteData[block][site].pZeqXmin);
					res = -1;
					break;
				} else if(data->siteData[block][site].pZeqXmaj + data->siteData[block][site].pZeqXmin > 1.0) {
					data->siteData[block][site].pZeqXmaj /= (data->siteData[block][site].pZeqXmaj + data->siteData[block][site].pZeqXmin);
					data->siteData[block][site].pZeqXmin /= (data->siteData[block][site].pZeqXmaj + data->siteData[block][site].pZeqXmin);
				}
			} else {
				data->siteData[block][site].pZeqXmin = 0.0;
			}

			site++;
		}
	}// end of while(inFile)
	
	if(res < 0.0) {
		freeModelParameters(params);
		freeCollectionData(data);
		fclose(inFile);		
		return -1;
	}

	data->numBlocks = numBlocks;

	fclose(inFile);
	return 0;
}
/** end of processElementSiteFile **/



#define NUM_COUNTS 12
#define NUM_CHROMS 108
/** processCountFile
    Reads data from count file used in previous version of the model
    Initializes neutral parameters according to flank count
    Summarizes element data in CollectionData data structure
    @param countFileName a file name for element data
    @param data Empty structure to be initialized
    @param params Empty structure to be initialized
    @return 0 if all OK, and -1 otherwise
*/
int processCountFile(const char* countFileName, CollectionData* data, ModelParameters* params) {
	char line[101];
	char* linePointer;
	int block,numBlocks,site,i;
	
	double WattersonsArray[NUM_CHROMS+1];
	double wattersonsA;
	
	int numSites, numFlankLsites, numFlankHsites;
	
	double minThetaLambda = 0.00001;		// minimal theta/lambda, in case relevant counts are zero

//	double minTheta, minLambda; FIXME - GET RID OF THIS PART IF NOT USED
	
	/*--- data entries in each line of the count file ---*/
	int  lineData[NUM_COUNTS];
	int elemTotal   = 0;
	int flankTotal  = 1;
	/*--- UNUSED	int elemLHdiv   = 2; ---*/
	int elemHnodiv  = 3;
	int elemLnodiv  = 4;
	int elemMdiv    = 5;
	int elemMnodiv  = 6;
	/*--- UNUSED	int flankLHdiv  = 7; ---*/
	int flankHnodiv = 8;
	int flankLnodiv = 9;
	int flankMdiv   = 10;
	/*--- UNUSED	int flankMnodiv = 11; ---*/

	if(data == NULL || params == NULL) {
		fprintf(stdout, "*** Error: cannot initialize null CollectionData or ModelParameters in readCountFile().\n");
		return -1;
	}
	
	FILE* countFile = fopen (countFileName, "r");
	if(countFile == NULL) {
		fprintf(stdout, "*** Error: cannot open count file '%s' in readCountFile().\n",countFileName);
		return -1;
	}
	
	numBlocks = -1; // ignore header line
	while ( fgets(line, 100, countFile) != NULL) {
		numBlocks++;
	}
	
	/*--- allocate memory for data and parameters ---*/
	
	data->numSitesPerBlock = NULL;
	data->blockNames = NULL;
	data->siteData   = NULL;
	params->theta    = NULL;
	params->lambda   = NULL;
	
	data->numSitesPerBlock = (int*)malloc(numBlocks *  sizeof(int));
	data->blockNames = (char**)malloc(numBlocks *  sizeof(char*));
	if(data->blockNames != NULL) {
		data->blockNames[0] = NULL;
		data->blockNames[0] = (char*)malloc(numBlocks * ID_LEN * sizeof(char));
	}
	data->siteData = (SiteData**)malloc(numBlocks *  sizeof(SiteData*));
	params->theta  = (double*) malloc(numBlocks *  sizeof(double));
	params->lambda = (double*) malloc(numBlocks *  sizeof(double));
	
	if(data->siteData == NULL || data->numSitesPerBlock == NULL || data->blockNames == NULL || data->blockNames[0] == NULL || params->theta == NULL || params->lambda == NULL) {
		fprintf(stdout, "*** Error: OOM allocating arrays for CollectionData and ModelParameters in readCountFile().\n");
		freeModelParameters(params);
		freeCollectionData(data);
		fclose(countFile);		
		return -1;
	}

	
	for(block=0; block<numBlocks; block++) {
		data->numSitesPerBlock[block] = 0;
		data->siteData[block] = NULL;
	}// end of for(block)
	
	/*--- read data ---*/

	computeWattersons_a(NUM_CHROMS, WattersonsArray);
	wattersonsA = WattersonsArray[NUM_CHROMS];
	numFlankLsites = 0;
	numFlankHsites = 0;
	data->numSites = 0;
	data->maxNumSamples = NUM_CHROMS;
	data->numBlocks = numBlocks;
	rewind(countFile);
	linePointer = fgets(line, 100, countFile);
//	minTheta = minLambda = 2.0;  FIXME - GET RID OF THIS PART IF NOT USED
	for(block=0; block<numBlocks; block++) {
		if(NULL == fgets(line, 100, countFile)) {
			fprintf(stdout, "*** Error: expecting line for block %d in readCountFile().\n",block);
			freeModelParameters(params);
			freeCollectionData(data);
			fclose(countFile);
			return -1;
		}
		
		/*--- read and save block identifier ---*/
		data->blockNames[block] = data->blockNames[0] + block*ID_LEN;
		linePointer=strtok(line, FIELD_TOKEN);
		strncpy(data->blockNames[block], linePointer,ID_LEN-1);
		strcat(data->blockNames[block],":");
		linePointer=strtok(NULL,FIELD_TOKEN);
		strcat(data->blockNames[block], linePointer);
		strcat(data->blockNames[block],"-");
		linePointer=strtok(NULL,FIELD_TOKEN);
		strcat(data->blockNames[block], linePointer);

		linePointer = strtok(NULL,FIELD_TOKEN);
		for(i=0; i<NUM_COUNTS; i++) {
		
			if(linePointer == NULL) {
				fprintf(stdout, "*** Error: bad format for line in count file found in readCountFile().\n Line: %s\n",line);
				freeModelParameters(params);
				freeCollectionData(data);
				fclose(countFile);
				return -1;
			}
			lineData[i] = atoi (linePointer);
			linePointer = strtok(NULL,FIELD_TOKEN);
		}// end of for(i)
		
		/*--- ignore blocks with no data in elements or flanks ---*/
		if(lineData[flankTotal] == 0 || lineData[elemTotal] == 0) {
			block--;
			numBlocks--;
			continue;
		}
		
		/*--- process flank counts and estimate local theta and lambda ---*/
		params->theta[block] = (lineData[flankHnodiv] + lineData[flankLnodiv])/((double)lineData[flankTotal]);
		if(params->theta[block] < minThetaLambda) {
			params->theta[block] = minThetaLambda;
		}
		params->theta[block] /= wattersonsA;

		
		params->lambda[block] = lineData[flankMdiv]/((double)lineData[flankTotal]);
		if(params->lambda[block] < minThetaLambda) {
			params->lambda[block] = minThetaLambda;
		}

/*--- DEBUG
		printf("block %4d:",block+1);
		printf("  %4d flank poly sites, %4d monomorphic divergent sites and %4d flank sites total --> ",
					lineData[flankHnodiv] + lineData[flankLnodiv], lineData[flankMdiv], lineData[flankTotal]);

		printf("theta = %8lf , lambda = %6lf\n",params->theta[block],params->lambda[block]);
DEBUG ---*/
		
		numFlankLsites += lineData[flankLnodiv];
		numFlankHsites += lineData[flankHnodiv];
		
		/*--- process element counts ---*/

		numSites = lineData[elemTotal];
		data->numSitesPerBlock[block] = numSites;
		data->numSites += numSites;
		data->siteData[block] = (SiteData*)malloc(numSites * sizeof(SiteData));
		if(data->siteData[block] == NULL) {
			fprintf(stdout, "*** Error: OOM allocating site data for block %d in readCountFile().\n", block+1);
			freeModelParameters(params);
			freeCollectionData(data);
			fclose(countFile);
			return -1;
		}
		site=0;
		
		if(numSites != lineData[elemMnodiv] + lineData[elemMdiv] + lineData[elemLnodiv] + lineData[elemHnodiv]) {
			fprintf(stdout, "*** Error: block %d element counts inconsistent with total count %d in readCountFile().\n", block+1, numSites);
			freeModelParameters(params);
			freeCollectionData(data);
			fclose(countFile);
			return -1;
		}

		/*--- Monomorphic non divergent sites (M1) ---*/
		for(i=0; i<lineData[elemMnodiv]; i++, site++) {
			data->siteData[block][site].polyClass  = MONOMORPHIC;
			data->siteData[block][site].numSamples = NUM_CHROMS;
			data->siteData[block][site].pZeqXmaj   = 1.0;
			data->siteData[block][site].pZeqXmin   = 0.0;
			data->siteData[block][site].siteId[0] = '\0';
		}// end of for(i) - flankMnodiv

		/*--- Monomorphic divergent sites (M0) ---*/
		for(i=0; i<lineData[elemMdiv]; i++, site++) {
			data->siteData[block][site].polyClass  = MONOMORPHIC;
			data->siteData[block][site].numSamples = NUM_CHROMS;
			data->siteData[block][site].pZeqXmaj   = 0.0;
			data->siteData[block][site].pZeqXmin   = 0.0;
			data->siteData[block][site].siteId[0] = '\0';
		}// end of for(i) - flankMnodiv

		/*--- Polymorphic sites with low MAF/DAF (L1) ---*/
		for(i=0; i<lineData[elemLnodiv]; i++, site++) {
			data->siteData[block][site].polyClass  = POLY_L;
			data->siteData[block][site].numSamples = NUM_CHROMS;
			/*--- sites with low MAF where Z != Xmaj are masked or tagged as H sites ---*/
			data->siteData[block][site].pZeqXmaj   = 1.0;
			data->siteData[block][site].pZeqXmin   = 0.0;
			data->siteData[block][site].siteId[0] = '\0';
		}// end of for(i) - flankMnodiv

		/*--- Polymorphic sites with high MAF/DAF (H1) ---*/
		for(i=0; i<lineData[elemHnodiv]; i++, site++) {
			data->siteData[block][site].polyClass  = POLY_H;
			data->siteData[block][site].numSamples = NUM_CHROMS;
			/*--- the only thing that matters here is the sum of these two probabilities, which we know is 1.0 ---*/
			data->siteData[block][site].pZeqXmaj   = 1.0;
			data->siteData[block][site].pZeqXmin   = 0.0;
			data->siteData[block][site].siteId[0] = '\0';
		}// end of for(i) - flankMnodiv

	}// end of for(block)
	
/* FIXME - GET RID OF THIS PART IF NOT USED
	for(block=0; block<numBlocks; block++) {
		if(params->theta[block] == 0.0) {
			params->theta[block] = minTheta;
		}
		if(params->lambda[block] == 0.0) {
			params->lambda[block] = minLambda;
		}
	}// end of for(block)
*/	
	data->numBlocks = numBlocks;

	/*--- set beta (beta[2]=0 since for all L sites, the minor allele is derived) ---*/
	params->beta[1] = numFlankHsites / ((double)numFlankLsites + numFlankHsites);
	params->beta[0] = 1 - params->beta[1];
	params->beta[2] = 0.0;

/*--- DEBUG
	printf("\nTotal %d flank poly H sites,and %d total flank poly sites --> ",
				numFlankHsites, numFlankLsites + numFlankHsites);
	printf("beta[0] = %8lf ,beta[1] = %8lf ,beta[2] = %8lf\n",params->beta[0],params->beta[1],params->beta[2]);
DEBUG ---*/

	fclose(countFile);
	return 0;
}
/** end of processCountFile **/



/** selectionEM
    Performs parameter inference through EM for parameters modeling selection - rho, eta, gamma.
    @param elementData Data summary for functional elements of interest
    @param params Model parameters. 
    		Neutral parameters (beta, theta, lambda) are kept constant.
    		Selection parameters (rho, eta, gamma) are used for initialization, and then updated throughout the procedure.
    @param emSetup struct with various EM setup arguments 
    @param emStatus used to return EM status (likelihood and other information)
	@return  0 if didn't reach convergence
	@return  1 if successfully converged
	@return  2 if converged but overshot (last update reduced likelihood)
	@return -1 if found error
*/
int selectionEM(CollectionData* elementData, ModelParameters* params, EMsetup* emSetup, EMstatus* emStatus) {
	
// DEBUG POSTERIORS	double cSel, cNeut;
	
	/*------------------------------------------------------------------------*/
	/*----------------------------   VARIABLES    ----------------------------*/
	/*------------------------------------------------------------------------*/
	
	FILE*  emLogFile        = NULL;
	char*  emLogFileName    = emSetup->logFileName;
	int    maxNumIterations = emSetup->maxNumIterations;
	int    logInterval      = emSetup->logInterval;
	double lnLd_epsilon     = emSetup->lnLd_epsilon;
	unsigned short fixRho   = emSetup->fixRho;
	unsigned short fixEta   = emSetup->fixEta;
	unsigned short fixGamma = emSetup->fixGamma;
	
	/*--- model parameters ---*/
	int numSites = elementData->numSites;
	int numBlocks = elementData->numBlocks;
	/*--- pre-estimated neutral params ---*/
	double* lambda = params->lambda;
	double* theta  = params->theta;
	double* beta   = params->beta;
	/*--- initial values for selection params ---*/
	double	rho    = emSetup->rhoInitial;
	double	eta    = emSetup->etaInitial;
	double	gamma  = emSetup->gammaInitial;

	double* WattersonsArray = NULL;		// stores values of Watterson's a
	
	/*--- arrays to store counts ---*/
	double* allCounts = NULL;			// used for memory allocation
	double* rhoCounts = NULL;			// stores counts for rho
	double* etaCounts = NULL;			// stores counts for eta
	double* gammaCounts = NULL;			// stores counts for gamma
	
	/*--- number of sites in the dataset in each polymorphism class and block---*/
	int* numSitesPerBlock	= elementData->numSitesPerBlock;	// total number of sites per block
	int* siteBoundaryArray	= NULL;		// used for memory allocation
	int* firstPatternInBlock	= NULL;	// index of first site pattern in block
	int* firstLpatternInBlock	= NULL;	// index of first L site pattern in block (boundary between L and M patterns)
	
	
	/*--- separately handle sites that cannot be under selection in the model (H sites in our case) ---*/
	/*--- pre-compute their contribution to the log likelihood function ---*/
	int numNeutralSites;				// number of restricted neutral sites
	int neutrals;						// number of restricted neutral sites in current block
	double lnLd_neutralSites;			// total log likelihood of restricted neutral sites (conditioning on them being neutral) 
	
	/*--- auxilliary values maintained for each site and each possible count ---*/
	/*--- these values are computed before the EM loop and do not change throughout the algorithm ---*/
	AuxStatPattern_SelectionEM* auxStats = NULL;	// used for memoray allocation

	/*--- for nmerical optimization of eta and gamma in M step ---*/
	SumLogsArguments gammaFunctionArgs;
	SumLogsArguments etaFunctionArgs;
	double* factorsArray = NULL;		// used for memory allocation
	double* gammaFactors = NULL;		// theta_b * average a_i across all sites in block (non high polymorphic)
	double* etaFactors   = NULL;		// lambda_b

	/*--- indexing variables ---*/
	int block, site;
	
	/*--- used in initialization ---*/
	double pZeqXmaj, pZeqXmin, wattersonsA;
	
	/*--- used in EM ---*/
	int res;		
	int emIteration;
	int sitePattern, numPatterns;
	double lnLd, lnLd_prev, lnLd_diff;
	double expectedLnLd, expectedLnLd_prev, expectedLnLd_diff;
	double expectedLnLd_rho, expectedLnLd_eta, expectedLnLd_gamma;
	double neutFactor, noDivFactor, divFactor;
	double qNeut, qNoDiv, qDiv, qSum;
	double pNeut, pNoDiv, pDiv, pSel;
	
	unsigned short zeroLikelihood;		// flag indicating that site with zero likelihood have been encountered
	
	
	/*-------------------------------------------------------------------------*/
	/*------------------------   ALLOCATION - START   -------------------------*/
	/*-------------------------------------------------------------------------*/
	
	WattersonsArray = (double*)malloc((elementData->maxNumSamples+1) * sizeof(double));
	allCounts       = (double*)malloc((2+(numBlocks+1)+(numBlocks+1)) * sizeof(double));
	auxStats        = (AuxStatPattern_SelectionEM*)malloc((numSites) * sizeof(AuxStatPattern_SelectionEM));
	factorsArray    = (double*)malloc((2*numBlocks) * sizeof(double));
	siteBoundaryArray  = (int*)malloc((2*numBlocks+1) * sizeof(int));
	
	if(WattersonsArray == NULL || allCounts == NULL || auxStats == NULL || factorsArray == NULL || siteBoundaryArray == NULL) {
		fprintf(stdout, "*** Error: OOM allocating space for arrays in selectionEM.\n");
		freeCheck(WattersonsArray);
		freeCheck(allCounts);
		freeCheck(auxStats);
		freeCheck(factorsArray);
		freeCheck(siteBoundaryArray);
		return -1;
	}

	rhoCounts   = allCounts;
	etaCounts   = rhoCounts + 2;
	gammaCounts = etaCounts + numBlocks+1;
	
	gammaFactors = factorsArray;
	etaFactors   = factorsArray + numBlocks;

	firstPatternInBlock = siteBoundaryArray;
	firstLpatternInBlock = siteBoundaryArray + numBlocks + 1;
	
	gammaFunctionArgs.numTerms     = numBlocks;
	gammaFunctionArgs.coefficients = gammaCounts;
	gammaFunctionArgs.factors      = gammaFactors;

	etaFunctionArgs.numTerms       = numBlocks;
	etaFunctionArgs.coefficients   = etaCounts;
	etaFunctionArgs.factors        = etaFactors;

	/*-----------------------------------------------------------------------*/
	/*------------------------   ALLOCATION - END   -------------------------*/
	/*-----------------------------------------------------------------------*/
	
	
	
	/*-----------------------------------------------------------------------------*/
	/*------------------------   INITIALIZATION - START  --------------------------*/
	/*-----------------------------------------------------------------------------*/
	
	/*--- pre-compute Watterson's a values for all possible numbers of samples ---*/
	computeWattersons_a(elementData->maxNumSamples, WattersonsArray);

	/*---------------------------------------------------*/
	/*---   INIT RESTRICTED NEUTRAL SITES (TYPE H)    ---*/
	/*---------------------------------------------------*/
	
	/*--- numNeutralSites holds the number of H sites - that cannot be under selection in the model  ---*/
	/*--- lnLd_neutralSites holds the contribution of these sites to the log likelihood function (conditioning on them being neutral)  ---*/
	numNeutralSites = 0;
	lnLd_neutralSites = 0.0;
	for(block=0; block<numBlocks; block++) {
		neutrals = 0;
		for(site=0; site<numSitesPerBlock[block]; site++) {
			if(POLY_H == elementData->siteData[block][site].polyClass) {
				wattersonsA = WattersonsArray[ elementData->siteData[block][site].numSamples ];
				pZeqXmaj    = elementData->siteData[block][site].pZeqXmaj;
				pZeqXmin    = elementData->siteData[block][site].pZeqXmin;

				if(emSetup->useSimplifiedModel) {
					qNeut = (1.0-lambda[block]) * beta[1]*theta[block]*wattersonsA;
				} else {
					qNeut = ((1.0-2.0*lambda[block]/3.0) * (pZeqXmaj+pZeqXmin) + (2.0*lambda[block]/3.0) * (1-pZeqXmaj-pZeqXmin) )
							* beta[1]*theta[block]*wattersonsA/3.0;
				}
				lnLd_neutralSites   += customLog(qNeut);
				
				neutrals++;
			}// end of if(POLY_H)
		}// end of for(site)
		numNeutralSites += neutrals;
	}// end of for(blocks) - neutral site initialization
	
	if(numNeutralSites == 0) {
		printf("*** Warning: no a-priori neutral sites in elements!\n");
	}
	/*----------------------------------------------------*/
	/*---   INIT AUXILLIARY EXPRESSIONS (TYPES M,L)    ---*/
	/*----------------------------------------------------*/
	
	/*--- group sites according to block, type (M/L), and site pattern  ---*/
	/*--- first process M sites and then L sites  ---*/
	numPatterns  = 0;
	firstPatternInBlock[0] = 0;
	for(block=0; block<numBlocks; block++) {

		/*-----------------------------*/
		/*---------- M sites ----------*/
		/*-----------------------------*/

		/*--- first identify identical patterns by recording P[Z = Xmaj] in auxStats[].noDiv ---*/
		for(site=0; site<numSitesPerBlock[block]; site++) {
			if(MONOMORPHIC == elementData->siteData[block][site].polyClass) {
				wattersonsA = WattersonsArray[ elementData->siteData[block][site].numSamples ];
				pZeqXmaj    = elementData->siteData[block][site].pZeqXmaj;
				
				for(sitePattern=firstPatternInBlock[block]; sitePattern<numPatterns; sitePattern++) {
					if(wattersonsA == auxStats[sitePattern].wattA && pZeqXmaj == auxStats[sitePattern].noDiv) {
						auxStats[sitePattern].numSites++;
						break;
					}
				}// end of for(sitePattern)
				
				if(sitePattern == numPatterns) {
					auxStats[numPatterns].wattA = wattersonsA;
					auxStats[numPatterns].noDiv = pZeqXmaj;
					auxStats[numPatterns].numSites = 1;
					numPatterns++;
				}// end of if(sitePattern == numPatterns)
				
			}// end of if(MONOMORPHIC)
		}// end of for(site)

		/*--- now compute stats for all distinct site patterns ---*/
		for(sitePattern=firstPatternInBlock[block]; sitePattern<numPatterns; sitePattern++) {
			wattersonsA = auxStats[sitePattern].wattA;
			pZeqXmaj    = auxStats[sitePattern].noDiv;
			
			if(emSetup->useSimplifiedModel) {
				auxStats[sitePattern].neut  = 
					(1-theta[block]*auxStats[sitePattern].wattA) * ( (1-lambda[block]) * pZeqXmaj + lambda[block] * (1-pZeqXmaj) );
			} else {
				auxStats[sitePattern].neut  = 
					(1-theta[block]*auxStats[sitePattern].wattA) * ( (1-lambda[block]) * pZeqXmaj + lambda[block]/3 * (1-pZeqXmaj) );
			}
			auxStats[sitePattern].noDiv = pZeqXmaj;
			if(emSetup->useSimplifiedModel) {
				auxStats[sitePattern].div   = lambda[block] * (1-pZeqXmaj);
			} else {
				auxStats[sitePattern].div   = lambda[block]/3 * (1-pZeqXmaj);
			}
		}// end of for(sitePattern)


		/*-----------------------------*/
		/*---------- L sites ----------*/
		/*-----------------------------*/

		/*--- first identify identical patterns by recording P[Z = Xmaj] in auxStats[].noDiv and P[Z = Xmin] in auxStats[].div ---*/
		firstLpatternInBlock[block] = numPatterns;
		for(site=0; site<numSitesPerBlock[block]; site++) {
			if(POLY_L == elementData->siteData[block][site].polyClass) {
				wattersonsA = WattersonsArray[ elementData->siteData[block][site].numSamples ];
				pZeqXmaj    = elementData->siteData[block][site].pZeqXmaj;
				pZeqXmin    = elementData->siteData[block][site].pZeqXmin;
				
				for(sitePattern=firstLpatternInBlock[block]; sitePattern<numPatterns; sitePattern++) {
					if(wattersonsA == auxStats[sitePattern].wattA && pZeqXmaj == auxStats[sitePattern].noDiv && pZeqXmin == auxStats[sitePattern].div) {
						auxStats[sitePattern].numSites++;
						break;
					}
				}// end of for(sitePattern)
				
				if(sitePattern == numPatterns) {
					auxStats[numPatterns].wattA = wattersonsA;
					auxStats[numPatterns].noDiv = pZeqXmaj;
					auxStats[numPatterns].div   = pZeqXmin;
					auxStats[numPatterns].numSites = 1;
					numPatterns++;
				}// end of if(sitePattern == numPatterns)
				
			}// end of if(POLY_L)
		}// end of for(site)

		/*--- now compute stats for all distinct site patterns ---*/
		for(sitePattern=firstLpatternInBlock[block]; sitePattern<numPatterns; sitePattern++) {
			wattersonsA = auxStats[sitePattern].wattA;
			pZeqXmaj    = auxStats[sitePattern].noDiv;
			pZeqXmin    = auxStats[sitePattern].div;
			
			if(emSetup->useSimplifiedModel) {
				auxStats[sitePattern].neut   = 
					(1-lambda[block]) * beta[0] * theta[block]* wattersonsA;
			} else {
				auxStats[sitePattern].neut   = 
					 lambda[block]/3  * (beta[0] * (1-pZeqXmaj) + beta[2] * (1-pZeqXmin) ) * theta[block]* wattersonsA / 3 +
					(1-lambda[block]) * (beta[0] *   pZeqXmaj   + beta[2] *   pZeqXmin   ) * theta[block]* wattersonsA / 3;
			}
			if(emSetup->useSimplifiedModel) {
				auxStats[sitePattern].noDiv  =  theta[block]* wattersonsA;
			} else {
				auxStats[sitePattern].noDiv  =  pZeqXmaj * theta[block]* wattersonsA / 3;
			}
			auxStats[sitePattern].div    =  0.0;
		}// end of for(sitePattern)


		/*--- set boundary with next block (or end of pattern list for last block) ---*/
		firstPatternInBlock[block+1] = numPatterns;
	}// end of for(block) - L/M site pattern initialization

/*--- DEBUG ---*
	printf("beta parameters: beta1=%.4lf, beta2=%.4lf, beta3=%.4lf\n", beta[0],beta[1],beta[2]);
	for(block=0; block<numBlocks; block++) {
		printf("block %d, theta %.4lf lambda %.4lf,  ", block+1, theta[block], lambda[block]);
		if(firstLpatternInBlock[block]-firstPatternInBlock[block] > 0) {
			printf(" M site patterns:");
			for(sitePattern=firstPatternInBlock[block]; sitePattern<firstLpatternInBlock[block]; sitePattern++) {
				printf(" %d sites with (%.4lf,%.4lf,%.4lf) a=%.4lf", 
					   auxStats[sitePattern].numSites, auxStats[sitePattern].neut, auxStats[sitePattern].noDiv, auxStats[sitePattern].div, auxStats[sitePattern].wattA);
			}// end of for(sitePattern)
		}
		if(firstPatternInBlock[block+1]-firstLpatternInBlock[block] > 0) {
			printf(" L site patterns:");
			for(sitePattern=firstLpatternInBlock[block]; sitePattern<firstPatternInBlock[block+1]; sitePattern++) {
				printf(" %d sites with (%.4lf,%.4lf,%.4lf) a=%.4lf", 
					   auxStats[sitePattern].numSites, auxStats[sitePattern].neut, auxStats[sitePattern].noDiv, auxStats[sitePattern].div, auxStats[sitePattern].wattA);
			}// end of for(sitePattern)
		}
		neutrals = 0;
		for(site=0; site< numSitesPerBlock[block]; site++) {
			if(POLY_H == elementData->siteData[block][site].polyClass) {
				neutrals++;
			}// end of if(POLY_H)
		}// end of for(site)
		printf(" %d H1 sites", neutrals);
		printf("\n");
	}// end of for(block) - DEBUG
*--- DEBUG ---*/
/*--- DEBUG OLD ---*
	printf("beta parameters: beta1=%.4lf, beta2=%.4lf, beta3=%.4lf\n", beta[0],beta[1],beta[2]);
	for(block=0; block<numBlocks; block++) {
		printf("block %d, theta %.4lf lambda %.4lf,  ", block+1, theta[block], lambda[block]);
		if(firstLpatternInBlock[block]-firstPatternInBlock[block] > 0) {
			sitePattern = firstPatternInBlock[block];
			if(auxStats[sitePattern].noDiv > 0.0) {
				printf(" %d M1 sites", auxStats[sitePattern].numSites);
				if(firstLpatternInBlock[block]-firstPatternInBlock[block] > 1) {
					sitePattern = firstPatternInBlock[block]+1;
					printf(" %d M0 sites", auxStats[sitePattern].numSites);
					if(firstLpatternInBlock[block]-firstPatternInBlock[block] > 2) {
						printf(" Error: found %d (>2) distinct monomprphic site patterns found for block %d.", firstLpatternInBlock[block]-firstPatternInBlock[block], block+1);
					}
				} else {
					printf(" %d M0 sites", 0);
				}
			} else {
				printf(" %d M1 sites %d M0 sites", 0, auxStats[sitePattern].numSites);
			}
		} else {
			printf(" %d M1 sites %d M0 sites", 0, 0);
		}
		if(firstPatternInBlock[block+1]-firstLpatternInBlock[block] > 0) {
			sitePattern = firstLpatternInBlock[block];
			printf(" %d L1 sites", auxStats[sitePattern].numSites);
			if(firstPatternInBlock[block+1]-firstLpatternInBlock[block] > 1) {
				printf(" Error: found %d (>1) distinct L site patterns found for block %d.", firstPatternInBlock[block+1]-firstLpatternInBlock[block], block+1);
			}
		} else {
			printf(" %d L1 sites", 0);
		}
		neutrals = 0;
		for(site=0; site< numSitesPerBlock[block]; site++) {
			if(POLY_H == elementData->siteData[block][site].polyClass) {
				neutrals++;
			}// end of if(POLY_H)
		}// end of for(site)
		printf(" %d H1 sites", neutrals);
		printf("\n");
	}// end of for(block) - DEBUG
		
*--- DEBUG OLD ---*/


	/*------------------------------------------*/
	/*---   INIT ETA/GAMMA FACTOR VECTORS    ---*/
	/*------------------------------------------*/
	gammaFunctionArgs.minFactor = 1000.0;
	etaFunctionArgs.minFactor   = 1000.0;
	gammaFunctionArgs.maxFactor = -1.0;
	etaFunctionArgs.maxFactor   = -1.0;
	for(block=0; block<numBlocks; block++) {
//		printf("theta %g, lambda %g wattA %g", theta[block], lambda[block], auxStats[ firstPatternInBlock[block] ].wattA);
		if(firstPatternInBlock[block+1] - firstPatternInBlock[block] > 0) {
			etaFactors[block] = lambda[block];
			if(etaFactors[block] > etaFunctionArgs.maxFactor) {
//				printf("eta max factor updated at block %d to %g.\n", block+1, etaFactors[block]);
				etaFunctionArgs.maxFactor = etaFactors[block];
//				printf(" lambda_max");
			}
			if(etaFactors[block] < etaFunctionArgs.minFactor) {
//				printf("eta min factor updated at block %d to %g.\n", block+1, etaFactors[block]);
				etaFunctionArgs.minFactor = etaFactors[block];
//				printf(" lambda_min");
			}
		} else {
			/*--- blocks with no M/L sites --> etaCounts[block] = 0  ---*/
//			printf(" no eta counts");
			etaFactors[block] = 0.0;
		}
		
		/*--- we now assume Watterson's a is constant across each block. Take value from first M site pattern  ---*/ 
		/*--- FIXME - fix to individual gamma factor per site ---*/
		if(firstLpatternInBlock[block] - firstPatternInBlock[block] > 0) {
			gammaFactors[block] = theta[block] * auxStats[ firstPatternInBlock[block] ].wattA;
			if(gammaFactors[block] > gammaFunctionArgs.maxFactor) {
/*--- DEBUG		printf("gamma max factor updated at block %d to %g.\n", block+1, gammaFactors[block]); ---*/
				gammaFunctionArgs.maxFactor = gammaFactors[block];
/*--- DEBUG		printf(" theta_max"); ---*/
			}
			if(gammaFactors[block] < gammaFunctionArgs.minFactor) {
/*--- DEBUG		printf("gamma min factor updated at block %d to %g.\n", block+1, gammaFactors[block]); ---*/
				gammaFunctionArgs.minFactor = gammaFactors[block];
/*--- DEBUG		printf(" theta_min"); ---*/
			}
		} else {
			/*--- blocks with no M sites   --> gammaCounts[block] = 0 ---*/
/*--- DEBUG		printf(" no gamma counts"); ---*/
			gammaFactors[block] = 0.0;
		}
/*--- DEBUG		printf("\n"); ---*/
	}// end of for(block)

	/*--- set inverse of max/min factors  ---*/ 
	etaFunctionArgs.maxFactorInverse = 1/etaFunctionArgs.maxFactor;
	etaFunctionArgs.minFactorInverse = 1/etaFunctionArgs.minFactor;
	gammaFunctionArgs.maxFactorInverse = 1/gammaFunctionArgs.maxFactor;
	gammaFunctionArgs.minFactorInverse = 1/gammaFunctionArgs.minFactor;
	
/*--- DEBUG		printf("gamma factors in [%g,%g], eta factors in [%g,%g].\n",
					gammaFunctionArgs.minFactor, gammaFunctionArgs.maxFactor,etaFunctionArgs.minFactor, etaFunctionArgs.maxFactor);  ---*/

	/*----------------------------------------------------------------------------*/
	/*----------------------------   INIT LOG FILE    ----------------------------*/
	/*----------------------------------------------------------------------------*/
	if(emLogFileName[0] != '\0') {
		emLogFile = fopen(emLogFileName,"w");
		if(emLogFile == NULL) {
			fprintf(stdout, "*** Error: could not open log file '%s' in selectionEM.\n", emLogFileName);
			freeCheck(WattersonsArray);
			freeCheck(allCounts);
			freeCheck(auxStats);
			freeCheck(factorsArray);
			freeCheck(siteBoundaryArray);
			return -1;
		}
		fprintf(emLogFile, 	"%7s %8s %8s %8s %10s %10s %10s %10s\n", 
							"iter", "rho", "eta", "gamma", "lnLd", "lnLd_df", "E[lnLd]", "E[lnLd]_df");
		fflush(emLogFile);
	}

	
	/*---------------------------------------------------------------------------*/
	/*------------------------   INITIALIZATION - END  --------------------------*/
	/*---------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------*/
	/*----------------------------   EM LOOP    ----------------------------*/
	/*----------------------------------------------------------------------*/
	
	emIteration = 0;
	zeroLikelihood = 0;
	lnLd = lnLd_prev = lnLd_diff = 0.0;
	expectedLnLd = expectedLnLd_prev = 0.0;
	strcpy(emStatus->endState,"running");
	printf("Progress: ");
	fflush(stdout);
	while(1) {
		
		/*-------------------------------------*/
		/*---------  E STEP - START  ----------*/
		/*-------------------------------------*/

		/*--- init counts and likelihood ---*/
		lnLd_prev = lnLd;
		/*--- contribution from sites restricted to be neutral ---*/
		/*--- check for rho == 1.0 ---*/
		if(rho < 1.0) {
			lnLd = lnLd_neutralSites + customLog(1-rho)*numNeutralSites;
		} else if(numNeutralSites == 0) {
			lnLd = 0.0;
		} else {
			lnLd = -DBL_MAX;
			strcpy(emStatus->endState,"zero-likelihood");
			break;
		}
		
// DEBUG POSTERIORS		cSel = 0.0;
// DEBUG POSTERIORS		cNeut = numNeutralSites;
		
		rhoCounts[0] = numNeutralSites;
		rhoCounts[1] = 0.0;
		
		for(block=0; block<=numBlocks; block++) {
			etaCounts[block]   = 0.0;
			gammaCounts[block] = 0.0;
		}// end of for(block)

		/*--- loop through all site patterns ---*/
		neutFactor  = 1-rho;
		for(block=0; block<numBlocks; block++) {
			
			/*-----------------------------*/
			/*---------- M sites ----------*/
			/*-----------------------------*/
			
			/*** neutFactor  = 1-rho;  <-- computed before looping on blocks ***/
			/*--- FIXME - this assumes constant Watterson's a in block ---*/
			noDivFactor = rho * (1 - eta*etaFactors[block]) * (1 - gamma*gammaFactors[block]);
			divFactor   = rho * eta;
			for(sitePattern=firstPatternInBlock[block]; sitePattern<firstLpatternInBlock[block]; sitePattern++) {
				/*--- likelihood components ---*/
				qNeut  = neutFactor  * auxStats[sitePattern].neut;
				qNoDiv = noDivFactor * auxStats[sitePattern].noDiv;
				qDiv   = divFactor   * auxStats[sitePattern].div;
				
				/*--- contribution of site pattern to likelihood ---*/
				qSum   = qNeut + qNoDiv + qDiv;
				if(qSum == 0.0) {
					lnLd = -DBL_MAX;
					strcpy(emStatus->endState,"zero-likelihood");
					zeroLikelihood = 1;
					break;
				}
				lnLd += customLog(qSum) * auxStats[sitePattern].numSites;

				/*--- contribution of site pattern to counts ---*/
				pNeut  = qNeut / qSum;
				pSel   = 1-pNeut;
				pNoDiv = qNoDiv / qSum;
				pDiv   = qDiv / qSum;
				/*--- multiply by number of sites ---*/
				pNeut  *= auxStats[sitePattern].numSites;
				pSel   *= auxStats[sitePattern].numSites;
				pNoDiv *= auxStats[sitePattern].numSites;
				pDiv   *= auxStats[sitePattern].numSites;
				/*--- update counts ---*/
				rhoCounts[0]         += pNeut;
				rhoCounts[1]         += pSel;
				etaCounts[block]     += pNoDiv;
				etaCounts[numBlocks] += pDiv;
				gammaCounts[block]   += pNoDiv;
				
// DEBUG POSTERIORS				cNeut += pNeut;
// DEBUG POSTERIORS				cSel  += pSel;
			}// end of for(sitePattern) - M sites
			
			if(zeroLikelihood) {
				break;
			}
			
			/*-----------------------------*/
			/*---------- L sites ----------*/
			/*-----------------------------*/
			
			/*** neutFactor  = 1-rho;  <-- computed before looping on blocks ***/
			noDivFactor = rho * (1 - eta*etaFactors[block]) * gamma;
			divFactor   = 0.0;		// <-- not really used
			for(sitePattern=firstLpatternInBlock[block]; sitePattern<firstPatternInBlock[block+1]; sitePattern++) {
				/*--- likelihood components ---*/
				qNeut  = neutFactor  * auxStats[sitePattern].neut;
				qNoDiv = noDivFactor * auxStats[sitePattern].noDiv;
				qDiv   = 0.0;		// <-- not really used
				
				/*--- contribution of site pattern to likelihood ---*/
				qSum   = qNeut + qNoDiv; /*** + qDiv;  <-- not really used ***/
				if(qSum == 0.0) {
					lnLd = -DBL_MAX;
					strcpy(emStatus->endState,"zero-likelihood");
					break;
				}
					
				lnLd += customLog(qSum) * auxStats[sitePattern].numSites;

				/*--- contribution of site pattern to counts ---*/
				pNeut  = qNeut / qSum;
				pSel   = 1-pNeut;
				/*--- multiply by number of sites ---*/
				pNeut  *= auxStats[sitePattern].numSites;
				pSel   *= auxStats[sitePattern].numSites;
				/*--- update counts ---*/
				rhoCounts[0]           += pNeut;
				rhoCounts[1]           += pSel;
				etaCounts[block]       += pSel;
				gammaCounts[numBlocks] += pSel;
				
// DEBUG POSTERIORS				cNeut += pNeut;
// DEBUG POSTERIORS				cSel  += pSel;
			}// end of for(sitePattern) - L sites
			
			if(zeroLikelihood) {
				break;
			}

/*--- DEBUG		printf("etaCounts[%d] = %g, gammaCounts[%d] = %g\n", block+1, etaCounts[block], block+1, gammaCounts[block]); ---*/
			
		}// end of for(block) - E step and current log-likelihood computation
		
/*--- DEBUG		printf("etaCounts[%d] = %g, gammaCounts[%d] = %g\n", block+1, etaCounts[block], block+1, gammaCounts[block]); ---*/
		/*-----------------------------------*/
		/*---------  E STEP - END  ----------*/
		/*-----------------------------------*/

		/*---------------------------------------------*/
		/*--------------   LOG PROGRESS  --------------*/
		/*---------------------------------------------*/

		lnLd_diff = lnLd - lnLd_prev;
		expectedLnLd_diff = expectedLnLd - expectedLnLd_prev;
		emIteration++;
		if (emIteration % logInterval == 0) {
			printf(".");
			if(emLogFile != NULL) {
				fprintf(emLogFile, "%7d %7lf %7g %7lf %10lf %10lf %10lf %10lf\n",
						emIteration, rho, eta, gamma, lnLd, lnLd_diff, expectedLnLd, expectedLnLd_diff);
				fflush(emLogFile);
			}
			if(emIteration%(logInterval*10) == 0) {
				printf(" ");
				if(emIteration%(logInterval*100) == 0) {
					printf("\n");
					printf("          ");
				}
			}
			fflush(stdout);
		}

		/*--------------------------------------------------------*/
		/*--------------   UPDATE STATUS AND BREAK  --------------*/
		/*--------------------------------------------------------*/

		if(zeroLikelihood) {
			break;
		}
		
		if(emIteration > 1) {
			if(lnLd_diff < lnLd_epsilon) {
				if(lnLd_diff < 0.0) {
					strcpy(emStatus->endState,"overshoot");
				} else {
					strcpy(emStatus->endState,"converged");
				}
				break;
			} else if(emIteration >= maxNumIterations){
				strcpy(emStatus->endState,"timeout");
				break;
			}
		}// end of if(emIteration > 0)


		/*-------------------------------------*/
		/*---------  M STEP - START  ----------*/
		/*-------------------------------------*/

		/*-------------------------------------------------------------------*/
		/*--------------   EXPECTED LOG LIKELIHOOD PRE UPDATE  --------------*/
		/*-------------------------------------------------------------------*/
		
		expectedLnLd_rho   = 0.0;
		if(rho < 1.0) {
			expectedLnLd_rho += rhoCounts[0]*customLog(1-rho);
		} else if(rhoCounts[0] > 0.0) {
			expectedLnLd_rho = -DBL_MAX;
		}
		if(rho > 0.0) {
			expectedLnLd_rho += rhoCounts[1]*customLog(rho);
		} else if(rhoCounts[1] > 0.0) {
			printf("\n\nError: rhoCounts[1] = %g with rho = %g.\n",rhoCounts[1],rho);
			expectedLnLd_rho = -DBL_MAX;
			strcpy(emStatus->endState,"error");
			break;
		}
		expectedLnLd_eta   = evalSumLogs(eta, &etaFunctionArgs);
		expectedLnLd_gamma = evalSumLogs(gamma, &gammaFunctionArgs);
		expectedLnLd_prev  = expectedLnLd_rho + expectedLnLd_eta + expectedLnLd_gamma;

		/*---------------------------------------------------*/
		/*--------------   UPDATE PARAMETERS   --------------*/
		/*---------------------------------------------------*/
		
		/*--- rho ---*/
		if(!fixRho) {
			rho = rhoCounts[1] / (rhoCounts[0]+rhoCounts[1]);
		}

		/*--- eta ---*/
/*--- DEBUG		printf("maximizing eta (%g): ", eta); ---*/
		if(!fixEta) {
			res = maximizeSumLogs(&etaFunctionArgs, &eta, &expectedLnLd_eta, OPT_METHOD);
			if(res < 0) {
				fprintf(stdout, "\n*** Error: when maximizing eta (%g) in selectionEM iteration %d.\n", eta, emIteration);
				strcpy(emStatus->endState,"error");
				break;
			} else if(res == 1) {
//				fprintf(stdout, "Could not find better value for eta in selectionEM iteration %d (rho = %g , eta = %g , gamma = %g ).\n", emIteration, rho, eta,gamma);
			}
				
		}

		/*--- gamma ---*/
/*--- DEBUG		printf("maximizing gamma (%g): ", gamma); ---*/
		if(!fixGamma) {
			res = maximizeSumLogs(&gammaFunctionArgs, &gamma, &expectedLnLd_gamma, OPT_METHOD);
			if(res < 0) {
				fprintf(stdout, "\n*** Error: when maximizing gamma (%g) in selectionEM iteration %d.\n", gamma, emIteration);
				strcpy(emStatus->endState,"error");
				break;
			} else if(res == 1) {
//				fprintf(stdout, "Could not find better value for gamma in selectionEM iteration %d (rho = %g , eta = %g , gamma = %g ).\n", emIteration, rho, eta,gamma);
			}
		}

		/*--------------------------------------------------------------------*/
		/*--------------   EXPECTED LOG LIKELIHOOD POST UPDATE  --------------*/
		/*--------------------------------------------------------------------*/
		
		expectedLnLd_rho = 0.0;
		if(rho < 1.0) {
			expectedLnLd_rho += rhoCounts[0]*customLog(1-rho);
		} else if(rhoCounts[0] > 0.0) {
			expectedLnLd_rho = -DBL_MAX;
		}
		if(rho > 0.0) {
			expectedLnLd_rho += rhoCounts[1]*customLog(rho);
		} else if(rhoCounts[1] > 0.0) {
			printf("\n*** Error: rhoCounts[1] = %g with rho = %g.\n",rhoCounts[1],rho);
			expectedLnLd_rho = -DBL_MAX;
			strcpy(emStatus->endState,"error");
			break;
		}

		//Note: expectedLnLd_eta  and expectedLnLd_gamma are evaluated during maximization
		expectedLnLd = expectedLnLd_rho + expectedLnLd_eta + expectedLnLd_gamma;

		/*-----------------------------------*/
		/*---------  M STEP - END  ----------*/
		/*-----------------------------------*/
		
	} // end of while(1)  -- EM update loop
	
	
	/*-----------------------------------------------------------------------*/
	/*----------------------------   FINALIZE    ----------------------------*/
	/*-----------------------------------------------------------------------*/
	
	printf("\n");
	
// DEBUG POSTERIORS			printf("cSel = %f, cNeut = %f, ratio = %f.\n", cSel, cNeut, cSel/(cSel+cNeut));

	/*--- write inferred parameters into params and return---*/
	params->rho   = rho;
	params->eta   = eta;
	params->gamma = gamma;
	
	emStatus->numIterations = emIteration;
	emStatus->lnLikelihood  = lnLd;
	emStatus->lnLd_diff     = lnLd_diff;

	if(emLogFile != NULL) {
		fclose(emLogFile);
	}
	freeCheck(WattersonsArray);
	freeCheck(siteBoundaryArray);
	freeCheck(auxStats);
	freeCheck(factorsArray);
	freeCheck(allCounts);

	if(emStatus->endState[0] == 't') {
		return 0;
	} else if(emStatus->endState[0] == 'c') {
		return 1;
	} else if(emStatus->endState[0] == 'o') {
		return 2;
	} 
	return -1;
}
/** end of selectionEM **/



/** beta1EM
    Performs inference of beta1/(beta1+beta3) using data at 'L' sites (low MAF) and
    pre-estimated lambda parameters (thetas are not needed)
    @param data Data summary for 'L' sites to consider
    @param lambda Pre-estimated lambda parameters
    @param beta1_p  At entry contains initial value for beta1/(beta1+beta3)
                    At exit contains MLE for beta1/(beta1+beta3)
    		lambda parameters should be set here. Also initial assignment for beta1Thetas and selection parameters are not used.
    @param emSetup struct with various EM setup arguments 
    @param emStatus used to return EM status (likelihood and other information)
	@return  0 if didn't reach convergence
	@return  1 if successfully converged
	@return  2 if converged but overshot (last update reduced likelihood)
	@return -1 if found error
*/
int beta1EM(CollectionData* data, double* lambda, double* beta1_p, EMsetup* emSetup, EMstatus* emStatus) {
	
	/*------------------------------------------------------------------------*/
	/*----------------------------   VARIABLES    ----------------------------*/
	/*------------------------------------------------------------------------*/
	
	FILE*  emLogFile        = NULL;
	char*  emLogFileName    = emSetup->logFileName;
	int    maxNumIterations = emSetup->maxNumIterations;
	int    logInterval      = emSetup->logInterval;
	double lnLd_epsilon     = emSetup->lnLd_epsilon;
	
	/*--- data ---*/
	int numSites = data->numSites;
	int numBlocks = data->numBlocks;

	/*--- set initial value ---*/
	double	beta1 = *beta1_p;
	
	/*--- expected counts associated with low DAF sites (among low MAF sites) ---*/
	double beta1Counts;
	double qBeta1, qBeta3, qSum;
	
	/*--- number of sites in the dataset in each polymorphism class and block---*/
	int* numSitesPerBlock	= data->numSitesPerBlock;	// total number of sites per block
	
	
	/*--- auxilliary stats maintained for each site ---*/
	/*--- these values are computed before the EM loop and do not change throughout the algorithm ---*/
	double pZeqXmaj, pZeqXmin;
	AuxStat_beta1EM* auxStats = NULL;	

	/*--- indexing variables ---*/
	int block, site, siteInBlock;
	
	
	/*--- used in EM ---*/
	int emIteration;
	double lnLd, lnLd_prev, lnLd_diff;
	double expectedLnLd, expectedLnLd_prev, expectedLnLd_diff;
	
	unsigned short zeroLikelihood;		// flag indicating that site with zero likelihood have been encountered
	
	
	/*-----------------------------------------------------------------*/
	/*------------------------   ALLOCATION   -------------------------*/
	/*-----------------------------------------------------------------*/
	
	auxStats        = (AuxStat_beta1EM*)malloc((numSites) * sizeof(AuxStat_beta1EM));
	
	if(auxStats == NULL) {
		fprintf(stdout, "*** Error: OOM allocating space for arrays in beta1EM.\n");
		freeCheck(auxStats);
		return -1;
	}

	/*-----------------------------------------------------------------------------*/
	/*------------------------   INITIALIZATION - START  --------------------------*/
	/*-----------------------------------------------------------------------------*/
	
	/*-------------------------------------*/
	/*---   INIT AUXILLIARY EXPRESSIONS ---*/
	/*-------------------------------------*/
	site=0;
	for(block=0; block<numBlocks; block++) {
		for(siteInBlock=0; siteInBlock<numSitesPerBlock[block]; siteInBlock++) {
         if(data->siteData[block][siteInBlock].polyClass == POLY_L) {
			   pZeqXmaj = data->siteData[block][siteInBlock].pZeqXmaj;
			   pZeqXmin = data->siteData[block][siteInBlock].pZeqXmin;
			
			   auxStats[site].beta1 = pZeqXmaj*(1-lambda[block]) + (1-pZeqXmaj)*lambda[block]/3;
			   auxStats[site].beta3 = pZeqXmin*(1-lambda[block]) + (1-pZeqXmin)*lambda[block]/3;
            site++;
         }
		}// end of for(site)
	}// end of for(block)
	
   numSites=site;

	/*----------------------------------------------------------------------------*/
	/*----------------------------   INIT LOG FILE    ----------------------------*/
	/*----------------------------------------------------------------------------*/
	if(emLogFileName[0] != '\0') {
		emLogFile = fopen(emLogFileName,"w");
		if(emLogFile == NULL) {
			fprintf(stdout, "*** Error: could not open log file '%s' in selectionEM.\n", emLogFileName);
			freeCheck(auxStats);
			return -1;
		}
		fprintf(emLogFile, 	"%7s %8s %10s %10s %10s %10s\n", 
							"iter", "beta1", "lnLd", "lnLd_df", "E[lnLd]", "E[lnLd]_df");
		fflush(emLogFile);
	}

	
	/*---------------------------------------------------------------------------*/
	/*------------------------   INITIALIZATION - END  --------------------------*/
	/*---------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------*/
	/*----------------------------   EM LOOP    ----------------------------*/
	/*----------------------------------------------------------------------*/
	
	emIteration = 0;
	zeroLikelihood = 0;
	lnLd = lnLd_prev = lnLd_diff = 0.0;
	expectedLnLd = expectedLnLd_prev = 0.0;
	strcpy(emStatus->endState,"running");
	printf("Progress: ");
	fflush(stdout);
	while(1) {
		
		lnLd_prev = lnLd;
		lnLd = 0.0;

		/*-------------------------------------*/
		/*---------  E STEP - START  ----------*/
		/*-------------------------------------*/

		beta1Counts= 0.0;
		
		/*--- loop through all sites ---*/
		// for(block=0; block<numBlocks; block++) {
			
		for(site=0; site<numSites; site++) {
			/*--- likelihood components ---*/
			qBeta1 =   beta1   * auxStats[site].beta1;
			qBeta3 = (1-beta1) * auxStats[site].beta3;
			
			/*--- contribution of site pattern to likelihood ---*/
			qSum   = qBeta1 + qBeta3;
			if(qSum == 0.0) {
				lnLd = -DBL_MAX;
				strcpy(emStatus->endState,"zero-likelihood");
				zeroLikelihood = 1;
				break;
			}
			lnLd += customLog(qSum);

			/*--- contribution of site pattern to counts ---*/
			beta1Counts += qBeta1 / qSum;
		}// end of for(site)
			
		if(zeroLikelihood) {
			break;
		}
		// }// end of for(block) - E step and current log-likelihood computation
	
		/*-----------------------------------*/
		/*---------  E STEP - END  ----------*/
		/*-----------------------------------*/

		/*---------------------------------------------*/
		/*--------------   LOG PROGRESS  --------------*/
		/*---------------------------------------------*/

		lnLd_diff = lnLd - lnLd_prev;
		expectedLnLd_diff = expectedLnLd - expectedLnLd_prev;
		emIteration++;
		if (emIteration % logInterval == 0) {
			printf(".");
			if(emLogFile != NULL) {
				fprintf(emLogFile, "%7d %7lf %10lf %10lf %10lf %10lf\n",
						emIteration, beta1, lnLd, lnLd_diff, expectedLnLd, expectedLnLd_diff);
				fflush(emLogFile);
			}
			if(emIteration%(logInterval*10) == 0) {
				printf(" ");
				if(emIteration%(logInterval*100) == 0) {
					printf("\n");
					printf("          ");
				}
			}
			fflush(stdout);
		}

		/*--------------------------------------------------------*/
		/*--------------   UPDATE STATUS AND BREAK  --------------*/
		/*--------------------------------------------------------*/

		if(zeroLikelihood) {
			break;
		}
		
		if(emIteration > 1) {
			if(lnLd_diff < lnLd_epsilon) {
				if(lnLd_diff < 0.0) {
					strcpy(emStatus->endState,"overshoot");
				} else {
					strcpy(emStatus->endState,"converged");
				}
				break;
			} else if(emIteration >= maxNumIterations){
				strcpy(emStatus->endState,"timeout");
				break;
			}
		}// end of if(emIteration > 0)


		/*-------------------------------------*/
		/*---------  M STEP - START  ----------*/
		/*-------------------------------------*/

		/*-------------------------------------------------------------------*/
		/*--------------   EXPECTED LOG LIKELIHOOD PRE UPDATE  --------------*/
		/*-------------------------------------------------------------------*/
		
		expectedLnLd_prev   = 0.0;
		if(beta1 > 0.0) {
			expectedLnLd_prev += beta1Counts*customLog(beta1);
		} else if(beta1Counts > 0) {
			expectedLnLd_prev = -DBL_MAX;
		}
		if(beta1 < 1.0) {
			expectedLnLd_prev += (numSites-beta1Counts)*customLog(1-beta1);
		} else if(beta1Counts < numSites) {
			expectedLnLd_prev = -DBL_MAX;
		}

		/*--------------------------------------------------*/
		/*--------------   UPDATE PARAMETER   --------------*/
		/*--------------------------------------------------*/
		
		beta1 = beta1Counts/numSites;

		/*--------------------------------------------------------------------*/
		/*--------------   EXPECTED LOG LIKELIHOOD POST UPDATE  --------------*/
		/*--------------------------------------------------------------------*/
		
		if(beta1 > 0.0) {
			expectedLnLd += beta1Counts*customLog(beta1);
		} else if(beta1Counts > 0) {
			printf("\n\nError: beta1Counts=%g with beta1 updated to %g.\n",beta1Counts,beta1);
			expectedLnLd_prev = -DBL_MAX;
		}
		if(beta1 < 1.0) {
			expectedLnLd_prev += (numSites-beta1Counts)*customLog(1-beta1);
		} else if(beta1Counts < numSites) {
			printf("\n\nError: numSites-beta1Counts=%g with beta1 updated to %g.\n",numSites-beta1Counts,beta1);
			expectedLnLd_prev = -DBL_MAX;
		}

		/*-----------------------------------*/
		/*---------  M STEP - END  ----------*/
		/*-----------------------------------*/
		
	} // end of while(1)  -- EM update loop
	
	
	/*-----------------------------------------------------------------------*/
	/*----------------------------   FINALIZE    ----------------------------*/
	/*-----------------------------------------------------------------------*/
	
	printf("\n");
	

	/*--- write inferred beta1 into pointer space ---*/
	*beta1_p = beta1;
	
	emStatus->numIterations = emIteration;
	emStatus->lnLikelihood  = lnLd;
	emStatus->lnLd_diff     = lnLd_diff;

	if(emLogFile != NULL) {
		fclose(emLogFile);
	}
	freeCheck(auxStats);

	if(emStatus->endState[0] == 't') {
		return 0;
	} else if(emStatus->endState[0] == 'c') {
		return 1;
	} else if(emStatus->endState[0] == 'o') {
		return 2;
	} 
	return -1;
}
/** end of beta1EM **/



/** computeFisherConfidence
    computes confidence intervals for model parameters based on 
    Fisher's Information Matrix.
    @param elementData data summary
    @param params estimated MLE model parameters
    @param emSetup EMsetup structure to see which parameters are fixed
    @param confidence_out empty 5-point array for output (order: rho, eta, gamma, alpha, tau)
    @return 0 if all OK, -1 otherwise
*/
int computeFisherConfidence(CollectionData* elementData, ModelParameters* params, EMsetup* emSetup, double* confidence_out){
	
	double  fisherIM_space[9], varianceMatrix_space[9];
	double* fisherIM[3];
	double* varianceMatrix[3];
	
	double alpha, tau, Dp, Pw, a, b, c, val;	// ALPHA-TAU
	
	int i, j;
	
	fisherIM[0] = fisherIM_space;
	fisherIM[1] = fisherIM_space+3;
	fisherIM[2] = fisherIM_space+6;
	varianceMatrix[0] = varianceMatrix_space;
	varianceMatrix[1] = varianceMatrix_space+3;
	varianceMatrix[2] = varianceMatrix_space+6;
	
	confidence_out[0] = confidence_out[1] = confidence_out[2] = confidence_out[3] = confidence_out[4] = -1.0;
	
	computeFisherIM(elementData, params, emSetup->useSimplifiedModel, fisherIM);
	
//	printf("Fisher's Information matrix:\n");
//	print3X3matrix(stdout, fisherIM);
	
	/*--- adjust for fixed parameters ---*/
	
	if(emSetup->fixRho) {
		fisherIM[0][0] = 1.0;
		fisherIM[0][1] = fisherIM[0][2] = fisherIM[1][0] = fisherIM[2][0] = 0.0;
	}
	if(emSetup->fixEta) {
		fisherIM[1][1] = 1.0;
		fisherIM[1][0] = fisherIM[1][2] = fisherIM[0][1] = fisherIM[2][1] = 0.0;
	}
	if(emSetup->fixGamma) {
		fisherIM[2][2] = 1.0;
		fisherIM[2][0] = fisherIM[2][1] = fisherIM[0][2] = fisherIM[1][2] = 0.0;
	}

	invertFisherIM(fisherIM, varianceMatrix);

//	printf("Variance-covariance matrix:\n");
//	print3X3matrix(stdout, varianceMatrix);

	/*--- adjust variance = 0.0 for fixed parameters ---*/
	if(emSetup->fixRho) {
		varianceMatrix[0][0] = 0.0;
	}
	if(emSetup->fixEta) {
		varianceMatrix[1][1] = 0.0;
	}
	if(emSetup->fixGamma) {
		varianceMatrix[2][2] = 0.0;
	}
	
	/*--- check that variance matrix is symmetric and all positive ---*/
	for(i=0; i<3; i++) {
		if(varianceMatrix[i][i] < 0.0) {
			fprintf(stdout, "*** Error:  when computing FIM-based confidence intervals:\n");
			fprintf(stdout, "     - entry [%d,%d] of variance-covariance matrix is negative.\n",i+1,i+1);
			print3X3matrix(stdout, varianceMatrix);
			return -1;
		}
		for(j=0; j<i; j++) {
			if(fabs(varianceMatrix[i][j]/varianceMatrix[j][i]-1) > 1e-7) {
				fprintf(stdout, "*** Error:  when computing FIM-based confidence intervals:\n");
				fprintf(stdout, "     - variance-covariance matrix found to be non-symmetric at entry [%d,%d] (diff=%g).\n",
								i+1,j+1, varianceMatrix[i][j] -  varianceMatrix[j][i]);
				print3X3matrix(stdout, varianceMatrix);
				return -1;
			}
		}// end of for(j)
	}// end of for(i)
	
	confidence_out[0] = sqrt(varianceMatrix[0][0]);
	confidence_out[1] = sqrt(varianceMatrix[1][1]);
	confidence_out[2] = sqrt(varianceMatrix[2][2]);
	
	/*--- compute confidence for derived parameters Dp and Pw ---*/
	// FIXME : fix for rho = 0, eta = 0, or gamma = 0
	Dp = params->rho * params->eta * params->meanLambda;
	a = 1 /params->rho;
	b = 1 / params->eta;

	val = a*a*varianceMatrix[0][0] + b*b*varianceMatrix[1][1] + 2*a*b*varianceMatrix[0][1];
	
	if(val < 0.0) {
		fprintf(stdout, "*** Error:  when computing FIM-based confidence intervals:\n");
		fprintf(stdout, "     - Var(Dp) estimated to be negative: %g.\n",val*Dp*Dp*Dp*Dp);
		print3X3matrix(stdout, varianceMatrix);
	} else {
		confidence_out[3] = Dp * sqrt(val);
	}
	
	Pw = params->rho * params->gamma * (params->meanTheta - params->eta * params->meanLamThe);
	a = 1 /params->rho;
	b = -params->meanLamThe / (params->meanTheta - params->eta * params->meanLamThe);
	c = 1 /params->gamma;
//	printf("CIs for tau: a = %g, b= %g, gamma=%g.\n", a, b, params->gamma);
	
	val =	a*a*varianceMatrix[0][0] + b*b*varianceMatrix[1][1] + c*c*varianceMatrix[2][2] + 
			2*a*b*varianceMatrix[0][1] + 2*a*c*varianceMatrix[0][2] + 2*b*c*varianceMatrix[1][2];
	
	if(val < 0.0) {
		fprintf(stdout, "*** Error:  when computing FIM-based confidence intervals:\n");
		fprintf(stdout, "     - Var(Pw) estimated to be negative: %g.\n",val*Pw*Pw*Pw*Pw);
		print3X3matrix(stdout, varianceMatrix);
	} else {
		confidence_out[4] = Pw * sqrt(val);
	}

/*-------------   ALPHA-TAU --------------------------------------------------------------------- */
	if(params->rho == 0.0 || params->eta == 0.0) {
		alpha = 0.0;
		confidence_out[5] = -1;
	} else {
		alpha = 1/(1+(1-params->rho)/(params->rho*params->eta));
		a = 1 / (params->rho * params->rho * params->eta);
		b = (1-params->rho) / (params->rho * params->eta * params->eta);

		val = a*a*varianceMatrix[0][0] + b*b*varianceMatrix[1][1] + 2*a*b*varianceMatrix[0][1];
	
		if(val < 0.0) {
			fprintf(stdout, "Error when computing FIM-based confidence intervals:\n");
			fprintf(stdout, " - Var(alpha) estimated to be negative: %g.\n",val*alpha*alpha*alpha*alpha);
			print3X3matrix(stdout, varianceMatrix);
		} else {
			confidence_out[5] = alpha * alpha * sqrt(val);
		}
	}
	
	if(params->rho == 0.0 || params->gamma == 0.0) {
		tau = 0.0;
		confidence_out[6] = -1;
	} else {
		tau   = 1/(1+(1-params->rho)/(params->rho*params->gamma));
		a = 1 / (params->rho * params->rho * params->gamma);
		b = (1-params->rho) / (params->rho * params->gamma * params->gamma);
//		printf("CIs for tau: a = %g, b= %g, gamma=%g.\n", a, b, params->gamma);
	
		val = a*a*varianceMatrix[0][0] + b*b*varianceMatrix[2][2] + 2*a*b*varianceMatrix[0][2];
	
		if(val < 0.0) {
			fprintf(stdout, "Error when computing FIM-based confidence intervals:\n");
			fprintf(stdout, " - Var(tau) estimated to be negative: %g.\n",val*tau*tau*tau*tau);
			print3X3matrix(stdout, varianceMatrix);
		} else {
			confidence_out[6] = tau * tau * sqrt(val);
		}
	}
/*-------------   ALPHA-TAU ---------------------------------------------------------------------*/

	return 0;
}
/** end of computeFisherConfidence **/



/** computeFisherIM
    computes Fisher's Information Matrix (rows indexed as rho, eta, gamma)
    @param elementData data summary
    @param params parameter estimates (estimated by EM)
    @param useSimplifiedModel boolean flag indicating whether or not to use simplified version of the model
    @param fisherIM empty 3X3 matrix for output
    @return 0 if all OK, -1 otherwise
    FIXME - will need to adapt to new model
*/
int computeFisherIM(CollectionData* elementData, ModelParameters* params, unsigned short useSimplifiedModel, double** fisherIM) {
	
	double* lambda = params->lambda;
	double* theta  = params->theta;
	double* beta   = params->beta;
	double	rho    = params->rho;
	double	eta    = params->eta;
	double	gamma  = params->gamma;
	
	int numBlocks = elementData->numBlocks;
	
	int block;
	int site;
	
	double thetaTimesA;
	double pZeqXmaj, pZeqXmin;
	
	/*--- core values recorded for each site to determine its contribution to the FIM ---*/
	double valA, valB, valC, valD, valE, valF, valP;
	/*--- auxilliary expressions ---*/
	double derivRho, derivEta, derivGamma;
	
	double* WattersonsArray = NULL;

		
	/*-----------------------------------------------------------------*/
	/*------------------------   ALLOCATION   -------------------------*/
	/*-----------------------------------------------------------------*/
	
	WattersonsArray = (double*)malloc((elementData->maxNumSamples+1) * sizeof(double));
	
	if(WattersonsArray == NULL) {
		fprintf(stdout, "*** Error: OOM allocating space for array in computeFisherIM.\n");
		return -1;
	}
	
		
	/*-----------------------------------------------------------------*/
	/*---------------------   INITIALIZE MATRIX   ---------------------*/
	/*-----------------------------------------------------------------*/
	
	/*--- pre-compute Watterson's a values for all possible numbers of samples ---*/
	computeWattersons_a(elementData->maxNumSamples, WattersonsArray);

	fisherIM[0][0] = fisherIM[0][1] = fisherIM[0][2] = 0.0;
	fisherIM[1][0] = fisherIM[1][1] = fisherIM[1][2] = 0.0;
	fisherIM[2][0] = fisherIM[2][1] = fisherIM[2][2] = 0.0;
		
	/*-----------------------------------------------------------------*/
	/*-------------------   MAIN LOOP ON ALL DATA   -------------------*/
	/*-----------------------------------------------------------------*/

	for(block=0; block<numBlocks; block++) {
		for(site=0; site<elementData->numSitesPerBlock[block]; site++) {
			thetaTimesA = theta[block]*WattersonsArray[ elementData->siteData[block][site].numSamples ];
			pZeqXmaj    = elementData->siteData[block][site].pZeqXmaj;
			pZeqXmin    = elementData->siteData[block][site].pZeqXmin;
			/*--- monomorphic sites ---*/
			if(MONOMORPHIC == elementData->siteData[block][site].polyClass) {
				if(useSimplifiedModel) {
					/*--- non-divergent ---*/
					if(pZeqXmaj == 1.0) {
						valA = 1.0;
						valB = 1.0;
						valC = -lambda[block];
						valD = -thetaTimesA;
						valE = thetaTimesA*lambda[block];
						valF = (1-lambda[block]) * (1-thetaTimesA);
					} 
					/*--- divergent ---*/
					else {
						valA = 1.0 / 3.0;
						valB = 0.0;
						valC = lambda[block];
						valD = 0.0;
						valE = 0.0;
						valF = lambda[block] * (1-thetaTimesA) / 3.0;
					}
				} else {
					valA = 1.0;
					valB = pZeqXmaj;
					valC = lambda[block]*((1-pZeqXmaj)/3.0 - pZeqXmaj);
					valD = -pZeqXmaj*thetaTimesA;
					valE = pZeqXmaj*thetaTimesA*lambda[block];
					valF = (1-thetaTimesA) * ( (1-lambda[block]) * pZeqXmaj + lambda[block]/3 * (1-pZeqXmaj) );
				}
			}
			/*--- polymorphic low ---*/
			else if(POLY_L == elementData->siteData[block][site].polyClass) {
				valA = thetaTimesA / 3.0;
				valB = 0.0;
				valC = 0.0;
				valD = 1.0;
				valE = -lambda[block];
				if(useSimplifiedModel) {
					valF = (1-lambda[block]) * beta[0] * thetaTimesA / 3.0;
				} else {
					valA *= pZeqXmaj;
					valF = 
						 lambda[block]/3  * (beta[0] * (1-pZeqXmaj) + beta[2] * (1-pZeqXmin) ) * thetaTimesA / 3 +
						(1-lambda[block]) * (beta[0] *   pZeqXmaj   + beta[2] *   pZeqXmin   ) * thetaTimesA / 3;
				}
			}
			/*--- polymorphic high ---*/
			else if(POLY_H == elementData->siteData[block][site].polyClass) {
				valA = 0.0;
				valB = 0.0;
				valC = 0.0;
				valD = 0.0;
				valE = 0.0;

				if(useSimplifiedModel) {
					valF = (1.0-lambda[block]) * beta[1]*thetaTimesA / 3.0;
				} else {
					valF = ((1.0-2.0*lambda[block]/3.0) * (pZeqXmaj+pZeqXmin) + (2.0*lambda[block]/3.0) * (1-pZeqXmaj-pZeqXmin) )
							* beta[1]*thetaTimesA/3.0;
				}
			}
			/*--- error ---*/
			else  {
				fprintf(stdout, "Error: unknown polymorphism class %d for site %d in block %d.\n", elementData->siteData[block][site].polyClass, site+1, block+1);
				freeCheck(WattersonsArray);
				return -1;
			}
			
			/*--- first parital derivatives of P[data at site] ---*/
			derivRho     = valA * (valB + valC*eta + valD*gamma + valE*eta*gamma) - valF;
			derivEta     = valA * (valC + valE*gamma);
			derivGamma   = valA * (valD + valE*eta);

			/*--- P[data at site] ---*/
			valP = derivRho*rho + valF;
			
			if(valP == 0.0) {
				fprintf(stdout, "Error: zero valP term for site %d of type %d in block %d.\n", site+1, elementData->siteData[block][site].polyClass, block+1);
				//fprintf(stdout, " P[Z=Maj] = %lf,  P[Z=Min] = %lf\n", pZeqXmaj, pZeqXmin);
				fprintf(stdout, " valA=%g, valB=%g, valC=%g, valD=%g, valE=%g, valF=%g, theta=%g, thetaTimesA = %g\n", 
									valA, valB, valC, valD, valE, valF,theta[block], thetaTimesA);
				freeCheck(WattersonsArray);
				return -1;
			}
			
			fisherIM[0][0] += (derivRho * derivRho) / (valP*valP);
			fisherIM[1][1] += (derivEta * derivEta * rho * rho) / (valP*valP);
			fisherIM[2][2] += (derivGamma * derivGamma * rho * rho) / (valP*valP);

			fisherIM[0][1] -= (derivEta * valF) / (valP*valP);
			fisherIM[0][2] -= (derivGamma * valF) / (valP*valP);
			fisherIM[1][2] -= (valA * valE * valF * rho * (1-rho)) / (valP*valP);
			fisherIM[1][2] -= (valA * valA * rho * rho * (valB*valE - valC*valD)) / (valP*valP);

		}// end of for(site)
	}// end of for(block)
	
	/*--- copy part above diagonal below diagonal ---*/
	
	fisherIM[1][0] = fisherIM[0][1];
	fisherIM[2][0] = fisherIM[0][2];
	fisherIM[2][1] = fisherIM[1][2];
			
	freeCheck(WattersonsArray);
	return 0;
}
/** end of computeFisherIM **/



/** computePosteriors
    computes posterior probabilities and counts of hidden 
    @param elementData data summary
    @param params parameter estimates (estimated by EM)
    @param useSimplifiedModel boolean flag indicating whether or not to use simplified version of the model
    @param outFileName name of file to which to write output
    @param totalCounts pre-allocated array (with 11 entries) for total counts
    @return 0 if all OK, -1 otherwise
*/
int computePosteriors(CollectionData* elementData, ModelParameters* params, unsigned short useSimplifiedModel, const char* outFileName, double* totalCounts) {
	
	double* lambda = params->lambda;
	double* theta  = params->theta;
	double* beta   = params->beta;
	double	rho    = params->rho;
	double	eta    = params->eta;
	double	gamma  = params->gamma;
	
	int numBlocks = elementData->numBlocks;
	
	int block;
	int site, i;
	
	double thetaTimesA;
	double pZeqXmaj, pZeqXmin;

	/*--- posterior probabilities and counts for M sites ---*/
	double pNeutDiv, pNeutNodiv, pSelDiv, pSelNodiv;
	double cNeutDiv, cNeutNodiv, cSelDiv, cSelNodiv;

	/*--- posterior probabilities for L sites ---*/
	double pSelLow, pNeutLow, pNeutHigh;
	double cSelLow, cNeutLow, cNeutHigh, cHigh;

	double pTotal;

	FILE* outFile;
	
	double* WattersonsArray = NULL;

		
	/*-----------------------------------------------------------------*/
	/*------------------------   ALLOCATION   -------------------------*/
	/*-----------------------------------------------------------------*/
	
	WattersonsArray = (double*)malloc((elementData->maxNumSamples+1) * sizeof(double));
	
	if(WattersonsArray == NULL) {
		fprintf(stdout, "OOM allocating space for array in computePosteriors.\n");
		return -1;
	}
	
	/*--- open output file and print header line ---*/	
	outFile = fopen(outFileName,"w");
	if(outFile == NULL) {
		fprintf(stdout, "Error: could not open output file '%s' in computePosteriors().\n", outFileName);
		freeCheck(WattersonsArray);
	}
	fprintf(outFile, "site\tNoDiv\tNoDiv_N\tNoDiv_S\tDiv\tDiv_N\tDiv_S\tPoly\tPolyL_N\tPolyH1_N\tPolyH2_N\tPolyL_S\n");

	/*--- pre-compute Watterson's a values for all possible numbers of samples ---*/
	computeWattersons_a(elementData->maxNumSamples, WattersonsArray);

	/*-----------------------------------------------------------------*/
	/*-------------------   MAIN LOOP ON ALL DATA   -------------------*/
	/*-----------------------------------------------------------------*/

	/*--- init counts ---*/
	cNeutDiv = cNeutNodiv = cSelDiv = cSelNodiv = 0.0;
	cSelLow = cNeutLow = cNeutHigh = cHigh = 0.0;
	for(i=0; i<11; i++) {
		totalCounts[i] = 0.0;
	}

	for(block=0; block<numBlocks; block++) {

		for(site=0; site<elementData->numSitesPerBlock[block]; site++) {
			thetaTimesA = theta[block]*WattersonsArray[ elementData->siteData[block][site].numSamples ];
			pZeqXmaj    = elementData->siteData[block][site].pZeqXmaj;
			pZeqXmin    = elementData->siteData[block][site].pZeqXmin;
			/*--- monomorphic sites ---*/
			if(MONOMORPHIC == elementData->siteData[block][site].polyClass) {
				if(useSimplifiedModel) {
					pNeutNodiv = (1-rho)*(1-thetaTimesA)*pZeqXmaj*(1-lambda[block]);
					pNeutDiv   = (1-rho)*(1-thetaTimesA)*(1-pZeqXmaj)*lambda[block];
					pSelNodiv  = rho*(1-gamma*thetaTimesA)*pZeqXmaj*(1-eta*lambda[block]);
					pSelDiv    = rho*(1-pZeqXmaj)*eta*lambda[block];
				} else {
					pNeutNodiv = (1-rho)*(1-thetaTimesA)*pZeqXmaj*(1-lambda[block]);
					pNeutDiv   = (1-rho)*(1-thetaTimesA)*(1-pZeqXmaj)*lambda[block]/3;
					pSelNodiv  = rho*(1-gamma*thetaTimesA)*pZeqXmaj*(1-eta*lambda[block]);
					pSelDiv    = rho*(1-pZeqXmaj)*eta*lambda[block]/3;
				}
				pTotal = pNeutNodiv + pNeutDiv + pSelNodiv + pSelDiv;
				cNeutNodiv += pNeutNodiv/pTotal;
				cNeutDiv   += pNeutDiv/pTotal;
				cSelNodiv  += pSelNodiv/pTotal;
				cSelDiv    += pSelDiv/pTotal;
			}
			/*--- polymorphic low ---*/
			else if(POLY_L == elementData->siteData[block][site].polyClass) {
				if(useSimplifiedModel) {
					pSelLow   = rho * gamma*thetaTimesA * (1 - eta*lambda[block]);
					pNeutLow  = (1-rho) * beta[0] * thetaTimesA * (1-lambda[block]);
					pNeutHigh = 0;
				} else {
					pSelLow   = rho * gamma*thetaTimesA/3 * pZeqXmaj * (1 - eta*lambda[block]);
					pNeutLow  = (1-rho) * beta[0] * thetaTimesA/3 * ((1-pZeqXmaj)*lambda[block]/3 + pZeqXmaj*(1-lambda[block]));
					pNeutHigh = (1-rho) * beta[2] * thetaTimesA/3 * ((1-pZeqXmin)*lambda[block]/3 + pZeqXmin*(1-lambda[block]));
				}
				pTotal = pSelLow + pNeutLow + pNeutHigh;
				cSelLow   += pSelLow/pTotal;
				cNeutLow  += pNeutLow/pTotal;
				cNeutHigh += pNeutHigh/pTotal;
			}
			/*--- polymorphic high ---*/
			else if(POLY_H == elementData->siteData[block][site].polyClass) {
				cHigh += 1.0;
			}
			/*--- error ---*/
			else  {
				fprintf(stdout, "*** Error: unknown polymorphism class %d for site %d in block %d.\n", elementData->siteData[block][site].polyClass, site+1, block+1);
				freeCheck(WattersonsArray);
				return -1;
			}
			if(!useSimplifiedModel) {
				fprintf(outFile, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
						elementData->siteData[block][site].siteId,
						cNeutNodiv+cSelNodiv, cNeutNodiv, cSelNodiv, cNeutDiv+cSelDiv, cNeutDiv, cSelDiv,
						cNeutLow+cNeutHigh+cHigh+cSelLow, cNeutLow,cHigh,cNeutHigh,cSelLow);

				totalCounts[0] += cNeutNodiv+cSelNodiv;
				totalCounts[1] += cNeutNodiv;
				totalCounts[2] += cSelNodiv;
				totalCounts[3] += cNeutDiv+cSelDiv;
				totalCounts[4] += cNeutDiv;
				totalCounts[5] += cSelDiv;
				totalCounts[6] += cNeutLow+cNeutHigh+cHigh+cSelLow;
				totalCounts[7] += cNeutLow;
				totalCounts[8] += cNeutHigh;
				totalCounts[9] += cHigh;
				totalCounts[10] += cSelLow;
				cNeutDiv = cNeutNodiv = cSelDiv = cSelNodiv = 0.0;
				cSelLow = cNeutLow = cNeutHigh = cHigh = 0.0;
			}
		}// end of for(site)
		if(useSimplifiedModel) {
			fprintf(outFile, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
					elementData->blockNames[block],
					cNeutNodiv+cSelNodiv, cNeutNodiv, cSelNodiv, cNeutDiv+cSelDiv, cNeutDiv, cSelDiv,
					cNeutLow+cNeutHigh+cHigh+cSelLow, cNeutLow,cHigh,cNeutHigh,cSelLow);

			totalCounts[0] += cNeutNodiv+cSelNodiv;
			totalCounts[1] += cNeutNodiv;
			totalCounts[2] += cSelNodiv;
			totalCounts[3] += cNeutDiv+cSelDiv;
			totalCounts[4] += cNeutDiv;
			totalCounts[5] += cSelDiv;
			totalCounts[6] += cNeutLow+cNeutHigh+cHigh+cSelLow;
			totalCounts[7] += cNeutLow;
			totalCounts[8] += cNeutHigh;
			totalCounts[9] += cHigh;
			totalCounts[10] += cSelLow;
			cNeutDiv = cNeutNodiv = cSelDiv = cSelNodiv = 0.0;
			cSelLow = cNeutLow = cNeutHigh = cHigh = 0.0;
		}
	}// end of for(block)
	
			
	freeCheck(WattersonsArray);
	fclose(outFile);
	return 0;
}
/** end of computePosteriors **/



/** computeFisherIM_bak - backup version for simplified model
    computes Fisher's Information Matrix (rows indexed as rho, eta, gamma)
    @param elementData data summary
    @param params parameter estimates (estimated by EM)
    @param fisherIM empty 3X3 matrix for output
    @return 0 if all OK, -1 otherwise
    FIXME - will need to adapt to new model
*/
int computeFisherIM_bak(CollectionData* elementData, ModelParameters* params, double** fisherIM) {
	
	double* lambda = params->lambda;
	double* theta  = params->theta;
	double* beta   = params->beta;
	double	rho    = params->rho;
	double	eta    = params->eta;
	double	gamma  = params->gamma;
	
	int numBlocks = elementData->numBlocks;
	
	int block;
	int site;
	
	double thetaTimesA;
	
	/*--- core values recorded for each site to determine its contribution to the FIM ---*/
	double valA, valB, valC, valD, valE, valF, valP;
	/*--- auxilliary expressions ---*/
	double gammaTerm, etaTerm, rhoTerm, rhoEtaTerm, rhoGammaTerm;
	
	double* WattersonsArray = NULL;

		
	/*-----------------------------------------------------------------*/
	/*------------------------   ALLOCATION   -------------------------*/
	/*-----------------------------------------------------------------*/
	
	WattersonsArray = (double*)malloc((elementData->maxNumSamples+1) * sizeof(double));
	
	if(WattersonsArray == NULL) {
		fprintf(stdout, "*** Error: OOM allocating space for arrays in computeFisherIM.\n");
		freeCheck(WattersonsArray);
		return -1;
	}
	
		
	/*-----------------------------------------------------------------*/
	/*---------------------   INITIALIZE MATRIX   ---------------------*/
	/*-----------------------------------------------------------------*/
	
	/*--- pre-compute Watterson's a values for all possible numbers of samples ---*/
	computeWattersons_a(elementData->maxNumSamples, WattersonsArray);

	fisherIM[0][0] = fisherIM[0][1] = fisherIM[0][2] = 0.0;
	fisherIM[1][0] = fisherIM[1][1] = fisherIM[1][2] = 0.0;
	fisherIM[2][0] = fisherIM[2][1] = fisherIM[2][2] = 0.0;
		
	/*-----------------------------------------------------------------*/
	/*-------------------   MAIN LOOP ON ALL DATA   -------------------*/
	/*-----------------------------------------------------------------*/

	for(block=0; block<numBlocks; block++) {
		for(site=0; site<elementData->numSitesPerBlock[block]; site++) {
			thetaTimesA = theta[block]*WattersonsArray[ elementData->siteData[block][site].numSamples ];
			/*--- monomorphic sites ---*/
			if(MONOMORPHIC == elementData->siteData[block][site].polyClass) {
				/*--- non-divergent ---*/
				if(elementData->siteData[block][site].pZeqXmaj == 1.0) {
					valA = 1.0;
					valB = 1.0;
					valC = -lambda[block];
					valD = 1.0;
					valE = -thetaTimesA;
					valF = (1-lambda[block]) * (1-thetaTimesA);
				} 
				/*--- divergent ---*/
				else {
					valA = 1.0 / 3.0;
					valB = 0.0;
					valC = lambda[block];
					valD = 1.0;
					valE = 0.0;
					valF = lambda[block] * (1-thetaTimesA) / 3.0;
				}
			}
			/*--- polymorphic low ---*/
			else if(POLY_L == elementData->siteData[block][site].polyClass) {
				valA = 1.0 / 3.0;
				valB = 1.0;
				valC = -lambda[block];
				valD = 0.0;
				valE = thetaTimesA;
				valF = (1-lambda[block]) * beta[0] * thetaTimesA / 3.0;
			}
			/*--- polymorphic high ---*/
			else if(POLY_H == elementData->siteData[block][site].polyClass) {
				valA = 0.0;
				valB = 0.0;
				valC = 0.0;
				valD = 0.0;
				valE = 0.0;
				valF = (1-lambda[block]) * beta[1] * thetaTimesA / 3.0;
			}
			/*--- error ---*/
			else  {
				fprintf(stdout, "*** Error: unknown polymorphism class %d for site %d in block %d.\n", elementData->siteData[block][site].polyClass, site+1, block+1);
				freeCheck(WattersonsArray);
				return -1;
			}
			
			etaTerm      = valB + valC*eta;
			gammaTerm    = valD + valE*gamma;
			rhoTerm      = valA * etaTerm * gammaTerm - valF;
			rhoEtaTerm   = valA * valC * gammaTerm;
			rhoGammaTerm = valA * valE * etaTerm;
			
			valP = rhoTerm*rho + valF;
			
			if(valP == 0.0) {
				fprintf(stdout, "*** Error: zero valP term for site %d of type %d in block %d.\n", site+1, elementData->siteData[block][site].polyClass, block+1);
				fprintf(stdout, "      valA=%g, valB=%g, valC=%g, valD=%g, valE=%g, valF=%g, theta=%g, thetaTimesA = %g\n", 
									valA, valB, valC, valD, valE, valF,theta[block], thetaTimesA);
				freeCheck(WattersonsArray);
				return -1;
			}
			
			fisherIM[0][0] += (rhoTerm * rhoTerm) / (valP*valP);
			fisherIM[1][1] += (rhoEtaTerm * rhoEtaTerm * rho * rho) / (valP*valP);
			fisherIM[2][2] += (rhoGammaTerm * rhoGammaTerm * rho * rho) / (valP*valP);

			fisherIM[0][1] -= (rhoEtaTerm * valF) / (valP*valP);
			fisherIM[0][2] -= (rhoGammaTerm * valF) / (valP*valP);
			fisherIM[1][2] -= (valA * valC * valE * valF * rho * (1-rho)) / (valP*valP);

		}// end of for(site)
	}// end of for(block)
	
	/*--- copy part above diagonal below diagonal ---*/
	
	fisherIM[1][0] = fisherIM[0][1];
	fisherIM[2][0] = fisherIM[0][2];
	fisherIM[2][1] = fisherIM[1][2];
			
	freeCheck(WattersonsArray);
	return 0;
}
/** end of computeFisherIM_bak **/



/******************************************************************************************************/
/******                                   AUXILLIARY FUNCTIONS                                   ******/
/******************************************************************************************************/



/** freeModelParameters
    frees allocated memory in ModelParameters (when initialized from file)
    @param params ModelParameters structure to be freed
    @return 0 always
*/
int freeModelParameters(ModelParameters* params) {
	
	freeCheck(params->theta);
	freeCheck(params->lambda);
		
	return 0;
}
/** end of freeModelParameters **/



/** freeCollectionData
    frees allocated memory in CollectionData (when initialized from file)
    @param data CollectionData structure to be freed
    @return 0 always
*/
int freeCollectionData(CollectionData* data) {
	
	int block;
	if(data->siteData != NULL) {
		for(block=0; block<data->numBlocks; block++) {
			freeCheck(data->siteData[block]);
		}// end of for(block)
	
		freeCheck(data->siteData);
	}

	freeCheck(data->numSitesPerBlock);
	if(data->blockNames != NULL) {
		freeCheck(data->blockNames[0]);
	}
	freeCheck(data->blockNames);
	
	return 0;
}
/** end of freeCollectionData **/



/** computeWattersons_a
    computes a series of values according Watterson's estimator of theta 
    to accommodate for missing data in consideration of polymorphism rate
    @param numValues Number of values to calculate
    @param valueArray Array in which to output calculated values (array length is numValues+1)
    @return 0 always
*/
int computeWattersons_a(int numValues, double* valueArray) {
	
	int k;
	
	if(numValues < 1) {
		return 0;
	}
	
	valueArray[0] = 0.0;
	valueArray[1] = 0.0;
	for(k=1; k<numValues;k++) {
		valueArray[k+1] = valueArray[k] + (1.0 / k);
	}// end of for(k)
	
	return 0;
}
/** end of computeWattersons_a **/



/** print3X3matrix
    prints a 3X3 matrix
    @param matrix3X3 3X3 matrix
*/
void print3X3matrix(FILE* filePointer, double** matrix3X3) {
	
	fprintf(filePointer, "  %4lf %4lf %4lf\n",matrix3X3[0][0],matrix3X3[0][1],matrix3X3[0][2]);
	fprintf(filePointer, "  %4lf %4lf %4lf\n",matrix3X3[1][0],matrix3X3[1][1],matrix3X3[1][2]);
	fprintf(filePointer, "  %4lf %4lf %4lf\n",matrix3X3[2][0],matrix3X3[2][1],matrix3X3[2][2]);

//	fprintf(filePointer, "  %g %g %g\n",matrix3X3[0][0],matrix3X3[0][1],matrix3X3[0][2]);
//	fprintf(filePointer, "  %g %g %g\n",matrix3X3[1][0],matrix3X3[1][1],matrix3X3[1][2]);
//	fprintf(filePointer, "  %g %g %g\n",matrix3X3[2][0],matrix3X3[2][1],matrix3X3[2][2]);
	
	return;
}
/** end of print3X3matrix **/



/** invertFisherIM
    inverts a 3X3 symmetric matrix (representing Fisher's Information Matrix)
    @param fisherIM 3X3 matrix for Fisher's Information Matrix
    @param invertedMatrix empty 3X3 matrix for output
    @return 0 if all OK, -1 otherwise
*/
int invertFisherIM(double** fisherIM, double** invertedMatrix) {
	
	double  matrixSpace[9];
	double* copiedMatrix[3];
	
	int ind1, ind2, ind3, res = 0;
	double factor;
	
	/*--- initialize copiedMatrix to fisherIM and invertedMatrix to unit matrix ---*/
	/*--- also check that FIM is symmetric matrix with positive values in diagonal ---*/
	for(ind1=0; ind1<3; ind1++) {
		copiedMatrix[ind1] = matrixSpace + 3*ind1;
		for(ind2=0; ind2<3; ind2++) {
			copiedMatrix[ind1][ind2] = fisherIM[ind1][ind2];
			if(ind1 == ind2) {
				invertedMatrix[ind1][ind2] = 1.0;
			} else {
				invertedMatrix[ind1][ind2] = 0.0;
			}				
		}// end of for(ind2)
	}// end of for(ind1)
	
	if(res < 0) {
		return res;
	}
	
/*---DEBUG
	printf("init:\n");
	printf(" mat:\n");
	print3X3matrix(copiedMatrix);
	printf(" invert:\n");	
	print3X3matrix(invertedMatrix);
---*/			
	
	/*--- convert copied matrix to unit matrix, and then the unit matrix transforms to the inverted matrix ---*/
	for(ind1=0; ind1<3; ind1++) {
		for(ind2=0; ind2<3; ind2++) {
			/*--- if ind1==ind2, this multiplies row ind1 by a factor to make copiedMatrix[ind1][ind1]==1 ---*/
			/*--- if ind1!=ind2, this adds to row ind2 a factor of row ind1 to make copiedMatrix[ind2][ind1]==0 ---*/
			if(ind1 == ind2) {
				factor = 1.0 - (1.0 / copiedMatrix[ind1][ind1]);
			} else {
				factor = copiedMatrix[ind2][ind1] / copiedMatrix[ind1][ind1];
			}
			for(ind3=0; ind3<3; ind3++) {
				copiedMatrix[ind2][ind3] -= copiedMatrix[ind1][ind3]*factor;
				invertedMatrix[ind2][ind3] -= invertedMatrix[ind1][ind3]*factor;
			}// end of for(ind3)
			if(ind1 == ind2 && fabs(copiedMatrix[ind1][ind1] - 1.0) > 1e-10) {
				fprintf(stdout, "*** Error when inverting Fisher's Information Matrix: in operation on row %d, diagonal element did not turn 1.0 (diff = %g)\n", ind1+1, copiedMatrix[ind1][ind1]-1.0);
				res = -1;
			} else if(ind1 != ind2 && fabs(copiedMatrix[ind2][ind1]) > 1e-10){
				fprintf(stdout, "*** Error when inverting Fisher's Information Matrix: in operation on rows %d,%d element did not nullify (%g)\n", ind1+1, ind2+1, copiedMatrix[ind2][ind1]);
				res = -1;
			}
/*---DEBUG
			printf("[%d,%d]:\n", ind1+1, ind2+1);
			printf(" mat:\n");
			print3X3matrix(copiedMatrix);
			printf(" invert:\n");	
			print3X3matrix(invertedMatrix);
---*/			
		}// end of for(ind2)
	}// end of for(ind1)
	
	return res;
}
/** end of invertFisherIM **/


/******************************************************************************************************/
/******                                         END OF FILE                                      ******/
/******************************************************************************************************/
