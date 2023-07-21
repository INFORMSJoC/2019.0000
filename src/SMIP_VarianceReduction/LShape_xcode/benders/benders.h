/*
 * benders.h
 *
 *  Created on: Feb 17, 2019
 *      Author: Jiajun Xu
 * Institution: University of Southern California
 *
 * Please send your comments or bug report to jiajunx (at) usc (dot) edu
 *
 */

#ifndef BENDERS_H_
#define BENDERS_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

//#undef ALGO_CHECK
//#undef STOCH_CHECK
#undef SAVE_DUALS

typedef struct{
	int		NUM_REPS;			/* Maximum number of replications that can be carried out. */
	long long *RUN_SEED;		/* seed used during optimization */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTER_TYPE;		/* type of master problem */
	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double	EPSILON;			/* Optimality gap */

	int		EVAL_FLAG;
	int		NUM_EVALS;
	long long *EVAL_SEED;
	int		EVAL_MIN_ITER;
	double	EVAL_ERROR;

	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	int		MULTICUT;			/* Set to 1 if multicut is to be solved */
	double	R1;
	double	R2;
	double	R3;

	int		MAX_OBS;			/* Maximum number of iterations before which SAA is invoked */
	int		SAA; 				/* Use SAA when continuous distribution in stoch file (1), or not (0) */

	int 	MULTIPLE_REP;		/* When multiple replications are needed, set this to (1), else (0) */
    double  Bag_num;
    int     reg;                /*when the master is MILP, if we use regularizer within the replication, 1 represent yes, 0 represent 0*/
}configType;

typedef struct {
	int		ck;					/* Iteration when the cut was generated */
	double  alpha;              /* scalar value for the right-hand side */
	vector  beta;               /* coefficients of the master problems's primal variables */
	BOOL	isIncumb;			/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;		/* right-hand side when using QP master, this is useful for quick updates */
	int 	rowNum;				/* row number for master problem in solver */
	string	name;
}oneCut;

typedef struct {
	int    	cnt;                    /* number of cuts */
	oneCut  **vals;					/* values which define the set of cuts */
}cutsType;

typedef struct {
	double	repTime;
	double 	iterTime;
	double 	masterIter;
	double 	subprobIter;
	double 	optTestIter;
	double 	iterAccumTime;
	double 	masterAccumTime;
	double 	subprobAccumTime;
	double 	optTestAccumTime;
}runTime;

typedef struct {
	int		numRV;					/* Number of random variables */
	int 	cnt;					/* Number of observations */
	vector	probs;					/* Probability of observation */
	vector	*vals;					/* Observation values */
} omegaType;

typedef struct {
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved. */

    oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */

	vector      candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */

	vector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	BOOL        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	vector		piM;

    int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType    *cuts;              /* optimality cuts */
	cutsType    *fCuts;             /* feasibility cuts */

	omegaType 	*omega;				/* all realizations observed during the algorithm */

    BOOL        optFlag;
    BOOL		optMode;
    BOOL		spFeasFlag;			/* Indicates whether the subproblem is feasible */
    int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	BOOL		infeasIncumb;		/* indicates if the incumbent solution is infeasbible */
    int         RepeatedTime;

	runTime		*time;				/* Run time structure */
}cellType;

#if defined(SAVE_DUALS)
typedef struct {
	int 	cnt;
	vector 	*vals;
	intvec	iter;
	intvec  obs;
}dualsType;
#endif

typedef struct {
    oneProblem    *sp;                /* compromise problem */
    int         cnt;                /* number of replications */
    intvec         ck;                    /* number of iterations for each replication */
    vector        objLB;                /* replication lower bound */
    vector        objUB;                /* replication upper bound, if batch solution is evaluated */
    double        objComp;            /* optimal value of compromise problem */
    double        quadScalar;            /* average proximal terms */
    vector        *incumbX;            /* batch incumbent solution */
    vector        compromiseX;        /* compromise solution */
    vector        avgX;                /* average solution across batches */
    runTime        *time;                /* Run time structure */
    double         Est;                 /*lower bound estimation*/
    vector         *distinctX;          /*distict solutions*/
    intvec         repeatTime;          /*appearance time of distinctX*/
    vector         distinctEst;         /*lower bound of the distict solutions*/
    int            distNum;              /*number of distict solutions*/
    int            bagging_idx;          /*index in distinctX that is bagging solution*/
    int            max_repeat;           /*maximum of repeat time*/
}batchSummary;

int parseCmdLine(int argc, char *argv[], string probName, string inputDir);
void createOutputDir(string outputDir, string algoName, string probName);
int readConfig();
void freeConfig();

/* algo.c */
int algo (oneProblem *orig, timeType *tim, stocType *stoc, string probName);
int solveBendersCell(stocType *stoc, probType **prob, cellType *cell);
BOOL optimal(cellType *cell);
void writeStatistic(FILE *soln, FILE *incumb, probType **prob, cellType *cell);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, batchSummary **batch, vector *meanSol);
cellType *newCell(stocType *stoc, probType **prob, vector xk);
int cleanCellType(cellType *cell, probType *prob, vector xk);
void freeCellType(cellType *cell);

/* masters.c */
int solveMaster(numType *num, sparseVector *dBar, cellType *cell);
int addCut2Master(cellType *cell, cutsType *cuts, oneCut *cut, int lenX);
int checkImprovement(probType *prob, cellType *cell, int candidCut);
int replaceIncumbent(probType *prob, cellType *cell);
int constructQP(probType *prob, cellType *cell, vector incumbX, double quadScalar);
int changeQPproximal(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell, vector xk);
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk);

int constructMILP(probType *prob, cellType *cell, vector incumbX, double quadScalar);
int changeMILPwithL1(LPptr lp, oneProblem *sp, int numCols);
int changeMILPproximal(LPptr lp, vector objx, int numCols, double sigma);
int changeMILPrhs(probType *prob, cellType *cell, vector xk);
int changeMILPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk, int offset);

oneProblem *newMaster(oneProblem *orig, double lb);

/* cuts.c */
int formOptCut(probType *prob, cellType *cell, vector Xvect, BOOL isIncumb);
double maxCutHeight(cutsType *cuts, vector xk, int betaLen);
double cutHeight(oneCut *cut, vector xk, int betaLen);
int reduceCuts(oneProblem *master, cutsType *cuts, vector vectX, vector piM, int betaLen, int *iCutIdx);
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, int *iCutIdx);
oneCut *newCut(int numX, int currentIter);
cutsType *newCuts(int maxCuts);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts, BOOL partial);

/* subprob.c */
int solveSubprob(probType *prob, oneProblem *subproblem, vector Xvect, vector obsVals, BOOL *spFeasFlag, double *subprobTime, vector piS, double *mubBar);
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs);
vector computeCostCoeff(numType *num, coordType *coord, sparseVector *dBar, vector obs);
int computeMU(LPptr lp, int numCols, double *mubBar);
oneProblem *newSubproblem(oneProblem *subprob);
void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, vector rhs, vector X);
int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, vector observ, vector spRHS, vector X);
int chgObjxwObserv(LPptr lp, numType *num, coordType *coord, vector cost, intvec indices, vector observ);
omegaType *newOmega(stocType *stoc);
void freeOmegaType(omegaType *omega, BOOL partial);

/* evaluate.c */
int evaluate(FILE *soln, stocType *stoc, probType **prob, cellType *cell, vector Xvect, double *ub_mean, double *ub_std);

/* compromise.c */
int buildCompromise(probType *prob, cellType *cell, batchSummary *batch);
int solveCompromise(probType *prob, batchSummary *batch);
int addL1Norm (probType *prob, batchSummary *batch);
int addBatchEquality (probType *prob, batchSummary *batch);
batchSummary *newBatchSummary(probType *prob, int numBatch);
void freeBatchType(batchSummary *batch);

BOOL InConvexHull(batchSummary *batch, int col);
void LowerBoundVariance(batchSummary *batch, double *mean, double *std);

#endif /* BENDERS_H_ */
