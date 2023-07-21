/*
 * setup.c
 *
 *  Created on: Feb 17, 2019
 *      Author: Jiajun Xu
 * Institution: University of Southern California
 *
 * Please send your comments or bug report to jiajunx (at) usc (dot) edu
 *
 */

#include "benders.h"

extern configType config;

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, batchSummary **batch, vector *meanSol) {
	vector	lb = NULL;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	(*meanSol) = meanProblem(orig, stoc, tim->col[1] - tim->col[0]);
	if ( (*meanSol) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		goto TERMINATE;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		goto TERMINATE;
	}

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(stdout, orig->name, tim, (*prob));
#endif

	/* ensure that we have a linear programs at all stages */
	t = 0;
	while ( t < tim->numStages ) {
		if ( (*prob)[t++]->sp->type  != PROB_LP )
			printf("Note :: Stage-%d problem is a mixed-integer program.\n", t);
	}

	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), (*meanSol));
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		goto TERMINATE;
	}
    
    if ( config.NUM_REPS > 1 )
        (*batch)  = newBatchSummary((*prob)[0], config.NUM_REPS);

	mem_free(lb);
	return 0;
	TERMINATE:
	mem_free(lb);
	return 1;
}//END setupAlgo()

/* This function is used to create cells used in the algorithm */
cellType *newCell(stocType *stoc, probType **prob, vector xk) {
	cellType    *cell;

	/* Allocate memory to all cells used in the algorithm. */
	if (!(cell = (cellType *) mem_malloc(sizeof(cellType))) )
		errMsg("Memory allocation", "new_cell", "failed to allocate memory to cell",0);
	cell->master = cell->subprob = NULL;
	cell->cuts = cell->fCuts = NULL;
	cell->omega = NULL; cell->piM = NULL;

	cell->candidX = cell->incumbX = NULL;

	/* Jiajun check, setup the master problem */
	cell->master = newMaster(prob[0]->sp, prob[0]->lb);
	if ( cell->master == NULL ) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}
	/* setup the subproblem */
	cell->subprob = newSubproblem(prob[1]->sp);

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to cell +-+-+-+-+-+-+-+-+-+- */
	cell->k 	= 0;
	cell->LPcnt = 0;
    cell->RepeatedTime = 0;

	/* candidate solution and estimates, Jiajun check */
	cell->candidX 	= duplicVector(xk, prob[0]->num->cols);
	cell->candidEst	= prob[0]->lb + vXvSparse(cell->candidX, prob[0]->dBar);

	/* stochastic elements */
	cell->omega  = newOmega(stoc);

	cell->optFlag 			= FALSE;

	cell->spFeasFlag = TRUE;
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= FALSE;

	/* incumbent solution and estimates */
	if (config.MASTER_TYPE == PROB_QP || config.reg == 1) {
		cell->incumbX   = duplicVector(xk, prob[0]->num->cols);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;     						/* The quadratic scalar, 'sigma'*/
		cell->iCutIdx   = 0;
		cell->incumbChg = TRUE;

		cell->maxCuts = config.CUT_MULT * prob[0]->num->cols + 3;

		if ( !(cell->piM = (vector) arr_alloc(prob[0]->num->rows + cell->maxCuts + 1, double)) )
			errMsg("allocation", "newMaster", "cell->piM", 0);
	}
	else {
		cell->incumbX   = NULL;
		cell->incumbEst = 0.0;
		cell->quadScalar= 0.0;
		cell->iCutIdx   = -1;
		cell->incumbChg = FALSE;

		cell->maxCuts = config.MAX_ITER;
	}
	cell->cuts 	  = newCuts(cell->maxCuts);
	cell->fCuts		= newCuts(cell->maxCuts);

	if ( !(cell->time = (runTime *) mem_malloc(sizeof(runTime)) ) )
		errMsg("setup", "newCell", "cell->runTime", 0);
	cell->time->repTime = cell->time->iterTime = cell->time->masterIter = cell->time->subprobIter = cell->time->optTestIter = 0.0;
	cell->time->iterAccumTime = cell->time->masterAccumTime = cell->time->subprobAccumTime = cell->time->optTestAccumTime = 0.0;

	/* construct the QP using the current incumbent */
	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob[0], cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}

		cell->incumbChg = FALSE;
#if defined(ALGO_CHECK)
		if ( writeProblem(cell->master->lp, "newQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return NULL;
		}
#endif
	}
    if (config.reg == 1){
        if ( constructMILP(prob[0], cell, cell->incumbX, cell->quadScalar) ) {
            errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
            return NULL;
    }
        
        cell->incumbChg = FALSE;
#if defined(ALGO_CHECK)
        if ( writeProblem(cell->master->lp, "newMILPMaster.lp") ) {
            errMsg("write problem", "new_master", "failed to write master problem to file",0);
            return NULL;
        }
#endif
    }

	return cell;
}//END newCell()

int cleanCellType(cellType *cell, probType *prob, vector xk) {
	int cnt;

	/* constants and arrays */
	cell->k = 0;
	cell->LPcnt = 0;
	cell->optFlag 		 = FALSE;
	cell->spFeasFlag 	 = TRUE;
    cell->RepeatedTime = 0;

	copyVector(xk, cell->candidX, prob->num->cols, TRUE);
	cell->candidEst	= prob->lb + vXvSparse(cell->candidX, prob->dBar);

	if (config.MASTER_TYPE == PROB_QP) {
		copyVector(xk, cell->incumbX, prob->num->cols, TRUE);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;
		cell->iCutIdx   = 0;
		cell->incumbChg = TRUE;
	}

	/* oneProblem structures and solver elements */
	for ( cnt = prob->num->rows+cell->cuts->cnt+cell->fCuts->cnt-1; cnt >= prob->num->rows; cnt-- )
		if (  removeRow(cell->master->lp, cnt, cnt) ) {
			errMsg("solver", "cleanCellType", "failed to remove a row from master problem", 0);
			return 1;
		}
	cell->master->mar = prob->num->rows;
	if( changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar)) {
		errMsg("algorithm", "cleanCellType", "failed to change the proximal term", 0);
		return 1;
	}

	/* cuts */
	if (cell->cuts) freeCutsType(cell->cuts, TRUE);
	if (cell->fCuts) freeCutsType(cell->fCuts, TRUE);
	cell->feasCnt 		= 0;
	cell->infeasIncumb 	= FALSE;

	/* stochastic components */
	if ( config.SAA == 1 ) {
		cnt = cell->omega->cnt;
		if (cell->omega) freeOmegaType(cell->omega, TRUE);
		cell->omega->cnt = cnt;
	}

	/* reset all the clocks */
	cell->time->repTime = cell->time->iterTime = cell->time->masterIter = cell->time->subprobIter = cell->time->optTestIter = 0.0;
	cell->time->iterAccumTime = cell->time->masterAccumTime = cell->time->subprobAccumTime = cell->time->optTestAccumTime = 0.0;

	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return 1;
		}

		cell->incumbChg = FALSE;

#if defined(ALGO_CHECK)
		if ( writeProblem(cell->master->lp, "cleanedQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return 1;
		}
#endif
	}

	return 0;
}//END cleanCellType()

void freeCellType(cellType *cell) {

	if ( cell ) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->candidX) mem_free(cell->candidX);
		if (cell->incumbX) mem_free(cell->incumbX);
		if (cell->cuts) freeCutsType(cell->cuts, FALSE);
		if (cell->fCuts) freeCutsType(cell->fCuts, FALSE);
		if (cell->omega) freeOmegaType(cell->omega, FALSE);
		if (cell->piM) mem_free(cell->piM);
		if (cell->time) mem_free(cell->time);
		mem_free(cell);
	}

}//END freeCellType()
