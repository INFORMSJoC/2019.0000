/*
 * algo.c
 *
 *  Created on: Feb 17, 2019
 *      Author: Jiajun Xu
 * Institution: University of Southern California
 *
 * Please send your comments or bug report to jiajunx (at) usc (dot) edu
 *
 */

#include "benders.h"

extern string outputDir;
extern configType config;

#if defined(SAVE_DUALS)
dualsType *duals = NULL;
#endif

int algo (oneProblem *orig, timeType *tim, stocType *stoc, string probName) {
	probType **prob = NULL;
	cellType *cell = NULL;
	vector 	 meanSol;
    double   std=0, compromise_time=0;
    batchSummary *batch = NULL;
	int 	 rep, m, n, out_idx=0;
	FILE 	*sFile, *iFile = NULL;
    char results_name[BLOCKSIZE];
    char incumb_name[BLOCKSIZE];
	clock_t	tic;


    while (TRUE) {
	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol) )
		goto TERMINATE;

	printf("Starting Benders decomposition.\n");
    sprintf(results_name, "results%d.txt", out_idx);
	sFile = openFile(outputDir, results_name, "w");
//    if ( config.MASTER_TYPE == PROB_QP )
    sprintf(incumb_name, "incumb%d.txt", out_idx);
    iFile = openFile(outputDir, incumb_name, "w");
	printDecomposeSummary(sFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);
    fprintf(sFile, "\n Number of observations in each replications: %d", config.MAX_OBS);

	for ( rep = 0; rep < config.NUM_REPS; rep++ ) {
		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "Replication-%d\n", rep+1);
		fprintf(stdout, "\n====================================================================================================================================\n");
		fprintf(stdout, "Replication-%d\n", rep+1);

		/* setup the seed to be used in the current iteration */
		config.RUN_SEED[0] = config.RUN_SEED[rep+1];
		config.EVAL_SEED[0] = config.EVAL_SEED[rep+1];

		if ( rep != 0 ) {
			/* clean up the cell for the next replication */
			if ( cleanCellType(cell, prob[0], meanSol) ) {
				errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
				goto TERMINATE;
			}
		}

		/* Update omega structure */
		if ( config.SAA ) {
			setupSAA(stoc, &config.RUN_SEED[0], &cell->omega->vals, &cell->omega->probs, &cell->omega->cnt, config.TOLERANCE);
			for ( m = 0; m < cell->omega->cnt; m++ )
				for ( n = 1; n <= stoc->numOmega; n++ )
					cell->omega->vals[m][n] -= stoc->mean[n-1];
		}

//        set initialized points
//        int idx;
//        for(idx=0;idx<=20;idx++){
//            cell->candidX[idx] = 0;
//        }
//        cell->candidX[2] = cell->candidX[3] = cell->candidX[4] = cell->candidX[5] = cell->candidX[8] = cell->candidX[12] = cell->candidX[14] = cell->candidX[19] = cell->candidX[6] = cell->candidX[9] = 1;
//        cell->candidX[0] = 10;
//        printVectorInSparse(cell->candidX, 20, NULL);
//
        
		tic = clock();
		/* Use two-stage algorithm to solve the problem */
		if ( solveBendersCell(stoc, prob, cell) ) {
			errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
			goto TERMINATE;
		}
		cell->time->repTime = ((double) (clock() - tic))/CLOCKS_PER_SEC;

		/* Write solution statistics for optimization process */
		writeStatistic(sFile, iFile, prob, cell);
		writeStatistic(stdout, NULL, prob, cell);

		/* evaluating the optimal solution*/
		if (config.EVAL_FLAG == 1) {
			if ( config.MASTER_TYPE == PROB_QP )
				evaluate(sFile, stoc, prob, cell, cell->incumbX);
			else
				evaluate(sFile, stoc, prob, cell, cell->candidX);
		}
        /* Save the batch details and build the compromise problem. */
        if ( config.MULTIPLE_REP ) {
            buildCompromise(prob[0], cell, batch);
        }
	}
    
    if ( config.MULTIPLE_REP ) {
        /* Solve the compromise problem. */
        tic = clock();
        if ( solveCompromise(prob[0], batch)) {
            errMsg("algorithm", "algo", "failed to solve the compromise problem", 0);
            goto TERMINATE;
        }
        compromise_time = ((double) (clock() - tic))/CLOCKS_PER_SEC;
        batch->time->repTime += compromise_time; 
        
        fprintf(sFile, "\n====================================================================================================================================\n");
        fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
        fprintf(sFile, "\n====================================================================================================================================\n");
        fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
        /* Evaluate the compromise solution */
        fprintf(sFile, "Incumbent solution, non-zero position: ");
        printVectorInSparse(batch->compromiseX, prob[0]->num->cols, sFile);
        evaluate(sFile, stoc, prob, cell, batch->compromiseX);
        
//        //evaluate the test solution if needed:
//        double test[] = {4, 0, 0, 0, 0,1, 0, 0, 0, 0,1,1, 0,1,0, 0, 0, 0,0, 0, 0};
//        evaluate(sFile, stoc, prob, cell, test);
        
        fprintf(sFile, "Time to solve compromise problem   : %f\n", compromise_time);
        fprintf(sFile, "Total time                         : %f\n", batch->time->repTime);
        fprintf(sFile, "Total time to solve master         : %f\n", batch->time->masterAccumTime);
        fprintf(sFile, "Total time to solve subproblems    : %f\n", batch->time->subprobAccumTime);
        
        fprintf(sFile, "Lower bound estimate               : %f\n", batch->Est);
        std = LowerBoundVariance(batch);
        fprintf(sFile, "Lower bound estimation std               : %f\n", std);
        
        fprintf(sFile, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
        fprintf(stdout, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
        /* Evaluate the average solution */
        printVectorInSparse(batch->avgX, prob[0]->num->cols, sFile);
        evaluate(sFile, stoc, prob, cell, batch->avgX);
        fclose(sFile);
        fprintf(iFile, "\n----------------------------------------- Compromise solution(1-indexed) --------------------------------------\n\n");
        printVectorInSparse(batch->compromiseX, prob[0]->num->cols, iFile);
        fprintf(iFile, "\n------------------------------------------- Average solution(1-indexed) ---------------------------------------\n\n");
        printVectorInSparse(batch->avgX, prob[0]->num->cols, iFile);
    }

	 fclose(iFile);
     
        
        /*outer loop condition*/
        
        if (InConvexHull(batch, prob[0]->num->cols) && (std  < config.std_tol)){
            printf("\nSuccessfully completed the L-shaped method.\n");
            break;
        }
        
        config.MAX_OBS = 2 * config.MAX_OBS;
        out_idx++;
    }
	

	/* free up memory before leaving */
	freeCellType(cell);
	freeProbType(prob, 2);
	mem_free(meanSol);
	return 0;

	TERMINATE:
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, 2);
	mem_free(meanSol);
	return 1;
}//END algo()

int solveBendersCell(stocType *stoc, probType **prob, cellType *cell) {
	int 	candidCut;
	clock_t	tic, mainTic;

#if defined(SAVE_DUALS)
	if ( duals == NULL ) {
		duals = (dualsType *) mem_malloc(sizeof(dualsType));
		duals->iter = (intvec) arr_alloc(config.MAX_ITER*cell->omega->cnt, int);
		duals->obs = (intvec) arr_alloc(config.MAX_ITER*cell->omega->cnt, int);
		duals->vals = (vector *) arr_alloc(config.MAX_ITER*cell->omega->cnt, vector);
		duals->cnt = 0;
	}
#endif

	mainTic = clock();
	/* Main loop of the algorithm */
	while (cell->k < config.MAX_ITER) {
		tic = clock();

		cell->k++;

#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
		printf("\nIteration-%d :: Incumbent estimate = %lf; Candidate estimate = %lf.\n", cell->k, cell->incumbEst, cell->candidEst);
#else
		if ( (cell->k-1) % 100 == 0) {
			printf("\nIteration-%4d: ", cell->k);
		}
#endif
		/******* 1a. Optimality tests *******/
//        if ( config.MASTER_TYPE == PROB_QP )
        if (optimal(cell))
            break;

		/******* 2. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formOptCut(prob[1], cell, cell->candidX, FALSE)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}

		/******* 4. Check improvement in predicted values at candidate solution *******/
		if ( config.MASTER_TYPE == PROB_QP ) {
			if ( cell->k > 1 ) {
				/* If the incumbent has not changed in the current iteration */
				checkImprovement(prob[0], cell, candidCut);
			}
			else
				cell->incumbEst = vXvSparse(cell->incumbX, prob[0]->dBar) + cutHeight(cell->cuts->vals[candidCut], cell->incumbX, prob[0]->num->cols);
		}
		else {
			cell->incumbEst = vXvSparse(cell->candidX, prob[0]->dBar) + cutHeight(cell->cuts->vals[candidCut], cell->candidX, prob[0]->num->cols);

			if (optimal(cell))
				break;
		}

		/******* 3. Solve the master problem to obtain the new candidate solution */
		if ( solveMaster(prob[0]->num, prob[0]->dBar, cell) ) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}
		cell->time->masterAccumTime += cell->time->masterIter; cell->time->subprobAccumTime += cell->time->subprobIter;
		cell->time->masterIter = cell->time->subprobIter = cell->time->optTestIter = 0.0;
		cell->time->iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time->iterAccumTime += cell->time->iterTime;

		if ( cell->k % 10 == 0 )
			printf("Time = %lf", ((double) (clock() - mainTic))/CLOCKS_PER_SEC);
	}

	return 0;
}//END solveCell()

BOOL optimal(cellType *cell) {

	if ( cell->RepeatedTime > 0 || cell->k > config.MIN_ITER ) {
        if (config.MASTER_TYPE == PROB_QP || config.reg == 1){
            return cell->optFlag = ((cell->incumbEst - cell->candidEst) < config.EPSILON);
        }
        return TRUE;
	}

	return FALSE;
}//END optimal()

void writeStatistic(FILE *soln, FILE *incumb, probType **prob, cellType *cell) {

	fprintf(soln, "\n------------------------------------------------------------ Optimization ---------------------------------------------------------\n");
	if ( config.MASTER_TYPE == PROB_QP )
		fprintf(soln, "Algorithm                          : Regularized Benders Decomposition\n");
	else
		fprintf(soln, "Algorithm                          : Benders Decomposition\n");
	fprintf(soln, "Number of iterations               : %d\n", cell->k);
	fprintf(soln, "Lower bound estimate               : %f\n", cell->incumbEst);
	if ( cell->k == config.MAX_ITER)
		fprintf(soln, "Maximum itertions reached with gap.");
	fprintf(soln, "Total time                         : %f\n", cell->time->repTime);
	fprintf(soln, "Total time to solve master         : %f\n", cell->time->masterAccumTime);
	fprintf(soln, "Total time to solve subproblems    : %f\n", cell->time->subprobAccumTime);
	fprintf(soln, "Total time in verifying optimality : %f\n", cell->time->optTestAccumTime);

    //print candidate for a moment here
	if ( incumb != NULL ) {
        printVectorInSparse(cell->candidX, prob[0]->num->cols, incumb);
		//printVector(cell->candidX, prob[0]->num->cols, incumb);
	}

//    if ( duals != NULL ) {
//        printf("\nNumber of duals discovered = %d", duals->cnt);
//    }

}//END WriteStat

