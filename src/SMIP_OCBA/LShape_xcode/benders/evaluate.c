/*/*
 * evaluate.c
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

int evaluate(FILE *soln, stocType *stoc, probType **prob, cellType *cell, vector Xvect) {
	vector 	observ, rhs, costTemp, cost;
	intvec	objxIdx;
	double 	obj, mean, variance, stdev, temp;
	int		cnt, m, status;

	if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
		errMsg("allocation", "evaluateOpt", "observ", 0);

	printf("\nStarting evaluation.\n");

	/* initialize parameters used for evaluations */
	cnt = 0.0; mean = 0.0; variance = 0.0; stdev = INFBOUND; cnt = 0;

	/* right-hand side */
	if (!(rhs =(vector) arr_alloc(prob[1]->num->rows+1, double)))
		errMsg("Allocation", "evaluate", "rhs",0);

	/* cost coefficients */
	if ( !(cost = (vector) arr_alloc(prob[1]->num->cols+1, double)) )
		errMsg("allocation", "evaluate", "cost", 0);
	if ( !(objxIdx = (intvec) arr_alloc(prob[1]->num->cols+1, int)) )
		errMsg("allocation", "evaluate", "objxIdx", 0);
	costTemp = expandVector(prob[1]->dBar->val, prob[1]->dBar->col, prob[1]->dBar->cnt, prob[1]->num->cols);
	for (m = 1; m <= prob[1]->num->rvdOmCnt; m++ ) {
		objxIdx[m] = prob[1]->coord->rvdOmCols[m] - 1;
		cost[m] = costTemp[objxIdx[m]+1];
	}
	mem_free(costTemp);

	/* change the right hand side with the solution */
	chgRHSwSoln(prob[1]->bBar, prob[1]->Cbar, rhs, Xvect);

	while (3.92 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
		/* use the stoc file to generate observations */
		generateOmega(stoc, observ, config.TOLERANCE, &config.EVAL_SEED[0]);

		for ( m = 0; m < stoc->numOmega; m++ )
			observ[m] -= stoc->mean[m];

		/* Change right-hand side with random observation */
		if ( chgRHSwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, observ-1, rhs, Xvect) ) {
			errMsg("algorithm", "evaluate", "failed to change right-hand side with random observations",0);
			return 1;
		}

		/* Change cost coefficients with random observations */
		if ( prob[1]->num->rvdOmCnt > 0 ) {
			if ( chgObjxwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, cost, objxIdx, observ-1) ) {
				errMsg("algorithm", "evaluate","failed to change cost coefficients with random observations", 0);
				return 1;
			}
		}

		if ( solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &status) ) {
			if ( status == STAT_INFEASIBLE ) {
				/* subproblem is infeasible */
				printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
				return 1;
			}
			else {
				errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
				return 1;
			}
		}

		/* use subproblem objective and compute evaluation statistics */
		obj = getObjective(cell->subprob->lp, PROB_LP);

#if defined(ALGO_CHECK)
		writeProblem(cell->subprob->lp, "evalSubprob.lp");
		printf("Evaluation objective function = %lf.\n", obj);
#endif


		if ( cnt == 0 )
			mean = obj;
		else {
			temp = mean;
			mean = mean + (obj - mean) / (double) (cnt + 1.0);
			variance  = cnt / (cnt + 1.0) * variance
					+ cnt * (mean - temp) * (mean - temp);
			stdev = sqrt(variance/ (double) cnt);
		}
		cnt++;

		/* Print the results every once in a while for long runs */
		if (!(cnt % 100)) {
			printf(".");
			fflush(stdout);
		}
		if (!(cnt % 10000))
			printf("\nObs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", cnt, mean, 3.92 * stdev / mean,  mean - 1.96  * stdev, mean + 1.96  * stdev);
	}//END while loop
	mean += vXvSparse(Xvect, prob[0]->dBar);

	printf("\n\nEvaluation complete. Final evaluation results :: \n");
	printf("Upper bound estimate               : %lf\n", mean);
	printf("Error in estimation                : %lf\n", 3.92 * stdev / mean);
	printf("Confidence interval at 95%%         : [%lf, %lf]\n", mean - 1.96  * stdev, mean + 1.96  * stdev);
	printf("Number of observations             : %d\n", cnt);

	if ( soln != NULL ) {
		/* Write the evaluation results to the summary file */
		fprintf(soln, "------------------------------------------------------------- Evaluation ----------------------------------------------------------\n");
		fprintf(soln, "Upper bound estimate           : %lf\n", mean);
		fprintf(soln, "Error in estimation            : %lf\n", 3.92 * stdev / mean);
		fprintf(soln, "Confidence interval at 95%%     : [%lf, %lf]\n", mean - 1.96  * stdev, mean + 1.96  * stdev);
		fprintf(soln, "Number of observations         : %d\n", cnt);
	}

	mem_free(observ); mem_free(rhs); mem_free(objxIdx); mem_free(cost);
	return 0;

}//END evaluate()

int evaluate_samples(FILE *soln, stocType *stoc, probType **prob, cellType *cell, vector Xvect) {
    vector     observ, rhs, costTemp, cost;
    intvec    objxIdx;
    double     obj, mean, variance, stdev, temp;
    int        cnt, m, status;

    if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
        errMsg("allocation", "evaluateOpt", "observ", 0);

    printf("\nStarting evaluation.\n");

    /* initialize parameters used for evaluations */
    cnt = 0.0; mean = 0.0; variance = 0.0; stdev = INFBOUND; cnt = 0;

    /* right-hand side */
    if (!(rhs =(vector) arr_alloc(prob[1]->num->rows+1, double)))
        errMsg("Allocation", "evaluate", "rhs",0);

    /* cost coefficients */
    if ( !(cost = (vector) arr_alloc(prob[1]->num->cols+1, double)) )
        errMsg("allocation", "evaluate", "cost", 0);
    if ( !(objxIdx = (intvec) arr_alloc(prob[1]->num->cols+1, int)) )
        errMsg("allocation", "evaluate", "objxIdx", 0);
    costTemp = expandVector(prob[1]->dBar->val, prob[1]->dBar->col, prob[1]->dBar->cnt, prob[1]->num->cols);
    for (m = 1; m <= prob[1]->num->rvdOmCnt; m++ ) {
        objxIdx[m] = prob[1]->coord->rvdOmCols[m] - 1;
        cost[m] = costTemp[objxIdx[m]+1];
    }
    mem_free(costTemp);

    /* change the right hand side with the solution */
    chgRHSwSoln(prob[1]->bBar, prob[1]->Cbar, rhs, Xvect);

    while (3.92 * stdev > config.EVAL_ERROR * DBL_ABS(mean) || cnt < config.EVAL_MIN_ITER ) {
        /* use the stoc file to generate observations */
        generateOmega(stoc, observ, config.TOLERANCE, &config.EVAL_SEED[0]);

        for ( m = 0; m < stoc->numOmega; m++ )
            observ[m] -= stoc->mean[m];

        /* Change right-hand side with random observation */
        if ( chgRHSwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, observ-1, rhs, Xvect) ) {
            errMsg("algorithm", "evaluate", "failed to change right-hand side with random observations",0);
            return 1;
        }

        /* Change cost coefficients with random observations */
        if ( prob[1]->num->rvdOmCnt > 0 ) {
            if ( chgObjxwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, cost, objxIdx, observ-1) ) {
                errMsg("algorithm", "evaluate","failed to change cost coefficients with random observations", 0);
                return 1;
            }
        }

        if ( solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &status) ) {
            if ( status == STAT_INFEASIBLE ) {
                /* subproblem is infeasible */
                printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
                return 1;
            }
            else {
                errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
                return 1;
            }
        }

        /* use subproblem objective and compute evaluation statistics */
        obj = getObjective(cell->subprob->lp, PROB_LP);

#if defined(ALGO_CHECK)
        writeProblem(cell->subprob->lp, "evalSubprob.lp");
        printf("Evaluation objective function = %lf.\n", obj);
#endif


        if ( cnt == 0 )
            mean = obj;
        else {
            temp = mean;
            mean = mean + (obj - mean) / (double) (cnt + 1.0);
            variance  = cnt / (cnt + 1.0) * variance
                    + cnt * (mean - temp) * (mean - temp);
            stdev = sqrt(variance/ (double) cnt);
        }
        cnt++;

        /* Print the results every once in a while for long runs */
        if (!(cnt % 100)) {
            printf(".");
            fflush(stdout);
        }
        if (!(cnt % 10000))
            printf("\nObs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", cnt, mean, 3.92 * stdev / mean,  mean - 1.96  * stdev, mean + 1.96  * stdev);
    }//END while loop
    mean += vXvSparse(Xvect, prob[0]->dBar);

    printf("\n\nEvaluation complete. Final evaluation results :: \n");
    printf("Upper bound estimate               : %lf\n", mean);
    printf("Error in estimation                : %lf\n", 3.92 * stdev / mean);
    printf("Confidence interval at 95%%         : [%lf, %lf]\n", mean - 1.96  * stdev, mean + 1.96  * stdev);
    printf("Number of observations             : %d\n", cnt);

    if ( soln != NULL ) {
        /* Write the evaluation results to the summary file */
        fprintf(soln, "------------------------------------------------------------- Evaluation ----------------------------------------------------------\n");
        fprintf(soln, "Upper bound estimate           : %lf\n", mean);
        fprintf(soln, "Error in estimation            : %lf\n", 3.92 * stdev / mean);
        fprintf(soln, "Confidence interval at 95%%     : [%lf, %lf]\n", mean - 1.96  * stdev, mean + 1.96  * stdev);
        fprintf(soln, "Number of observations         : %d\n", cnt);
    }

    mem_free(observ); mem_free(rhs); mem_free(objxIdx); mem_free(cost);
    return cnt;

}//END evaluate()


