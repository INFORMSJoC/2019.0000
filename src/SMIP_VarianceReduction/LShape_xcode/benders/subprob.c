/*
 * subprob.c
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
extern string outputDir;

int solveSubprob(probType *prob, oneProblem *subproblem, vector Xvect, vector obsVals, BOOL *spFeasFlag, double *subprobTime, vector piS, double *mubBar) {
	vector 	rhs, cost;
	intvec	indices;
	int  	status, n;
	clock_t	tic;

#ifdef DEBUG
    printf("prob->num->rows = %d, prob->num->cols = %d \n", prob->num->rows, prob->num->cols);
#endif
	if ( !(indices = (intvec) arr_alloc(max(prob->num->rows, prob->num->cols), int)) )
		errMsg("allocation", "solve_subporb", "indices", 0);
	for ( n = 0; n < max(prob->num->rows,prob->num->cols); n++ )
		indices[n] = n;

	/* (a) compute the right-hand side using current observation and first-stage solution */
	rhs = computeRHS(prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, obsVals);
	if ( rhs == NULL ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
		return 1;
	}

	/* (b) change the right-hand side in the solver */
	if ( changeRHS(subproblem->lp, prob->num->rows, indices, rhs + 1) ) {
		errMsg("solver", "solve_subprob", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	/* (c) compute the cost coefficients using current observation */
	cost = computeCostCoeff(prob->num, prob->coord, prob->dBar, obsVals);
	if ( cost == NULL ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem cost coefficients", 0);
		return 1;
	}

#ifdef DEBUG
    printf("cost vector (%d elements):\n", prob->num->cols);
    printVector(cost, prob->num->cols, NULL);
#endif
	/* (d) change cost coefficients in the solver */
	if ( changeObjx(subproblem->lp, prob->num->cols, indices, cost+1) ) {
		errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver",0);
		return 1;
	}

	tic = clock();
	/* (e) Solve the subproblem to obtain the optimal dual solution. */
	if ( solveProblem(subproblem->lp, subproblem->name, subproblem->type, &status) ) {
		if ( status == STAT_INFEASIBLE ) {
			printf("Subproblem is infeasible: need to create feasibility cut.\n");
			(*spFeasFlag) = FALSE;
		}
		else {
			errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
			return 1;
		}
	}
	(*subprobTime) += ((double) (clock() - tic))/CLOCKS_PER_SEC;

#if defined(ALGO_CHECK)
	writeProblem(subproblem->lp, "cellSubprob.lp");
#endif

#if defined(ALGO_CHECK)
	double obj;
	obj = getObjective(subproblem->lp, PROB_LP);
	printf("Objective value of Subproblem  = %lf;\t", obj);
#endif

	/* Record the dual and reduced cost on bounds. */
	if ( getDual(subproblem->lp, piS, prob->num->rows) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
		return 1;
	}

#if 0
	FILE *basisFile;
	intvec cstat, rstat;
	cstat = (intvec) arr_alloc( prob->num->cols+1, int);
	rstat = (intvec) arr_alloc( prob->num->rows+1, int);

	/* Obtain the status of columns and rows in the basis. */
	if ( getBasis(subproblem->lp, cstat, rstat) ) {
		errMsg("algorithm", "newBasis", "failed to get the basis column and row status", 0);
		return 1;
	}

	basisFile = openFile(outputDir, "cstat.txt", "a");
	printIntvec(cstat, prob->num->cols, basisFile);
	fclose(basisFile);

	basisFile = openFile(outputDir, "rstat.txt", "a");
	printIntvec(cstat, prob->num->rows, basisFile);
	fclose(basisFile);

	mem_free(cstat); mem_free(rstat);
#endif


	if ( computeMU(subproblem->lp, prob->num->cols, mubBar) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to compute mubBar for subproblem", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(cost);
	mem_free(indices);
	return 0;
}// END solveSubprob()

/* This function computes the right hand side of the subproblem, based on a given X vector and a given observation of omega.
 * It is defined as:
 * 			rhs = R(omega) - T(omega) x X
 * and is calculated as:
 * 			rhs = (Rbar - Tbar x X) + (Romega - Tomega x X)
 *
 * where the "bar" denotes the fixed or mean value, and the "omega" denotes a random variation from this mean. The function allocates an array
 * for the vector, which must be freed by the customer.  Also, the zeroth position of this rhs vector is reserved, and the actual values begin at rhs[1].
 * R is b, and T is C
 \***********************************************************************/
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector observ) {
	int cnt;
	vector rhs;
	sparseVector bomega;
	sparseMatrix Comega;

	

	/* Start with the values of b -- fixed */
	rhs = expandVector(bBar->val, bBar->col, bBar->cnt, num->rows);

	/* (cumulatively) subtract values of C x X -- fixed */
	rhs = MSparsexvSub(Cbar, X, rhs);
    
    /*get b(omega) and C(omega)*/
    if(num->rvRowCnt > 0){
        bomega.cnt = num->rvbOmCnt;    bomega.col = coord->rvbOmRows; bomega.val = coord->rvOffset[0] + observ;
        Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols; Comega.row = coord->rvCOmRows; Comega.val = coord->rvOffset[1] + observ;
        for (cnt = 1; cnt <= bomega.cnt; cnt++)
            rhs[bomega.col[cnt]] += bomega.val[cnt];
        rhs = MSparsexvSub(&Comega, X, rhs);
    }

	return rhs;
}//END computeRHS()

vector computeCostCoeff(numType *num, coordType *coord, sparseVector *dBar, vector observ) {
	vector cost;
	sparseVector dOmega;
	int	cnt;

	dOmega.cnt = num->rvdOmCnt; dOmega.col = coord->rvdOmCols; dOmega.val = coord->rvOffset[2] + observ;

	cost = expandVector(dBar->val, dBar->col, dBar->cnt, num->cols);
    for (cnt = 1; cnt <= dOmega.cnt; cnt++){
		cost[dOmega.col[cnt]] += dOmega.val[cnt];
//#ifdef DEBUG
//        printf("cnt = %d, dOmega.val = %f, cost[dOmega.col[cnt]]=%f, dOmega.col[cnt]=%d \n", cnt, dOmega.val[cnt], cost[dOmega.col[cnt]], dOmega.col[cnt]);
//#endif
    }

	return cost;
}//END computeCostCoeff()

/* This function compute the reduced cost of every second stage variables. They will be used to calculate the \mu x b and then added to the \pi x b. */
int computeMU(LPptr lp, int numCols, double *mubBar) {
	vector	dj, u;
	intvec	cstat;
	int		n;

	(*mubBar) = 0.0;

	if ( !(dj = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "dual slacks", 0);
	if ( !(u = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "TDA solutions", 0);
	if ( !(cstat = (intvec) arr_alloc( numCols+1, int)))
		errMsg("allocation", "stochasticUpdates", "cstat", 0);

	if ( getPrimal(lp, u, numCols) ) {
		errMsg("solver", "forOptPass", "failed to obtain primal solution", 0);
		return 1;
	}
	if (getDualSlacks(lp, dj, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	for (n = 1; n <= numCols;  n++) {
		switch (cstat[n]) {
		case AT_LOWER:
			(*mubBar) += dj[n]*u[n];
			break;
		case AT_UPPER:
			(*mubBar) += dj[n]*u[n];
			break;
		default:
			break;
		}
	}

	mem_free(u); mem_free(dj); mem_free(cstat);
	return 0;
}//END compute_mu()

oneProblem *newSubproblem(oneProblem *subprob) {

	/* since the basic structure of subproblem is not modified during the course of the algorithm, we just load it onto the solver */
	subprob->lp = setupProblem(subprob->name, subprob->type, subprob->mac, subprob->mar, subprob->objsen, subprob->objx, subprob->rhsx, subprob->senx,subprob->matbeg, subprob->matcnt, subprob->matind, subprob->matval, subprob->bdl, subprob->bdu, NULL, subprob->cname, subprob->rname, subprob->ctype);
	if ( subprob->lp == NULL ) {
		errMsg("Problem Setup", "new_subprob", "subprob",0);
		return NULL;
	}

#if defined(SETUP_CHECK)
	if (writeProblem(subprob->lp, "newSubproblem.lp") ){
		errMsg("solver", "newSubproblem", "failed to write subproblems to file", 0);
		return NULL;
	}
#endif

	return subprob;
}//END new_subprob

void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, vector rhs, vector X) {
	int cnt;

	/* copy the original right-hand side */
	for (cnt = 1; cnt <= bBar->cnt; cnt++)
		rhs[bBar->col[cnt]] = bBar->val[cnt];

	/* change the right-hand side with first stage solution */
	rhs = MSparsexvSub(Cbar, X, rhs);

}//END chgRHSwMean()

int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, vector observ, vector spRHS, vector X) {
	sparseVector bomega;
	sparseMatrix Comega;
	vector 	rhs;
	intvec	indices;
	int		cnt, stat1;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->rvbOmRows; bomega.val = coord->rvOffset[0]+observ;
	Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols; Comega.row = coord->rvCOmRows; Comega.val = coord->rvOffset[2] + observ;

	if ( !(indices = (intvec) arr_alloc(num->rows, int)) )
		errMsg("allocation", "chgRHSwRand", "indices", 0);
	if ( !(rhs = (vector) arr_alloc(num->rows+1, double)) )
		errMsg("allocation", "chgRHSwRand", "rhs", 0);

	/* copy right-hand side modified with mean information */
	for ( cnt = 1; cnt <= num->rows; cnt++ ) {
		rhs[cnt] = spRHS[cnt];
		indices[cnt-1] = cnt-1;
	}

	/* change right-hand side with randomness in b */
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* change right-hand side with randomness in transfer matrix */
	rhs = MSparsexvSub(&Comega, X, rhs);

	/* change the right-hand side in the solver */
	stat1 = changeRHS(lp, num->rows, indices, rhs + 1);
	if ( stat1 ) {
		errMsg("solver", "chgRHSwRand", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	mem_free(rhs); mem_free(indices);
	return 0;

}//END chgRHSwRand()

int chgObjxwObserv(LPptr lp, numType *num, coordType *coord, vector cost, intvec indices, vector observ) {
	vector vals;
	int n;

	if ( !(vals = (vector) arr_alloc(num->rvdOmCnt+1, double)) )
		errMsg("allocation", "chgObjwObserv", "vals", 0);

	for ( n = 1; n <= num->rvdOmCnt; n++ )
		vals[n] = cost[n] + observ[coord->rvOffset[2]+n];

	if ( changeObjx(lp, num->rvdOmCnt, indices+1, vals+1) ) {
		errMsg("solver", "chgObjswObserv", "failed to change the cost coefficients in the solver",0);
		return 1;
	}

	mem_free(vals);
	return 0;
}//END chgObjwObserv()

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a vector to hold an array of
 * observation and the probability associated with it. */
omegaType *newOmega(stocType *stoc) {
	omegaType *omega;
	int cnt, i, base, idx;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation","newOmega", "omega", 0);
	if ( !(omega->probs = (vector) arr_alloc(config.MAX_OBS, double)) )
		errMsg("allocation", "newOmega", "omega->probs", 0);
	if ( !(omega->vals = (vector *) arr_alloc(config.MAX_OBS, vector)) )
		errMsg("allocation", "newOmega", "omega->vals", 0);
	omega->cnt = 0; omega->numRV = stoc->numOmega;

	if ( config.SAA == 1 ) {
		omega->cnt = config.MAX_OBS;
		return omega;
	}

	if ( strstr(stoc->type, "BLOCKS") != NULL ) {
		if ( (omega->cnt = stoc->numVals[0]) <= config.MAX_OBS) {
			omega->vals = (vector *) mem_realloc(omega->vals, omega->cnt*sizeof(vector));
			omega->probs = (vector) mem_realloc(omega->probs, omega->cnt*sizeof(double));
			for ( cnt = 0; cnt < omega->cnt; cnt++) {
				omega->probs[cnt]= stoc->probs[0][cnt];
				if ( !(omega->vals[cnt] = (vector) arr_alloc(omega->numRV+1, double)) )
					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
				for (i = 0; i < omega->numRV; i++)
					omega->vals[cnt][i+1]=stoc->vals[i][cnt]-stoc->mean[i];
				omega->vals[cnt][0] = oneNorm(omega->vals[cnt]+1, omega->numRV);
			}
		}
		else {
			omega->cnt = config.MAX_OBS;
			config.SAA = 1;
		}
	}
	else if ( strstr(stoc->type, "INDEP") != NULL ) {
		omega->cnt = 1; i = 0;
		while ( i < stoc->numOmega ) {
			omega->cnt *= stoc->numVals[i];
			if (omega->cnt > config.MAX_OBS) {
				omega->cnt = config.MAX_OBS;
				config.SAA = 1;
				break;
			}
			i++;
		}

		if ( !config.SAA ){
			omega->vals = (vector *) mem_realloc(omega->vals, omega->cnt*sizeof(vector));
			omega->probs = (vector) mem_realloc(omega->probs, omega->cnt*sizeof(double));
			for ( cnt = 0; cnt < omega->cnt; cnt++) {
				if ( !(omega->vals[cnt] = (vector) arr_alloc(omega->numRV+1, double)) )
					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
				omega->probs[cnt] = 1; base = omega->cnt;
				for ( i = 0; i < omega->numRV; i++ ) {
					base /= stoc->numVals[i];
					idx = (int)((double) cnt / (double) base) % stoc->numVals[i];
					omega->vals[cnt][i+1] = stoc->vals[i][idx]-stoc->mean[i];
					omega->probs[cnt] *= stoc->probs[i][idx];
				}
			}
		}
	}
	else {
		omega->cnt = config.MAX_OBS;
		config.SAA = 1;
	}

	return omega;
}//END newOmega()

void freeOmegaType(omegaType *omega, BOOL partial) {
	int n;

	if ( omega->vals ) {
		for ( n = 0; n < omega->cnt; n++ )
			if ( omega->vals[n] )
				mem_free(omega->vals[n]);
		if ( partial ) {
			omega->cnt = 0;
			return;
		}
		mem_free(omega->vals);
	}
	if ( omega->probs ) mem_free(omega->probs);
	mem_free(omega);

}//END freeOmegaType()
