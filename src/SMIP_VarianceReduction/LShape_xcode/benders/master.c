/*
 * master.c
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

/* In this his function the master problem is solved after the newest cut is added to master problem, the incumbent cut is updated if necessary.
 * Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveMaster(numType *num, sparseVector *dBar, cellType *cell) {
	int 	status, i;
	clock_t	tic;
    vector new_X;
    BOOL sameIndicator = TRUE;

#if defined(ALGO_CHECK)
	writeProblem(cell->master->lp,"cellMaster.lp");
#endif

	tic = clock();
	/* solve the master problem */
	if ( solveProblem(cell->master->lp, cell->master->name, config.MASTER_TYPE, &status) ) {
		writeProblem(cell->master->lp, "error.lp");
		errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
		return 1;
	}
	cell->time->masterIter = ((double) (clock() - tic))/CLOCKS_PER_SEC;

	/* increment the number of problems solved during algorithm */
	cell->LPcnt++;

    if (!(new_X =(vector) arr_alloc(num->cols+1, double)))
        errMsg("Allocation", "solveMaster", "last_x",0);
   
	/* Get the most recent optimal solution to master program */
	if ( getPrimal(cell->master->lp, new_X, num->cols) ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}
    for(i = 0; i <= num->cols; i++){
        if (DBL_ABS(cell->candidX[i] - new_X[i]) > config.TOLERANCE){
            sameIndicator = FALSE;
            break;
        }
    }
    if (sameIndicator){
        cell->RepeatedTime++;
    }
    else{
        cell->RepeatedTime = 0;
        for(i = 0; i <= num->cols; i++){
            cell->candidX[i] = new_X[i];
        }
    }
    

	if ( cell->master->type == PROB_QP ) {
		/* Get the dual solution too */
		if ( getDual(cell->master->lp, cell->piM, cell->master->mar) ) {
			errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
			return 1;
		}

		addVectors(cell->candidX, cell->incumbX, NULL, num->cols);
		cell->candidEst = vXvSparse(cell->candidX, dBar) + getPrimalPoint(cell->master->lp, num->cols);
	}
    else if(cell->master->type == PROB_MILP){
        cell->candidEst = getObjective(cell->master->lp, PROB_MILP);
    }
    else{
		cell->candidEst = getObjective(cell->master->lp, PROB_LP);
    }
    free(new_X);
	return 0;
}//END solveMaster()


int addCut2Master(cellType *cell, cutsType *cuts, oneCut *cut, int lenX) {
	intvec 	indices;
	int 	cnt;
	static int cummCutNum = 0;

	/* If it is optimality cut being added, check to see if there is room for the candidate cut, else drop a cut */
	if (cuts->cnt == cell->maxCuts) {
		/* make room for the latest cut */
		if( reduceCuts(cell->master, cuts, cell->candidX, cell->piM, lenX, &cell->iCutIdx) < 0 ) {
			errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
			return -1;
		}
	}

//    if ( config.MASTER_TYPE == PROB_QP )
//        cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell->incumbX, NULL, lenX);

	if (!(indices = arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;
	indices[0] = lenX;

	/* Add the cut to the cell cuts structure and assign a row number. */
	cuts->vals[cuts->cnt] = cut;
	cut->rowNum = cell->master->mar++;

	/* Set up the cut name */
	sprintf(cut->name, "cut_%04d", cummCutNum++);

	/* Add the row in the solver */
	if ( addRow(cell->master->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta, cut->name) ) {
		errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
		return -1;
	}

	mem_free(indices);
	return cuts->cnt++;
}//END addCuts2Master()


int checkImprovement(probType *prob, cellType *cell, int candidCut) {

	/* Calculate height at new candidate x with newest cut included */
	cell->candidEst = vXvSparse(cell->candidX, prob->dBar) + cutHeight(cell->cuts->vals[candidCut], cell->candidX, prob->num->cols);

#if defined(ALGO_CHECK)
	printf("Incumbent estimate = %lf; Candidate estimate = %lf\n", cell->incumbEst, cell->candidEst);
#endif

	/* If we see considerable improvement, then change the incumbent */
	if ( cell->candidEst <= cell->incumbEst ) {
		/* when we find an improvement, then we need to replace the incumbent x with candidate x */
		if ( replaceIncumbent(prob, cell) ) {
			errMsg("algorithm", "checkImprovement", "failed to replace incumbent solution with candidate", 0);
			return 1;
		}
		cell->iCutIdx = candidCut;
		cell->incumbChg = FALSE;
		printf("+"); fflush(stdout);
	}

	return 0;
}//END checkImprovement()

int replaceIncumbent(probType *prob, cellType *cell) {

	/* replace the incumbent solution with the candidate solution */
	copyVector(cell->candidX, cell->incumbX, prob->num->cols, 1);
	cell->incumbEst = cell->candidEst;

	/* update the right-hand side and the bounds with new incumbent solution */
	if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
		return 1;
	}

	/* update the candidate cut as the new incumbent cut */
	cell->incumbChg = TRUE;

	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	cell->infeasIncumb = FALSE;

	return 0;
}//END replaceIncumbent()


int constructQP(probType *prob, cellType *cell, vector incumbX, double quadScalar) {
	int status;

	status = changeQPproximal(cell->master->lp, prob->num->cols, quadScalar);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
		return 1;
	}
	status = changeQPrhs(prob, cell, incumbX);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the right-hand side to convert the problem into QP", 0);
		return 1;
	}
	status = changeQPbds(cell->master->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, incumbX);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into QP", 0);
		return 1;
	}

	return status;
}//END constructQP()

//Jiajun function for construction MILP, MIP with L1 penalty in the master problem
int constructMILP(probType *prob, cellType *cell, vector incumbX, double L1Scalar) {
    
    //    if ( changeMILPwithL1(cell->master->lp, prob->sp, prob->num->cols) ) {
    //        errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
    //        return 1;
    //    }
    
    //#ifdef DEBUG
    //    writeProblem(cell->master->lp, "After_changeMILPwithL1.lp");
    //#endif
    
    if ( changeMILPproximal(cell->master->lp, prob->sp->objx, prob->num->cols, L1Scalar) ) {
        errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
        return 1;
    }
    
    if ( changeMILPrhs(prob, cell, incumbX) ) {
        errMsg("algorithm", "algoIntSD", "failed to change the right-hand side to convert the problem into QP", 0);
        return 1;
    }
    
    //    writeProblem(cell->master->lp, "After_constructMILP.lp");
    
    if ( changeMILPbds(cell->master->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, incumbX, 0) ) {
        errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into MILP", 0);
        return 1;
    }
    return 0;
}//END constructMILP()

int changeMILPwithL1(LPptr lp, oneProblem *sp, int numCols){
    int n, status;
    vector lb, ub;
    vector new_objx, new_matval;
    int indices[numCols];
    char ctype[numCols];
    
    //    printf("Number of column in probType %d, number of column in lp %d ", numCols, );
    
    if (!(new_objx = arr_alloc(numCols+1, double)))
        errMsg("Allocation", "changeMILPwithL1", "new_objx",0);
    if (!(new_matval = arr_alloc(numCols+1, double)))
        errMsg("Allocation", "changeMILPwithL1", "new_matval",0);
    if (!(lb = arr_alloc(numCols+1, double)))
        errMsg("Allocation", "changeMILPwithL1", "lb",0);
    if (!(ub = arr_alloc(numCols+1, double)))
        errMsg("Allocation", "changeMILPwithL1", "ub",0);
    
    //    status = getObj(lp, objs, 0, numCols-1);
    //    if (status)    {
    //        errMsg("solver", "changeMILPwithL1", "failed to get the objective cost coeffs in the solver", 0);
    //        return 1;
    //    }
    /*construct d_minus*/
    for(n = 0; n < numCols; n++){
        lb[n] = 0;
        ub[n] = 1;
        new_objx[n] = -sp->objx[n];
        new_matval[n] = -sp->matval[n];
        indices[n] = numCols + 1 + n;
        ctype[n] = 'B';
        // printf("n= %d, new_matbeg, sp->matind, new_matval=%d, %d, %f \n", n, new_matbeg[n], sp->matind[n], new_matval[n]);
    }
    
    
    status = addcols(lp, numCols, sp->numnz, new_objx, sp->matbeg, sp->matind, new_matval, lb, ub, NULL);
    
    if (status)    {
        errMsg("solver", "changeMILPwithL1", "failed to add d- columns in the solver", 0);
        return 1;
    }
    
    
    changeCtype(lp, numCols, indices, ctype);
    
#ifdef DEBUG
    writeProblem(lp, "changeMILPwithL1.lp");
#endif
    
    mem_free(new_objx);
    mem_free(lb);
    mem_free(ub);
    mem_free(new_matval);
    return 0;
}

//Jiajun todo: add L1 here, replace copyQPseparable()
/* Jiajun: In this function, we consider the first stage only contains binary variables.
 update the variables, and add L1 penalty. Replace the x = \hat x + d, where \hat x is the incubent solution, d is the change. The original objective function is c^T*x, and it changes to c^T*d + sigma*|d|  (we minus the constant term c^T * \hat x, and add the L1 penalty). Then we replace d = d_plus - d_minus, then the absolute value of d, |d| = d_plus + d_minus, where d_plus and d_minus are binary variables. To formulate this, we update the cost coeff of x (new cost coeff = old cost coeff + sigma), and use it to represent d_plus. We add new columns for d_minus, and use the cost coeff is -c, constraint coeff -a.  */
int changeMILPproximal(LPptr lp, vector objx, int numCols, double sigma) {
    int    n;
    vector objx_with_L1;  //the array of d+ and d-;
    intvec indices;
    
    if (!(objx_with_L1 = arr_alloc(2 * numCols, double)))
        errMsg("Allocation", "changeMILPproximal", "objx_with_L1",0);
    if (!(indices = arr_alloc(2 * numCols, int)))
        errMsg("Allocation", "changeMILPproximal", "indices",0);
    
    /*0 ~ numCols - 1  is d+; numCols + 1 ~ 2 * numCols + 1 is d-*/
    for(n = 0; n < numCols; n++){
        objx_with_L1[n] = objx[n] + sigma;
        indices[n] = n;
    }
    for(n = 0; n < numCols; n++){
        objx_with_L1[numCols + n] = -objx[n] + sigma;
        indices[numCols + n] = numCols + n + 1;
    }
    
    changeObjx(lp, 2 * numCols, indices, objx_with_L1);
    
    mem_free(objx_with_L1);
    mem_free(indices);
    return 0;
}//END constructMILPproximal

/* The following function is the same with changeQPrhs() */
int changeMILPrhs(probType *prob, cellType *cell, vector xk) {
    int     status = 0, cnt;
    vector     rhs;
    intvec     indices;
    
    if (!(rhs =(vector) arr_alloc(prob->num->rows+cell->cuts->cnt+1, double)))
        errMsg("Allocation", "changeRhs", "rhs",0);
    if (!(indices =(intvec) arr_alloc(prob->num->rows+cell->cuts->cnt, int)))
        errMsg("Allocation", "changeRhs", "indices",0);
    /* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
     vector, and the x vector. */
    for (cnt = 0; cnt < prob->num->rows; cnt++) {
        rhs[cnt + 1] = prob->sp->rhsx[cnt];
        indices[cnt] = cnt;
    }
    
    /* b - A * xbar */
    rhs = MSparsexvSub(prob->Dbar, xk, rhs);
    
    /*** new rhs = alpha - beta * xbar (benders cuts)***/
    for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
        rhs[prob->num->rows+cnt+1] = cell->cuts->vals[cnt]->alpha - vXv(cell->cuts->vals[cnt]->beta, xk, NULL, prob->sp->mac);
        indices[prob->num->rows+cnt] = cell->cuts->vals[cnt]->rowNum;
        
        cell->cuts->vals[cnt]->alphaIncumb = rhs[prob->num->rows+cnt+1];
    }
    
    /* Now we change the right-hand of the master problem. */
    status = changeRHS(cell->master->lp, prob->num->rows + cell->cuts->cnt, indices, rhs+1);
    if (status)    {
        errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
        return 1;
    }
    
    mem_free(rhs);
    mem_free(indices);
    return 0;
}//END changeMILPrhs()

/* Jiajun: this function only consider binary variables: -\hat x <= d <= 1- \hat x*/
int changeMILPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk, int offset) {
    int     status = 0, cnt;
    vector    bounds;
    intvec    indices;
    char     *lu;
    
    if (!(bounds = arr_alloc(2 * numCols, double)))
        errMsg("Allocation", "changeBounds", "bounds",0);
    if (!(indices = arr_alloc(2 * numCols, int)))
        errMsg("Allocation", "change_bounds", "indices",0);
    if (!(lu = arr_alloc(2 * numCols, char)))
        errMsg("Allocation", "changeBounds", "lu",0);
    
    //printf("xk: ");
    /* Change the Bound, now we have 2 * numCols of variables(numCols of d+/-). If \hat x[i] = 0: d-[i]=0; \hat x[i] = 1: d+[i]=0  */
    for (cnt = 0; cnt < numCols; cnt++) {
        lu[cnt] = 'B';
        bounds[cnt] = 0;
        if(DBL_ABS(xk[cnt+1] - 0) < 0.0001)
            indices[cnt] = numCols + 1 + cnt;
        else
            indices[cnt] = cnt;
        //printf("%f ", xk[cnt]);
    }
    
    status = changeBDS(lp, numCols, indices, lu, bounds);
    if (status) {
        errMsg("algorithm", "changeMILP", "failed to change the bound in the solver", 0);
        return 1;
    }
    
    mem_free(bounds); mem_free(indices); mem_free(lu);
    return 0;
}//END changeMILPbds()

/* Construct the Q diagonal matrix and copy it for quadratic problem. */
int changeQPproximal(LPptr lp, int numCols, double sigma) {
	int    n;
	vector qsepvec;

	if (!(qsepvec = arr_alloc(numCols+1, double)))
		errMsg("Allocation", "constructQP", "qsepvec",0);

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (n = 0; n < numCols; n++)
		qsepvec[n] = 0.5 * sigma;
	qsepvec[n] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	if ( copyQPseparable(lp, qsepvec) ) {
		errMsg("solver", "constructQP", "failed to copy Q matrix", 0);
		return 1;
	}

	mem_free(qsepvec);
	return 0;
}//END changeQPproximal()

/* In the regularized QP method, we need to change the rhs of x to d. The
 * 		 A * x 			= b
 * 		 eta + beta * x >= alpha
 * Since x = xbar + d, the corresponding changes will therefore be:
 * 		 A * d = b - A * xbar
 * 		 eta + beta * d >= alpha - beta * xbar
 * But as long as the incumbent sulotion does not change, b - A * xbar and alpha - beta * xbar (for the existing cuts) won't change. So we only need
 * to change it when the incumbent changes.
 *
 * On the other hand, in each iteration, a new cut will be added (and/or some cuts may be dropped) and therefore we need to shift the rhs of the
 * added cut from _alpha_ to _alpha - beta * xbar_, which has taken care of in the routine addCut() in cuts.c. We do not need to worry about the shift
 * of rhs for the dropped cuts.
 * This function performs the change of rhs when the incumbent changes, as described above. */
int changeQPrhs(probType *prob, cellType *cell, vector xk) {
	int 	status = 0, cnt;
	vector 	rhs;
	intvec 	indices;

	if (!(rhs =(vector) arr_alloc(prob->num->rows+cell->cuts->cnt+1, double)))
		errMsg("Allocation", "changeRhs", "rhs",0);
	if (!(indices =(intvec) arr_alloc(prob->num->rows+cell->cuts->cnt, int)))
		errMsg("Allocation", "changeRhs", "indices",0);
	/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
	 vector, and the x vector. */
	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, xk, rhs);

	/*** new rhs = alpha - beta * xbar (benders cuts)***/
	for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
		rhs[prob->num->rows+cnt+1] = cell->cuts->vals[cnt]->alpha - vXv(cell->cuts->vals[cnt]->beta, xk, NULL, prob->sp->mac);
		indices[prob->num->rows+cnt] = cell->cuts->vals[cnt]->rowNum;

		cell->cuts->vals[cnt]->alphaIncumb = rhs[prob->num->rows+cnt+1];
	}

	/* Now we change the right-hand of the master problem. */
	status = changeRHS(cell->master->lp, prob->num->rows + cell->cuts->cnt, indices, rhs+1);
	if (status)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrhs()

/* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
 * (incumbent solution). */
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk) {
	int 	status = 0, cnt;
	vector	lbounds, ubounds;
	intvec	lindices, uindices;
	char 	*llu, *ulu;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);
	if (!(llu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "llu",0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);
	if (!(ulu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "ulu",0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt] - xk[cnt + 1];
		uindices[cnt] = cnt;
		ulu[cnt] = 'U';
	}

	status = changeBDS(lp, numCols, uindices, ulu, ubounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the upper bound in the solver", 0);
		return 1;
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt] - xk[cnt + 1];
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
	}

	status = changeBDS(lp, numCols, lindices, llu, lbounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the lower bound in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices); mem_free(llu);
	mem_free(ubounds); mem_free(uindices); mem_free(ulu);

	return 0;
}//END changeQPbds()

/* This subroutine initializes the master problem by copying information from the decomposed prob[0](type: oneProblem) and adding a column for
 * theta for modified benders decomposition. */
/*Jiajun, check add L1?*/
oneProblem *newMaster(oneProblem *orig, double lb) {
	oneProblem 	*master;
	int         r, i, j, idx, cnt;
	//long        colOffset, rowOffset;
	//char        *q;

	if (!(master = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "new_master", "Faile to allocate memory to mcell->sp", 0);

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Allocating memory to master -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	master->type 	= config.MASTER_TYPE;               /* type of problem: LP, QP, MIP or MIQP . Jiajun L1 here? */
	master->objsen 	= orig->objsen;                 	/* sense of the objective: 1 for minimization and -1 for maximization */
	master->mar 	= orig->mar;                       	/* number of rows */
	master->numInt 	= orig->numInt;                 	/* number of integer variables in the problem  */
	master->numnz 	= orig->numnz;                   	/* number of non-zero elements in constraint matrix */
	master->matsz 	= orig->matsz;                   	/* extended matrix size */
	master->marsz 	= orig->marsz;                   	/* extended row size */
	master->rstorsz = orig->rstorsz;               		/* memory size for storing row names */
	master->mac 	= orig->mac+1;           			/* number of columns + etas */
	master->macsz 	= orig->macsz + 1;       			/* extended column size */
	master->cstorsz = orig->cstorsz + NAMESIZE;    		/* memory size for storing column names */

	/* Allocate memory to the information whose type is string */
	if (!(master->name = (string) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->name",0);
	if (!(master->senx = (string) arr_alloc(master->marsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->senx",0);
	if (!(master->ctype = (string) arr_alloc(master->macsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->ctype",0);
	if (!(master->objname = (string) arr_alloc(NAMESIZE,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objname",0);
	if (!(master->cname = (string*) arr_alloc(master->macsz,string)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cname",0);
	if (!(master->cstore = (string) arr_alloc(master->cstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cstore",0);
	if ( master->mar > 0 ) {
		if (!(master->rname = (string *) arr_alloc(master->marsz,string)))
			errMsg("Allocation", "new_master", "Fail to allocate memory to master->rname",0);
		if (!(master->rstore = (string) arr_alloc(master->rstorsz, char)))
			errMsg("Allocation", "new_master", "Fail to allocate memory to master->rstore",0);
	}
	else {
		master->rname = NULL; master->rstore = NULL;
	}

	/* Allocate memory to the information whose type is vector */
	if (!(master->objx = (vector) arr_alloc(master->macsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objx",0);
	if (!(master->rhsx = (vector) arr_alloc(master->marsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rhsx",0);
	if (!(master->matval = (vector) arr_alloc(master->matsz, double)))
		errMsg("allocation", "new_master", "master->matval",0);
	if (!(master->bdl = (vector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdl",0);
	if (!(master->bdu = (vector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdu",0);

	/* Allocate memory to the information whose type is intvec */
	if (!(master->matbeg = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matbeg",0);
	if (!(master->matcnt = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matcnt",0);
	if (!(master->matind = (intvec) arr_alloc(master->matsz, int)))
		errMsg("allocation", "new_master", "master->matind",0);

	strcpy(master->name, orig->name);           /* Copy problem name */
	strcpy(master->objname, orig->objname);     /* Copy objective name */

	/* Copy problem's column and row names, and calculate the pointers for master/copy row and column names. */
    for (i = 0; i < orig->cstorsz; i++){
        master->cstore[i] = orig->cstore[i];
    }
    //colOffset = master->cstore - orig->cstore;

    if ( master->mar > 0 ) {
        for (i = 0; i < orig->rstorsz; i++){
            master->rstore[i] = orig->rstore[i];
        }
        //rowOffset = master->rstore - orig->rstore;
    }

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < orig->mac; j++) {
		/* Copy objective function coefficients */
		master->objx[j] = orig->objx[j];
		/* Copy the decision variable type */
		master->ctype[j] = orig->ctype[j];
		/* Copy the upper bound and lower bound */
		master->bdu[j] = orig->bdu[j];
		master->bdl[j] = orig->bdl[j];
		/* Copy column names, offset by length */
        if (j == 0)
            master->cname[0] = master->cstore;
        else
            master->cname[j] = master->cname[j - 1] + (orig->cname[j] - orig->cname[j - 1]);
        //master->cname[j] = orig->cname[j] + colOffset;
		/* Copy the master sparse matrix beginning position of each column */
		master->matbeg[j] = cnt;
		/* Copy the sparse matrix non-zero element count */
		master->matcnt[j] = orig->matcnt[j];
		/* Loop through all non-zero elements in this column */
		for (idx = orig->matbeg[j]; idx < orig->matbeg[j] + orig->matcnt[j]; idx++) {
			/* Copy the non-zero coefficient */
			master->matval[cnt] = orig->matval[idx];
			/* Copy the row entry of the non-zero elements */
			master->matind[cnt] = orig->matind[idx];
			cnt++;
		}
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < orig->mar; r++) {
		/* Copy the right hand side value */
		master->rhsx[r] = orig->rhsx[r];
		/* Copy the constraint sense */
		master->senx[r] = orig->senx[r];
		/* Copy row names, offset by length */
        if (r == 0)
            master->rname[0] = master->rstore;
        else
            master->rname[r] = orig->rname[r - 1] + (orig->rname[r] - orig->rname[r - 1]);
        //master->rname[r] = orig->rname[r] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	//colOffset = orig->cstorsz;
	strcpy(master->cstore + orig->cstorsz, "eta");
	master->cname[orig->mac] = master->cstore + orig->cstorsz;
	master->objx[orig->mac] = 1.0;			// orig->mac is the last column in the original master
	master->ctype[orig->mac] = 'C';
	master->bdu[orig->mac] = INFBOUND;
	master->bdl[orig->mac] = lb;
	master->matbeg[orig->mac] = orig->numnz;	// Beginning point in matval/matind in eta columns. every eta column begins at the same address
	master->matcnt[orig->mac] = 0;               // Only optimality cuts has eta

//    printf("----------------\n");
//    for (i = 0; i <= orig->mar; i++){
//        printf("master->rhsx[i]=%f\n", master->rhsx[i]);fflush(NULL);
//        printf("master->rname[i]=%s\n", master->rname[i]);fflush(NULL);
//        printf("master->senx[i]=%c\n", master->senx[i]);fflush(NULL);
//    }
//    
//    for (i = 0; i <= orig->mac; i++){
//        printf("master->objx[i]=%f\n", master->objx[i]);fflush(NULL);
//        printf("master->lb[i]=%f\n", master->bdl[i]);fflush(NULL);
//        printf("master->ub[i]=%f\n", master->bdu[i]);fflush(NULL);
//        printf("master->cname[i]=%s\n", master->cname[i]);fflush(NULL);
//        printf("master->matbeg[i]=%d\n", master->matbeg[i]);fflush(NULL);
//        printf("master->matcnt[i]=%d\n\n", master->matcnt[i]);fflush(NULL);
//    }
    
	/* Load the copy into CPLEX, Jiajun master->type is input to setupProblem */
	master->lp = setupProblem(master->name, master->type, master->mac, master->mar, master->objsen, master->objx, master->rhsx, master->senx, master->matbeg, master->matcnt,master->matind, master->matval, master->bdl, master->bdu, NULL, master->cname, master->rname, master->ctype);
	if ( master->lp == NULL ) {
		errMsg("Problem Setup", "new_master", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(ALGO_CHECK)
	if ( writeProblem(master->lp, "newMaster.lp") ) {
		errMsg("solver", "newMaster", "failed to write master problem to file", 0);
		return NULL;
	}
#endif

	return master;

}//END newMaster()

