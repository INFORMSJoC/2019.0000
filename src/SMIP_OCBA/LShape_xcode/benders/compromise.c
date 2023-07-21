/*
 * compromise.c
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

int buildCompromise(probType *prob, cellType *cell, batchSummary *batch) {
	vector	coef;
	intvec	indices;
	int 	i, idx, cnt, cOffset, rOffset1=0, rOffset2=0;
	char 	*q, tempName[NAMESIZE], batchNameSuffix[NAMESIZE];

	sprintf(batchNameSuffix, "_B%02d", batch->cnt);

	batch->ck[batch->cnt] 	 = cell->k;
	//batch->objLB[batch->cnt] = cell->incumbEst;
    
    if (config.MASTER_TYPE == PROB_QP){
        batch->objLB[batch->cnt] = cell->incumbEst;
        batch->incumbX[batch->cnt] = duplicVector(cell->incumbX, prob->num->cols);
    }
    else{
        batch->objLB[batch->cnt] = cell->candidEst;
        batch->incumbX[batch->cnt] = duplicVector(cell->candidX, prob->num->cols);
    }
	batch->cnt++;

	/* a. Setup or update the batch problem */
	batch->sp->numInt   += prob->sp->numInt;
	batch->sp->matsz	+= prob->sp->matsz;
	batch->sp->marsz 	+= prob->sp->marsz;
	batch->sp->macsz	+= prob->sp->macsz+1;
	batch->sp->cstorsz  += prob->sp->cstorsz + strlen(batchNameSuffix) * batch->sp->macsz;
	batch->sp->rstorsz  += prob->sp->rstorsz + strlen(batchNameSuffix) * batch->sp->marsz;
    
#if defined(BATCH_CHECK)
    printf("batch->cnt=%d, batch->sp->mar = %d, batch->sp->mac = %d, batch->sp->macsz= %d, batch->sp->marsz= %d, batch->sp->matsz= %d\n", batch->cnt, batch->sp->mar, batch->sp->mac,batch->sp->macsz, batch->sp->marsz, batch->sp->matsz);
#endif

	batch->sp->rhsx 	= (vector) mem_realloc(batch->sp->rhsx, batch->sp->marsz*sizeof(double));
	batch->sp->senx		= (string) mem_realloc(batch->sp->senx, batch->sp->marsz*sizeof(char));
	batch->sp->rstore 	= (string) mem_realloc(batch->sp->rstore, (batch->sp->rstorsz+prob->sp->marsz*NAMESIZE)*sizeof(char));
	batch->sp->rname 	= (string *) mem_realloc(batch->sp->rname, batch->sp->marsz*sizeof(string));
	batch->sp->cstore 	= (string) mem_realloc(batch->sp->cstore, (batch->sp->cstorsz+prob->sp->macsz*NAMESIZE)*sizeof(char));
	batch->sp->cname 	= (string *) mem_realloc(batch->sp->cname, batch->sp->macsz*sizeof(string));

	batch->sp->objx 	= (vector) mem_realloc(batch->sp->objx, batch->sp->macsz*sizeof(double));
	batch->sp->ctype 	= (string) mem_realloc(batch->sp->ctype, batch->sp->macsz*sizeof(char));
	batch->sp->matval 	= (vector) mem_realloc(batch->sp->matval, batch->sp->matsz*sizeof(double));
	batch->sp->bdl 		= (vector) mem_realloc(batch->sp->bdl, batch->sp->macsz*sizeof(double));
	batch->sp->bdu 		= (vector) mem_realloc(batch->sp->bdu, batch->sp->macsz*sizeof(double));

	batch->sp->matbeg 	= (intvec) mem_realloc(batch->sp->matbeg, batch->sp->macsz*sizeof(int));
	batch->sp->matcnt 	= (intvec) mem_realloc(batch->sp->matcnt, batch->sp->macsz*sizeof(int));
	batch->sp->matind 	= (intvec) mem_realloc(batch->sp->matind, batch->sp->matsz*sizeof(int));

	/* Copy/append problem's column names */
	cnt = batch->sp->cstorsz; cOffset = idx = batch->sp->mac;
	batch->sp->cname[idx] = batch->sp->cstore + cnt;
	for (q = prob->sp->cname[0]; q < prob->sp->cname[0] + prob->sp->cstorsz; q++) {
		if ( *q == '\0' ) {
			for ( i = 0; i < strlen(batchNameSuffix); i++ )
				batch->sp->cstore[cnt++] = batchNameSuffix[i];
			if ( idx < (batch->sp->macsz-1) )
				batch->sp->cname[idx+1] = batch->sp->cstore + cnt + 1;
			idx++;
		}
		batch->sp->cstore[cnt++] = *q;
	}
	batch->sp->cstorsz = cnt; batch->sp->mac = idx;

	/* Copy/append problems's row names */
	if ( batch->sp->rstorsz > 0 ) {
		cnt = batch->sp->rstorsz;
		rOffset1 = batch->sp->mar; rOffset2 = idx = (batch->cnt-1)*prob->num->rows;
#if defined(BATCH_CHECK)
        printf("rOffset1 = %d, rOffset2 = %d\n", rOffset1, rOffset2);
#endif
		batch->sp->rname[idx] = batch->sp->rstore + cnt;
		for (q = prob->sp->rname[0]; q < prob->sp->rname[0] + prob->sp->rstorsz; q++) {
			if ( *q == '\0' ) {
				for ( i = 0; i < strlen(batchNameSuffix); i++ )
					batch->sp->rstore[cnt++] = batchNameSuffix[i];
				if ( idx < (batch->sp->marsz-1) )
					batch->sp->rname[idx+1] = batch->sp->rstore + cnt + 1;
				idx++;
			}
			batch->sp->rstore[cnt++] = *q;
		}
		batch->sp->rstorsz = cnt; batch->sp->mar += (idx - rOffset2);
	}
#if defined(BATCH_CHECK)
    printf("batch->sp->mar = %d, batch->sp->mac = %d, batch->sp->macsz= %d, batch->sp->marsz= %d, batch->sp->matsz= %d\n", batch->sp->mar, batch->sp->mac,batch->sp->macsz, batch->sp->marsz, batch->sp->matsz);
#endif

	/* Copy all column information from the batch master problem */
	cnt = batch->sp->numnz;
	for (i = 0; i < prob->sp->mac; i++) {
		batch->sp->objx[cOffset+i] 	 = prob->sp->objx[i];		/* Copy objective function coefficients */
		batch->sp->ctype[cOffset+i]  = prob->sp->ctype[i];		/* Copy the decision variable type */
		batch->sp->bdu[cOffset+i] 	 = prob->sp->bdu[i];		/* Copy the upper bound and lower bound */
		batch->sp->bdl[cOffset+i] 	 = prob->sp->bdl[i];
		batch->sp->matbeg[cOffset+i] = cnt;						/* Copy the master sparse matrix beginning position of each column */
		batch->sp->matcnt[cOffset+i] = prob->sp->matcnt[i];		/* Copy the sparse matrix non-zero element count */
		batch->sp->ctype[cOffset+i]  = prob->sp->ctype[i];
		/* Loop through all non-zero elements in this column */
		for (idx = prob->sp->matbeg[i]; idx < prob->sp->matbeg[i] + prob->sp->matcnt[i]; idx++) {
			batch->sp->matval[cnt] = prob->sp->matval[idx];				/* Copy the non-zero coefficient */
			batch->sp->matind[cnt] = rOffset1+prob->sp->matind[idx];		/* Copy the row entry of the non-zero elements */
			cnt++;
		}
	}
	batch->sp->numnz = cnt;

	/* Copy all information concerning rows of batch master problem */
	for (i = 0; i < prob->sp->mar; i++) {
		batch->sp->rhsx[rOffset2+i]  = prob->sp->rhsx[i];		/* Copy the right hand side value */
		batch->sp->senx[rOffset2+i]  = prob->sp->senx[i];		/* Copy the constraint sense */
	}

	/* Initialize information for the extra column in the new master. */
	idx = batch->sp->mac++;
	sprintf(tempName, "eta%s", batchNameSuffix);
	strcpy(batch->sp->cstore+batch->sp->cstorsz, tempName);
	batch->sp->cname[idx] 	= batch->sp->cstore + batch->sp->cstorsz;
	batch->sp->objx[idx] 	= 1.0;
	batch->sp->ctype[idx] 	= 'C';
	batch->sp->bdu[idx]		= INFBOUND;
	batch->sp->bdl[idx] 	= prob->lb;
	batch->sp->matbeg[idx] 	= prob->sp->numnz;
	batch->sp->matcnt[idx] 	= 0;

	/* If this is the first time then we need to load the batch problem solver, else append the new information to the already
	 * existent batch problem. */
	if ( batch->cnt == 1 ) {
		/* Load the copy into CPLEX */
		batch->sp->lp = setupProblem(batch->sp->name, batch->sp->type, batch->sp->mac, batch->sp->mar, batch->sp->objsen,
				batch->sp->objx, batch->sp->rhsx, batch->sp->senx, batch->sp->matbeg, batch->sp->matcnt,batch->sp->matind,
				batch->sp->matval, batch->sp->bdl, batch->sp->bdu, NULL, batch->sp->cname, batch->sp->rname, batch->sp->ctype);
		if ( batch->sp->lp == NULL ) {
			errMsg("Batch problem setup", "buildCompromise", "failed to setup master problem in the solver",0);
			return 1;
		}
	}
	else {
		/* Add new rows to the problem */
		for  ( i = 0; i < prob->num->rows; i++ ) {
			sprintf(tempName, "%s%s", prob->sp->rname[i], batchNameSuffix);
			if ( addRow(batch->sp->lp, 0, batch->sp->rhsx[i+rOffset2], batch->sp->senx[i+rOffset2],
					0, NULL, NULL, tempName) ) {
				errMsg("solver", "buildCompromise", "failed to add new row to problem in solver", 0);
				return 1;
			}
		}

		/* Add new columns to the problem */
		for ( i = 0; i < prob->num->cols; i++ ) {
			if ( addCol(batch->sp->lp, batch->sp->matcnt[i+cOffset], batch->sp->objx[i+cOffset],
					0, batch->sp->matind+batch->sp->matbeg[i+cOffset], batch->sp->matval+batch->sp->matbeg[i+cOffset],
					batch->sp->bdu[i+cOffset], batch->sp->bdl[i+cOffset], batch->sp->cname[i+cOffset])) {
				errMsg("Batch problem setup", "buildCompromise", "failed to setup master problem in the solver",0);
				return 1;
			}
		}
		if ( addCol(batch->sp->lp, batch->sp->matcnt[i+cOffset], batch->sp->objx[i+cOffset],
				0, batch->sp->matind, batch->sp->matval,
				batch->sp->bdu[i+cOffset], batch->sp->bdl[i+cOffset], batch->sp->cname[i+cOffset])) {
			errMsg("Batch problem setup", "buildCompromise", "failed to setup master problem in the solver",0);
			return 1;
		}
	}

	/* b. Change the right-hand side with the incumbent solution of the current batch */
	indices = (intvec) arr_alloc(max(prob->num->cols,prob->num->rows)+1, int);
	coef = (vector) arr_alloc(prob->num->rows+1, double);
	for (i = 0; i < prob->num->rows; i++) {
		coef[i+1]  = prob->sp->rhsx[i];
		indices[i] = i+rOffset1;
	}

//    /* b - A * xbar */
//    coef = MSparsexvSub(prob->Dbar, cell->incumbX, coef);

	/* change the right-hand of the problem on the solver. */
	if ( changeRHS(batch->sp->lp, prob->num->rows, indices, coef+1) ) {
		errMsg("solver", "buildCompromise", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	/* c. Add the cuts in the problem */
	for ( i = 0; i < max(prob->num->cols,prob->num->rows); i++ )
		indices[i+1] = i+cOffset;
	indices[0] = idx;

	/* Optimality minorants */
	batch->sp->mar += cell->cuts->cnt;
	for ( cnt = 0; cnt < cell->cuts->cnt; cnt++ ) {
		sprintf(tempName, "%s%s", cell->cuts->vals[cnt]->name, batchNameSuffix);
        //TODO:cell->cuts->vals[cnt]->beta[0] = (double) (batch->ck[batch->cnt-1]) / (double) cell->cuts->vals[cnt]->numSamples;
		if ( addRow(batch->sp->lp, prob->num->cols+1, cell->cuts->vals[cnt]->alphaIncumb, GE, 0, indices,
				cell->cuts->vals[cnt]->beta, tempName) ) {
			errMsg("solver", "buildCompromise", "failed to add new cut row to problem in solver", 0);
			return 1;
		}
	}
    
    /*record the time information from replications*/
    
    batch->time->repTime += cell->time->repTime;
    batch->time->iterTime += cell->time->iterTime;
    batch->time->masterIter += cell->time->masterIter ;
    batch->time->subprobIter += cell->time->subprobIter;
    batch->time->optTestIter += cell->time->optTestIter;
    batch->time->iterAccumTime += cell->time->iterAccumTime;
    batch->time->masterAccumTime += cell->time->masterAccumTime;
    batch->time->subprobAccumTime += cell->time->subprobAccumTime;
    batch->time->optTestAccumTime += cell->time->optTestAccumTime;

	/* Feasibility cuts, if any */
//    batch->sp->mar += cell->fcuts->cnt;
//    for ( cnt = 0; cnt < cell->fcuts->cnt; cnt++ ) {
//        sprintf(tempName, "%s%s", cell->fcuts->vals[cnt]->name, batchNameSuffix);
//        if ( addRow(batch->sp->lp, prob->num->cols+1, cell->fcuts->vals[cnt]->alphaIncumb, GE, 0, indices,
//                cell->fcuts->vals[cnt]->beta, tempName) ) {
//            errMsg("solver", "buildCompromise", "failed to add new feasibility cut row to problem in solver", 0);
//            return 1;
//        }
//    }

	/* d. Change the bounds on the batch variables */
//    if (changeQPbds(batch->sp->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, cell->incumbX, cOffset) ) {
//        errMsg("algorithm", "buildCompromise", "failed to change the bounds to convert the problem into QP", 0);
//        return 1;
//    }

	/* e. Change the proximal term */
    //Jiajun BUG:
	//for ( i = 0; i < batch->cnt; i++ )
    //batch->quadScalar = ((batch->cnt-1)*batch->quadScalar + cell->quadScalar)/batch->cnt;
    batch->quadScalar = config.MIN_QUAD_SCALAR;

//    qsepvec = (vector) arr_alloc(batch->sp->mac, double); idx = 0;
//    for ( cnt = 0; cnt < batch->cnt; cnt++ ) {
//        for (i = 0; i < prob->num->cols; i++)
//            qsepvec[idx++] = 0.5 * batch->quadScalar;
//        qsepvec[idx++] = 0.0;
//    }
//
//    if ( copyQPseparable(batch->sp->lp, qsepvec) ) {
//        errMsg("solver", "changeQPproximal", "failed to copy Q matrix", 0);
//        return 1;
//    }

#if defined(BATCH_CHECK)
	sprintf(tempName, "master_B%02d.lp", batch->cnt);
	if ( writeProblem(cell->master->lp, tempName) ) {
		errMsg("solver", "buildCompromise", "failed to write master problem to file", 0);
		return 1;
	}
	if ( writeProblem(batch->sp->lp, "batch.lp") ) {
		errMsg("solver", "buildCompromise", "failed to write master problem to file", 0);
		return 1;
	}
#endif

	mem_free(indices); mem_free(coef);
    //mem_free(qsepvec);

	return 0;
}//END buildCompromise()

int solveCompromise(probType *prob, batchSummary *batch) {
	int b, j, status, tempX;
    clock_t    tic;
    
	/* Add the equality constraints that tie the replications together */
	if ( addBatchEquality(prob, batch) ) {
		errMsg("algorithm", "solveCompromise", "failed to add the equality constraint", 0);
		return 1;
	}
    
    if ( addL1Norm(prob, batch) ) {
        errMsg("algorithm", "solveCompromise", "failed to add L1 norm", 0);
        return 1;
    }
    


	/* solve the compromise problem */
    if(config.MASTER_TYPE == PROB_QP){
        changeQPSolverType(ALG_CONCURRENT);
    }
    else if(config.MASTER_TYPE == PROB_MILP){
        changeMILPSolverType(ALG_CONCURRENT);
    }
	
    tic = clock();
	if ( solveProblem(batch->sp->lp, batch->sp->name, config.MASTER_TYPE, &status) ) {
		writeProblem(batch->sp->lp, "error.lp");
		errMsg("algorithm", "solveCompromise", "failed to solve the compromise problem", 0);
		return 1;
	}
    batch->time->masterIter += ((double) (clock() - tic))/CLOCKS_PER_SEC;
    
    if(config.MASTER_TYPE == PROB_MILP){
        batch->Est = getObjective(batch->sp->lp, PROB_MILP);
    }
    else{
        batch->Est = getObjective(batch->sp->lp, PROB_LP);
    }
    batch->Est = batch->Est / config.NUM_REPS;

	/* Get the primal solution to the compromise problem */
	batch->compromiseX = (vector) arr_alloc(prob->num->cols+1, double);
//    if(config.MASTER_TYPE == PROB_MILP)
//        getMIPPrimal(batch->sp->lp, batch->compromiseX, prob->num->cols);
//    else
    getPrimal(batch->sp->lp, batch->compromiseX, prob->num->cols);
    
//    for ( j = 1; j <= prob->num->cols; j++ )
//        batch->compromiseX[j] += batch->incumbX[0][j];

	/* Evaluate the average solution */
	batch->avgX = (vector) arr_alloc(prob->num->cols+1, double);
	for (j = 1; j <= prob->num->cols; j++) {
		batch->avgX[j] = 0.0;
		for ( b = 0; b < batch->cnt; b++) {
			batch->avgX[j] += batch->incumbX[b][j];
		}
		batch->avgX[j] /= (double) batch->cnt;
        tempX = (int) batch->avgX[j];
        if ((batch->avgX[j] - tempX) >= 0.5)
            batch->avgX[j] = tempX + 1;
        else
            batch->avgX[j] = tempX;
	}

	return 0;
}//END solveCompromise()

int addBatchEquality (probType *prob, batchSummary *batch) {
	double coef[2] = {1.0, -1.0}, rhs;
	int	j, b, indices[2];
	char tempName[NAMESIZE];

	for ( j = 0; j < prob->num->cols; j++ ) {
		indices[0] = j;
		for ( b = 1; b < batch->cnt; b++ ) {
			sprintf(tempName, "%s_B%02d_%02d", prob->sp->cname[j], 1, b);
			indices[1] = b*(prob->num->cols+1) + j;
			//rhs = batch->incumbX[b][j+1] - batch->incumbX[0][j+1];
            rhs = 0;
			if ( addRow(batch->sp->lp, 2, rhs, EQ, 0, indices, coef, tempName) ) {
				errMsg("solver", "buildCompromise", "failed to add new feasibility cut row to problem in solver", 0);
				return 1;
			}
		}
	}

#if defined(BATCH_CHECK)
	if ( writeProblem(batch->sp->lp, "batchwEqConstraints.lp") ) {
		errMsg("solver", "buildCompromise", "failed to write master problem to file", 0);
		return 1;
	}
#endif

	return 0;
}//END addBatchEquality()

batchSummary *newBatchSummary(probType *prob, int numBatch) {
	batchSummary *batch;

	/* Setup batch summary structure */
	batch = (batchSummary *) mem_malloc(sizeof(batchSummary));
	batch->ck = (intvec) arr_alloc(numBatch, int);
	batch->objLB = (vector) arr_alloc(numBatch, double);
	batch->objUB = (vector) arr_alloc(numBatch, double);
	batch->incumbX = (vector *) arr_alloc(numBatch, vector);
	batch->cnt = 0;
	batch->quadScalar = config.MIN_QUAD_SCALAR;
	batch->avgX = batch->compromiseX = NULL;
    batch->Est = 0;

	/* Setup the elements of the batch problem */
	batch->sp = (oneProblem *) mem_malloc(sizeof(oneProblem));
	batch->sp->lp = NULL;
	batch->sp->type 	= config.MASTER_TYPE;
	batch->sp->objsen 	= prob->sp->objsen;

	batch->sp->mar = batch->sp->numInt = batch->sp->numnz =  batch->sp->mac = 0;
	batch->sp->matsz = batch->sp->marsz = batch->sp->rstorsz = batch->sp->macsz = batch->sp->cstorsz = 0;

	batch->sp->name 	= (string) arr_alloc(NAMESIZE, char);	strcpy(batch->sp->name, "batchProblem");
	batch->sp->objname 	= (string) arr_alloc(NAMESIZE, char);	strcpy(batch->sp->objname, prob->sp->objname);

	batch->sp->rhsx 	= (vector) arr_alloc(1, double);
	batch->sp->senx		= (string) arr_alloc(1, char);
	batch->sp->rname 	= (string *) arr_alloc(1, string);
	batch->sp->rstore 	= (string) arr_alloc(1, char);
	batch->sp->cname 	= (string *) arr_alloc(1, string);
	batch->sp->cstore 	= (string) arr_alloc(1, char);

	batch->sp->objx 	= (vector) arr_alloc(1, double);
	batch->sp->ctype 	= (string) arr_alloc(1, char);
	batch->sp->matval 	= (vector) arr_alloc(1, double);
	batch->sp->bdl 		= (vector) arr_alloc(1, double);
	batch->sp->bdu 		= (vector) arr_alloc(1, double);

	batch->sp->matbeg 	= (intvec) arr_alloc(1, int);
	batch->sp->matcnt 	= (intvec) arr_alloc(1, int);
	batch->sp->matind 	= (intvec) arr_alloc(1, int);
    
    if ( !(batch->time = (runTime *) mem_malloc(sizeof(runTime)) ) ){
        errMsg("setup", "newCell", "cell->runTime", 0);}
    batch->time->repTime = batch->time->iterTime = batch->time->masterIter = batch->time->subprobIter = batch->time->optTestIter = 0.0;
    batch->time->iterAccumTime = batch->time->masterAccumTime = batch->time->subprobAccumTime = batch->time->optTestAccumTime = 0.0;

	return batch;
}//END newBatchSummary()

void freeBatchType(batchSummary *batch) {
	int n;

	if ( batch ) {
		if (batch->ck) mem_free(batch->ck);
		if (batch->objLB) mem_free(batch->objLB);
		if (batch->objUB) mem_free(batch->objUB);
		if (batch->incumbX) {
			for ( n = 0; n < batch->cnt; n++ ) {
				if (batch->incumbX[n]) mem_free(batch->incumbX[n]);
			}
			mem_free(batch->incumbX);
		}
		if (batch->compromiseX) mem_free(batch->compromiseX);
		if (batch->avgX) mem_free(batch->avgX);
		if ( batch->sp ) freeOneProblem(batch->sp);
        if (batch->time) mem_free(batch->time);
		mem_free(batch);
	}

}

int addL1Norm (probType *prob, batchSummary *batch) {
    int i, j, cnt;
    int status;
    int NumABS;
    intvec rmatbeg, rmatind;
    vector obj, rhs, rmatval, lb, ub;
    string sense;
    
    NumABS = prob->sp->macsz * batch->cnt;
    
    obj = (vector) arr_alloc(NumABS, double);
    lb = (vector) arr_alloc(NumABS, double);
    ub = (vector) arr_alloc(NumABS, double);
    
    for (i = 0; i < NumABS; i++){
        obj[i] = batch->quadScalar;
        lb[i] = 0;
        ub[i] = 1;
    }
    /*add new variables: d_i = abs(X_i - incumbX_i)*/
    status = Newcols (batch->sp->lp, NumABS, obj, lb, ub, NULL, NULL);
    if(status)
        errMsg("addL1Norm", "Newcols", "failed to add new cols", 0);
    
    
#if defined(BATCH_CHECK)
    if ( writeProblem(batch->sp->lp, "batchwNorm_addcols.lp") ) {
        errMsg("solver", "buildCompromise", "failed to write master problem to file", 0);
        return 1;
    }
#endif
    
    /*add the constraints: d_i >= X_i - incumbX_i, d_i >= -X_i +incumbX_i, which are
     -x_i +d_i >= -incumbX_i,
     x_i +d_i >= incumbX_i,
     Thus,
     rhs = [-incumbX_i, incumbX_i ...]
     rmatbeg = [0, 2, 4, ...]
     rmatind = [0, col_number * batch_cnt, 0, col_number * batch_cnt, 1, col_number * batch_cnt+1, 1, col_number * batch_cnt+1...], col_number =first stage variable number + 1
     rmatval = [-1, 1, 1, 1] * NumABS; */
    rhs = (vector) arr_alloc(2*NumABS, double);
    sense = (string) arr_alloc(2*NumABS, char);
    
    rmatbeg = (intvec) arr_alloc(2*NumABS, int);
    rmatind = (intvec) arr_alloc(4*NumABS, int);
    rmatval = (vector) arr_alloc(4*NumABS, double);
    
    for (i = 0; i < batch->cnt; i++){
        for (j = 0; j < prob->sp->macsz; j++){
            cnt = i * prob->sp->macsz + j;
            rhs[2 * cnt] = -batch->incumbX[i][j+1]; rhs[2 * cnt+ 1] = batch->incumbX[i][j+1];
            sense[2 * cnt] = 'G'; sense[2 * cnt + 1] = 'G';
            rmatbeg[2 * cnt] = 4 * cnt; rmatbeg[2 * cnt + 1] = 4 * cnt + 2;
            
            rmatind[4 * cnt] = rmatind[4 * cnt + 2] = i * (prob->sp->macsz + 1) + j;
            rmatind[4 * cnt + 1] = rmatind[4 * cnt + 3] = batch->cnt * (prob->sp->macsz + 1) + cnt;
            
            rmatval[4 * cnt] = -1;
            rmatval[4 * cnt + 1] = rmatval[4 * cnt + 2] = rmatval[4 * cnt + 3] = 1;
            
        }
    }

#ifdef BATCH_CHECK
    printf("\n rhs:\n");
    for(i=0;i < 2*NumABS;i++)
        printf("%f ", rhs[i]);
    printf("\n rmatbeg:\n");
    for(i=0;i < 2*NumABS;i++)
        printf("%d ", rmatbeg[i]);
    printf("\n rmatind:\n");
    for(i=0;i < 4*NumABS;i++)
        printf("%d ", rmatind[i]);
    printf("\n rmatval:\n");
    for(i=0;i < 4*NumABS;i++)
        printf("%f ", rmatval[i]);
#endif
    
    
    status = Addrows (batch->sp->lp, 0, 2*NumABS, 4*NumABS, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
    if(status)
        errMsg("addL1Norm", "Addrows", "failed to add rows", 0);
    
    #if defined(BATCH_CHECK)
        if ( writeProblem(batch->sp->lp, "batchwNorm.lp") ) {
            errMsg("solver", "buildCompromise", "failed to write master problem to file", 0);
            return 1;
        }
    #endif

    mem_free(obj);
    mem_free(rhs);
    mem_free(rmatval);
    mem_free(lb);
    mem_free(ub);
    mem_free(sense);
    mem_free(rmatbeg); mem_free(rmatind);
    
    return 0;
}

BOOL InConvexHull(batchSummary *batch, int col){
    int i, j;
    BOOL flag=TRUE;
    int counter=0;
    
    for(i = 0; i < batch->cnt; i++){
        flag = TRUE;
        for(j = 0; j < col; j++)
            if(DBL_ABS(batch->compromiseX[j] - batch->incumbX[i][j]) > config.TOLERANCE){
                flag = FALSE;
                break;
            }
        if(flag) counter++;
    }
    if (counter > 14){
        return TRUE;
    }
    return FALSE;
}

void LowerBoundVariance(batchSummary *batch, double *mean, double *std){
    int cnt;
    double temp, variance=0;
    
    *mean = batch->objLB[0] ;
    for(cnt = 1; cnt < batch->cnt; cnt++){
        temp = *mean;
        *mean = *mean + (batch->objLB[cnt] - *mean) / (double) (cnt + 1.0);
        variance  = cnt / (cnt + 1.0) * variance
        + cnt * (*mean - temp) * (*mean - temp);
    }
    *std = sqrt(variance/ (double) cnt);
}
