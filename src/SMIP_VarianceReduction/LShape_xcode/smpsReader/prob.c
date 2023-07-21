/*
 * prob.c
 *
 *  Created on: Feb 17, 2019
 *      Author: Jiajun Xu
 * Institution: University of Southern California
 *
 * Please send your comments or bug report to jiajunx (at) usc (dot) edu
 */

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

/* Decomposes the problem _orig_ into subproblems as well as decomposes the stochastic information _stoc_ into stage stochastic information. The decomposition
 * is carried out using information specified in _tim_. The function also stores stage lower bound information provided in _Lb_. It returns an array of
 * probType structures, each probType corresponds to a particular stage */
probType **newProb(oneProblem *orig, stocType *stoc, timeType *tim, vector lb, double TOLERANCE) {
	probType **prob;
	char	 *q;
    int		 i, k, m, t;
    int rOffset = 0, cOffset = 0;

	/* allocate memory to elements of probType */
	if ( !(prob = (probType **) arr_alloc(tim->numStages, probType *)) )
		errMsg("allocation", "newProb", "prob", 0);

	/* allocate memory to members of probType for stagewise subProblems, and allocate values to static fields*/
	for ( t = 0; t < tim->numStages; t++ ) {
		if ( !(prob[t] = (probType *) mem_malloc(sizeof(probType))) )
			errMsg("allocation", "newProb", "prob[t]", 0);
		if ( !(prob[t]->name = (string) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "newProb", "stage name", 0);
		if( !(prob[t]->sp = (oneProblem *) mem_malloc (sizeof(oneProblem))))
			errMsg("allocation", "newProb", "stage problem", 0);
		prob[t]->sp->lp = NULL; prob[t]->lb = 0.0;

		if (  !(prob[t]->dBar = (sparseVector *) mem_malloc(sizeof(sparseVector))) )
			errMsg("allocation", "newProb", "stage problem cost coefficients", 0);
		if ( !(prob[t]->dBar->col = (intvec) arr_alloc(orig->mac+1, int)) )
			errMsg("allocation", "newProb", "stage problem cost coefficients columns", 0);
		if ( !(prob[t]->dBar->val = (vector) arr_alloc(orig->mac+1, double)) )
			errMsg("allocation", "newProb", "stage problem cost coefficients values", 0);
		prob[t]->dBar->cnt = 0;

		if ( !(prob[t]->bBar = (sparseVector *) mem_malloc(sizeof(sparseVector))) )
			errMsg("allocation", "newProb", "stage problem right hand side", 0);
		if ( !(prob[t]->bBar->col = (intvec) arr_alloc(orig->mar+1, int)) )
			errMsg("allocation", "newProb", "stage problem right-hand rows", 0);
		if ( !(prob[t]->bBar->val = (vector) arr_alloc(orig->mar+1, double)) )
			errMsg("allocation", "newProb", "stage problem right-hand values", 0);
		prob[t]->bBar->cnt = 0;

		strcpy(prob[t]->name, tim->stgNames[t]);

		if ( t < tim->numStages - 1 ) {
			prob[t]->sp->mar = prob[t]->sp->marsz = tim->row[t+1] - tim->row[t];
			prob[t]->sp->mac = prob[t]->sp->macsz = tim->col[t+1] - tim->col[t];
			prob[t]->sp->rstorsz = (int) (orig->rname[tim->row[t+1]] - orig->rname[tim->row[t]]);
			prob[t]->sp->cstorsz = (int) (orig->cname[tim->col[t+1]] - orig->cname[tim->col[t]]);
			rOffset += prob[t]->sp->rstorsz;
			cOffset += prob[t]->sp->cstorsz;
            prob[t]->sp->type = PROB_MILP; //Jiajun, first stage MIP
		}
		else {
			prob[t]->sp->mar = prob[t]->sp->marsz = orig->mar - tim->row[t];
			prob[t]->sp->mac = prob[t]->sp->macsz = orig->mac - tim->col[t];
			prob[t]->sp->rstorsz = orig->rstorsz - rOffset;
			prob[t]->sp->cstorsz = orig->cstorsz - cOffset;
            prob[t]->sp->type = PROB_LP; //Jiajun, second stage LP
		}
		prob[t]->sp->numInt = 0;
		prob[t]->sp->numBin = 0;
		prob[t]->sp->matsz = 0;
		prob[t]->sp->numnz = 0;
		prob[t]->sp->objsen = orig->objsen;
		

		/* stage oneProblem */
		if(!(prob[t]->sp->name = (string) arr_alloc(NAMESIZE, char)))
			errMsg("allocation", "newProb", "stage problem name", 0);
		if(!(prob[t]->sp->objname = (string) arr_alloc(NAMESIZE, char)))
			errMsg("allocation", "newProb", "stage problem objname", 0);
		if(!(prob[t]->sp->objx = (vector) arr_alloc(prob[t]->sp->macsz, double)))
			errMsg("allocation", "newProb", "stage problem objx", 0);
		if(!(prob[t]->sp->bdl = (vector) arr_alloc(prob[t]->sp->macsz, double)))
			errMsg("allocation", "newProb", "stage problem bdl", 0);
		if(!(prob[t]->sp->bdu = (vector) arr_alloc(prob[t]->sp->macsz, double)))
			errMsg("allocation", "newProb", "stage problem bdu", 0);
		if(!(prob[t]->sp->ctype = (string) arr_alloc(prob[t]->sp->macsz, char)))
			errMsg("allocation", "newProb", "stage problem column type", 0);
		if(!(prob[t]->sp->rhsx = (vector) arr_alloc(prob[t]->sp->marsz, double)))
			errMsg("allocation", "newProb", "stage problem rhsx", 0);
		if(!(prob[t]->sp->senx = (string) arr_alloc(prob[t]->sp->marsz, char)))
			errMsg("allocation", "newProb", "stage problem senx", 0);
		if(!(prob[t]->sp->matbeg = (intvec) arr_alloc(prob[t]->sp->macsz, int)))
			errMsg("allocation", "newProb", "stage problem matbeg", 0);
		if(!(prob[t]->sp->matcnt = (intvec) arr_alloc(prob[t]->sp->macsz, int)))
			errMsg("allocation", "newProb", "stage problem matcnt", 0);
		if(!(prob[t]->sp->cname = (string *) arr_alloc(prob[t]->sp->macsz, string)))
			errMsg("allocation", "newProb", "stage problem cname", 0);
		if(!(prob[t]->sp->cstore = (string) arr_alloc(prob[t]->sp->cstorsz, char)))
			errMsg("allocation", "newProb", "stage problem cstore", 0);
		if(!(prob[t]->sp->rname = (string *) arr_alloc(prob[t]->sp->marsz, string)))
			errMsg("allocation", "newProb", "stage problem rname", 0);
		if(!(prob[t]->sp->rstore = (string) arr_alloc(prob[t]->sp->rstorsz, char)))
			errMsg("allocation", "newProb", "stage problem rstore", 0);
		if(!(prob[t]->sp->matval = (vector) arr_alloc(orig->matsz, double)))
			errMsg("allocation", "newProb", "stage problem matval", 0);
		if(!(prob[t]->sp->matind = (intvec) arr_alloc(orig->matsz, int)))
			errMsg("allocation", "newProb", "stage problem matind", 0);
		strcpy(prob[t]->sp->objname, orig->objname);
		sprintf(prob[t]->sp->name, "%s_%d", orig->name, t);

		/* stage transfer matrix */
		if ( t == 0 )
			prob[t]->Cbar = NULL;
		else {
			if ( !(prob[t]->Cbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
				errMsg("allocation", "newProb", "stage transfer matrix", 0);
			if(!(prob[t]->Cbar->row = (intvec) arr_alloc(orig->matsz + 1, int)))
				errMsg("allocation", "newProb", "transfer matrix rows", 0);
			if(!(prob[t]->Cbar->col = (intvec) arr_alloc(orig->matsz + 1, int)))
				errMsg("allocation", "newProb", "transfer matrix columns", 0);
			if(!(prob[t]->Cbar->val = (vector) arr_alloc(orig->matsz + 1, double)))
				errMsg("allocation", "newProb", "transfer matrix values", 0);
			prob[t]->Cbar->cnt = 0;
		}

		/* TODO (HG): stage dynamics: include dynamics to the model */
		if ( t == 0) {
			prob[t]->Abar = NULL;
			prob[t]->Bbar = NULL;
			prob[t]->aBar = NULL;
			prob[t]->cBar = NULL;
		}
		else {
			prob[t]->Abar = NULL;
			prob[t]->Bbar = NULL;
			prob[t]->aBar = NULL;
			prob[t]->cBar = NULL;
		}

		/* Stage recourse matrix: Dbar is used to setup the QP in the forward pass. Since the QP is used only for non-terminal
		 * stages, Dbar for terminal stage is set to NULL */
		if ( t == tim->numStages-1 )
			prob[t]->Dbar = NULL;
		else {
			if ( !(prob[t]->Dbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
				errMsg("allocation", "newProb", "stage constraint matrix", 0);
			if(!(prob[t]->Dbar->row = (intvec) arr_alloc(orig->matsz+1, int)))
				errMsg("allocation", "newProb", "stage constraint matrix rows", 0);
			if(!(prob[t]->Dbar->col = (intvec) arr_alloc(orig->matsz+1, int)))
				errMsg("allocation", "newProb", "stage constraint matrix columns", 0);
			if(!(prob[t]->Dbar->val = (vector) arr_alloc(orig->matsz+1, double)))
				errMsg("allocation", "newProb", "stage constraint matrix values", 0);
			prob[t]->Dbar->cnt = 0;
		}
	}

	for ( t = 0; t < tim->numStages-1; t++ ) {
		/* lower bound on cost-to-go function */
		if ( lb != NULL )
			prob[t]->lb = lb[t];

		/* copy column names of non-terminal stage*/
		m = 0;
		for ( q = orig->cname[tim->col[t]]; q < orig->cname[tim->col[t+1]]; q++ )
			prob[t]->sp->cstore[m++] = *q;

		/* copy column information for non-terminal stage */
		//cOffset = prob[t]->sp->cstore - orig->cname[tim->col[t]];
		for ( m = tim->col[t]; m < tim->col[t+1]; m++ ) {
			k = m - tim->col[t];
			prob[t]->dBar->val[prob[t]->dBar->cnt+1] = orig->objx[m];
			prob[t]->dBar->col[prob[t]->dBar->cnt+1] = m - tim->col[t]+1;
			prob[t]->dBar->cnt++;
			prob[t]->sp->objx[k] = orig->objx[m];
			prob[t]->sp->bdl[k] = orig->bdl[m];
			prob[t]->sp->bdu[k] = orig->bdu[m];
			if ( orig->ctype[m] == 'I')
				prob[t]->sp->numInt++;
			else if ( orig->ctype[m] == 'B' )
				prob[t]->sp->numBin++;
			prob[t]->sp->ctype[k] = orig->ctype[m];
			prob[t]->sp->cname[k] = orig->cname[m] + (prob[t]->sp->cstore - orig->cname[tim->col[t]]);
			prob[t]->sp->matcnt[k] = 0;
			for ( i = orig->matbeg[m]; i < orig->matbeg[m]+orig->matcnt[m]; i++ ) {
				if (orig->matind[i] < tim->row[t+1]) {
					/* The coefficient is part of the current stage constraint matrix */
					if ( k == 0 )
						prob[t]->sp->matbeg[k] = 0;
					else
						prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1] + prob[t]->sp->matcnt[k-1];
					prob[t]->sp->matval[prob[t]->sp->matsz] = orig->matval[i];
					prob[t]->sp->matind[prob[t]->sp->matsz] = orig->matind[i] - tim->row[t];
					++prob[t]->sp->matcnt[k];
					++prob[t]->sp->matsz;
					++prob[t]->sp->numnz;
					prob[t]->Dbar->val[prob[t]->Dbar->cnt+1] = orig->matval[i];
					prob[t]->Dbar->col[prob[t]->Dbar->cnt+1] = m-tim->col[t]+1;
					prob[t]->Dbar->row[prob[t]->Dbar->cnt+1] = orig->matind[i] - tim->row[t]+1;
					++prob[t]->Dbar->cnt;
				}
				else {
					/* The coefficient is part of the next stage transfer matrix */
					prob[t+1]->Cbar->val[prob[t+1]->Cbar->cnt+1] = orig->matval[i];
					prob[t+1]->Cbar->row[prob[t+1]->Cbar->cnt+1] = orig->matind[i] - tim->row[t+1]+1;
					prob[t+1]->Cbar->col[prob[t+1]->Cbar->cnt+1] = m-tim->col[t]+1;
					++prob[t+1]->Cbar->cnt;
				}
			}
		}

		/* if integer or binary variables are encountered, then label the stage problem as a mixed integer LP */
//        if ( prob[t]->sp->numInt + prob[t]->sp->numBin > 0 )
//            prob[t]->sp->type = PROB_MILP;

		/* copy row name of non-terminal stage */
		m = 0;
		for ( q = orig->rname[tim->row[t]]; q < orig->rname[tim->row[t+1]]; q++ )
			prob[t]->sp->rstore[m++] = *q;

		/* copy row information for non-terminal stage */
		//rOffset = prob[t]->sp->rstore - orig->rname[tim->row[t]];
		for ( m = tim->row[t]; m < tim->row[t+1]; m++ ) {
			k = m - tim->row[t];
			prob[t]->sp->rhsx[k] = orig->rhsx[m];
			prob[t]->sp->senx[k] = orig->senx[m];
			prob[t]->sp->rname[k] = orig->rname[m] + (prob[t]->sp->rstore - orig->rname[tim->row[t]]);
			prob[t]->bBar->val[prob[t]->bBar->cnt+1] = orig->rhsx[m];
			prob[t]->bBar->col[prob[t]->bBar->cnt+1] = m - tim->row[t]+1;
			prob[t]->bBar->cnt++;
		}
	}

	/* Now copy the terminal stage problem */
	/* copy column information for terminal stage*/
	m = 0;
	for ( q = orig->cname[tim->col[t]]; q < orig->cname[0] + orig->cstorsz; q++ )
		prob[t]->sp->cstore[m++] = *q;

	//cOffset = prob[t]->sp->cstore - orig->cname[tim->col[t]];
	for ( m = tim->col[t]; m < orig->mac; m++ ) {
		k = m - tim->col[t];
		prob[t]->dBar->val[prob[t]->dBar->cnt+1] = orig->objx[m];
		prob[t]->dBar->col[prob[t]->dBar->cnt+1] = m-tim->col[t]+1;
		prob[t]->dBar->cnt++;
		prob[t]->sp->objx[k] = orig->objx[m];
		prob[t]->sp->bdl[k] = orig->bdl[m];
		prob[t]->sp->bdu[k] = orig->bdu[m];
        //Jiajun, don't check for 'C'
//        if (orig->ctype[m] != 'C') {
//            errMsg("setup", "newProb", "integer variable in non-root stage",0);
//            return NULL;
//        }
//        else
//            prob[t]->sp->ctype[k] = orig->ctype[m];
		prob[t]->sp->cname[k] = orig->cname[m] + (prob[t]->sp->cstore - orig->cname[tim->col[t]]);
		prob[t]->sp->matcnt[k] = 0;
		if ( orig->matcnt[m] > 0 )
			for ( i = orig->matbeg[m]; i < orig->matbeg[m]+orig->matcnt[m]; i++ ) {
				/* The coefficient is part of the current stage constraint matrix */
				if ( k == 0 )
					prob[t]->sp->matbeg[k] = 0;
				else
					prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1] + prob[t]->sp->matcnt[k-1];
				prob[t]->sp->matval[prob[t]->sp->matsz] = orig->matval[i];
				prob[t]->sp->matind[prob[t]->sp->matsz] = orig->matind[i]-tim->row[t];
				++prob[t]->sp->matcnt[k];
				++prob[t]->sp->matsz;
				++prob[t]->sp->numnz;
			}
		else {
			if ( k == 0 )
				prob[t]->sp->matbeg[k] = 0;
			else
				prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1];
		}
	}

	/* copy row information for terminal stage */
	m = 0;
	for ( q = orig->rname[tim->row[t]]; q < orig->rname[0] + orig->rstorsz; q++ )
		prob[t]->sp->rstore[m++] = *q;

	//rOffset = prob[t]->sp->rstore - orig->rname[tim->row[t]];
	for ( m = tim->row[t]; m < orig->mar; m++ ) {
		k = m - tim->row[t];
		prob[t]->sp->rhsx[k] = orig->rhsx[m];
		prob[t]->sp->senx[k] = orig->senx[m];
		prob[t]->sp->rname[k] = orig->rname[m] + (prob[t]->sp->rstore - orig->rname[tim->row[t]]);
		prob[t]->bBar->val[prob[t]->bBar->cnt+1] = orig->rhsx[m];
		prob[t]->bBar->col[prob[t]->bBar->cnt+1] = m - tim->row[t]+1;
		prob[t]->bBar->cnt++;
	}

#ifdef DECOMPOSE_CHECK
//	/* write stage problems in LP format to verify decomposition */
//	char fname[BLOCKSIZE];
//	for ( t = 0; t < tim->numStages; t++) {
//		if ( !(prob[t]->sp->lp = setupProblem(prob[t]->sp->name, prob[t]->sp->type, prob[t]->sp->mac, prob[t]->sp->mar, prob[t]->sp->objsen, prob[t]->sp->objx, prob[t]->sp->rhsx, prob[t]->sp->senx,
//				prob[t]->sp->matbeg, prob[t]->sp->matcnt, prob[t]->sp->matind, prob[t]->sp->matval, prob[t]->sp->bdl, prob[t]->sp->bdu, NULL, prob[t]->sp->cname, prob[t]->sp->rname, prob[t]->sp->ctype)) ) {
//			errMsg("solver", "newProb", "failed to setup stage problem in solver", 0);
//			return prob;
//		}
//		sprintf(fname, "stageProb%d.lp", t);
//		if ( writeProblem(prob[t]->sp->lp, fname) ) {
//			errMsg("solver", "newProb", "failed to write stage problem", 0);
//			return prob;
//		}
//	}
#endif

	/* save size information in numType */
	for ( t = 0; t < tim->numStages; t++ ) {
		if ( !(prob[t]->num = (numType *) mem_malloc(sizeof(numType))) )
			errMsg("allocation", "newProb", "prob[t]->num",0);

		prob[t]->num->cols = prob[t]->sp->mac;
		prob[t]->num->rows = prob[t]->sp->mar;
		prob[t]->num->intCols = prob[t]->sp->numInt;
		prob[t]->num->binCols = prob[t]->sp->numBin;
		prob[t]->num->numRV = prob[t]->num->rvColCnt = prob[t]->num->rvRowCnt = 0;
		prob[t]->num->rvAOmCnt = prob[t]->num->rvBOmCnt = prob[t]->num->rvCOmCnt = prob[t]->num->rvDOmCnt = 0;
		prob[t]->num->rvaOmCnt = prob[t]->num->rvbOmCnt = prob[t]->num->rvcOmCnt = prob[t]->num->rvdOmCnt = 0;

		if ( t == 0 ) {
			prob[t]->mean = NULL;
			prob[t]->coord = NULL;
			prob[t]->num->prevCols = prob[t]->num->prevRows = prob[t]->num->cntCcols = prob[t]->num->cntCrows = 0;
		}
		else {
			if ( !(prob[t]->coord = (coordType *) mem_malloc(sizeof(coordType))) )
				errMsg("allocation", "newProb", "prob[t]->coord",0);
			prob[t]->num->prevCols = prob[t-1]->num->cols;
			prob[t]->num->prevRows = prob[t-1]->num->rows;
			prob[t]->coord->CCols = findElems(prob[t]->Cbar->col, prob[t]->Cbar->cnt, &prob[t]->num->cntCcols);
			prob[t]->coord->CRows = findElems(prob[t]->Cbar->row, prob[t]->Cbar->cnt, &prob[t]->num->cntCrows);
			prob[t]->coord->allRVCols = prob[t]->coord->allRVRows = prob[t]->coord->rvCols = prob[t]->coord->rvRows = NULL;
			prob[t]->coord->rvCOmCols = prob[t]->coord->rvCOmRows = prob[t]->coord->rvbOmRows = prob[t]->coord->rvdOmCols = NULL;
		}
	}

	/* decompose the stochastic elements of the problem. Go through the list of random variable and assign them to
	 * appropriate parts (right-hand side and objective coefficients). */
	for ( m = 0; m < stoc->numOmega; m++ ) {
		if ( stoc->col[m] == -1 ) {
			/* randomness in right-hand side */
			t = 0;
			while ( t < tim->numStages ) {
				if ( stoc->row[m] < tim->row[t] )
					break;
				t++;
			}
			t--;
		}
		else {
			/* randomness in either objective function coefficients or the transfer matrix */
			t = 0;
			while ( t < tim->numStages ) {
				if ( stoc->col[m] < tim->col[t] )
					break;
				t++;
			}
			t--;
		}

		if ( t == 0 ) {
			errMsg("setup", "newProb", "encountered randomness in root-stage", 0);
			return NULL;
		}

		if ( prob[t]->num->numRV == 0) {
			if ( !(prob[t]->coord->allRVCols = (intvec) arr_alloc(stoc->numOmega+1, int)) )
				errMsg("allocation", "newProb", "prob->coord->allRVCols", 0);
			if ( !(prob[t]->coord->allRVRows= (intvec) arr_alloc(stoc->numOmega+1, int)) )
				errMsg("allocation", "newProb", "prob->coord->allRVRows", 0);
			if ( !(prob[t]->coord->rvOffset = (intvec) arr_alloc(3, int)))
				errMsg("allocation", "newProb", "prob->coord->rOffset", 0);
			if ( !(prob[t]->mean = (vector) arr_alloc(stoc->numOmega+1, double)) )
				errMsg("allocation", "newProb", "prob->mean", 0);
			prob[t]->omBeg = m;
		}

		prob[t]->num->numRV++;
		prob[t]->mean[prob[t]->num->numRV] = stoc->mean[m];

		if ( stoc->col[m] == -1 && stoc->row[m] != -1 ) {
			/* Right-hand side */
			if ( prob[t]->num->rvbOmCnt == 0 ) {
				prob[t]->coord->rvbOmRows = (intvec) arr_alloc(stoc->numOmega+1, int);
				prob[t]->coord->rvOffset[0] = m;
			}
			prob[t]->coord->allRVCols[prob[t]->num->numRV] = -1;
			prob[t]->coord->allRVRows[prob[t]->num->numRV] = stoc->row[m]-tim->row[t]+1;
			prob[t]->coord->rvbOmRows[++prob[t]->num->rvbOmCnt] = prob[t]->coord->allRVRows[prob[t]->num->numRV];
		}
		else if ( stoc->col[m] != -1 && stoc->row[m] != -1 ) {
			/* Transfer matrix */
			if ( prob[t]->num->rvCOmCnt == 0 ) {
				prob[t]->coord->rvCOmCols = (intvec) arr_alloc(stoc->numOmega, int);
				prob[t]->coord->rvCOmRows = (intvec) arr_alloc(stoc->numOmega, int);
				prob[t]->coord->rvOffset[1] = m;
			}
			prob[t]->coord->allRVCols[prob[t]->num->numRV] = stoc->col[m]-tim->col[t]+1;
			prob[t]->coord->allRVRows[prob[t]->num->numRV] = stoc->row[m]-tim->row[t]+1;
			prob[t]->num->rvCOmCnt++;
			prob[t]->coord->rvCOmCols[prob[t]->num->rvCOmCnt] = prob[t]->coord->allRVCols[prob[t]->num->numRV];
			prob[t]->coord->rvCOmRows[prob[t]->num->rvCOmCnt] = prob[t]->coord->allRVRows[prob[t]->num->numRV];
		}
		else {
			/* Cost coefficients */
			if ( prob[t]->num->rvdOmCnt == 0 ) {
				prob[t]->coord->rvdOmCols = (intvec) arr_alloc(stoc->numOmega+1,int);
				prob[t]->coord->rvOffset[2] = m;
			}
			prob[t]->coord->allRVCols[prob[t]->num->numRV] = stoc->col[m]-tim->col[t]+1;
			prob[t]->coord->allRVRows[prob[t]->num->numRV] = -1;
			prob[t]->coord->rvdOmCols[++prob[t]->num->rvdOmCnt] = prob[t]->coord->allRVCols[prob[t]->num->numRV];
		}
	}

	for ( t = 1; t < tim->numStages; t++ ) {
		prob[t]->coord->rvCols = findElems(prob[t]->coord->allRVCols, prob[t]->num->numRV, &prob[t]->num->rvColCnt);
		prob[t]->coord->rvRows = findElems(prob[t]->coord->allRVRows, prob[t]->num->numRV, &prob[t]->num->rvRowCnt);
	}

	/* Modify the dBar, bBar and Cbar with mean values computed from stoch file */
	for ( t = 1; t < tim->numStages; t++ ) {
		/* Right-hand side */
		for ( m = 1; m <= prob[t]->num->rvbOmCnt; m++ ) {
			i = 1;
			while ( i <= prob[t]->bBar->cnt ) {
				if ( prob[t]->bBar->col[i] == prob[t]->coord->rvbOmRows[m])
					break;
				i++;
			}
			prob[t]->bBar->val[i] = stoc->mean[prob[t]->coord->rvOffset[0]+m-1];
		}

		/* Transfer matrix */
		for ( m = 1; m <= prob[t]->num->rvCOmCnt; m++ ) {
			i = 1;
			while ( i <= prob[t]->num->rvCOmCnt ) {
				if ( prob[t]->Cbar->col[i] == prob[t]->coord->rvCOmCols[m] && prob[t]->Cbar->row[i] == prob[t]->coord->rvCOmRows[m] )
					break;
				i++;
			}
			prob[t]->Cbar->val[i] = stoc->mean[prob[t]->coord->rvOffset[1]+m-1];
		}

		/* Cost coefficients */
		for ( m = 1; m <= prob[t]->num->rvdOmCnt; m++ ) {
			i = 1;
			while ( i <= prob[t]->dBar->cnt ) {
				if ( prob[t]->dBar->col[i] == prob[t]->coord->rvdOmCols[m])
					break;
				i++;
			}
			prob[t]->dBar->val[i] = stoc->mean[prob[t]->coord->rvOffset[2]+m-1];
		}
	}

#if defined(DECOMPOSE_CHECK)
	t = 1;
	prob[t]->sp->lp = setupProblem(prob[t]->sp->name, prob[t]->sp->type, prob[t]->sp->mac, prob[t]->sp->mar, prob[t]->sp->objsen, prob[t]->sp->objx, prob[t]->sp->rhsx, prob[t]->sp->senx,
			prob[t]->sp->matbeg, prob[t]->sp->matcnt, prob[t]->sp->matind, prob[t]->sp->matval, prob[t]->sp->bdl, prob[t]->sp->bdu, NULL, prob[t]->sp->cname, prob[t]->sp->rname,
			prob[t]->sp->ctype);
	if ( prob[t]->sp->lp == NULL ) {
		errMsg("solver", "newSubprob", "subprob",0);
		return NULL;
	}
#endif

	return prob;
}//END newProb()

/* setup and solve the original problem _orig_ with expected values for all random variables provided in _stoc_. If the problem is an mixed-integer program,
 *  then a relaxed problem is solved. The function returns a vector of mean value solutions, if there is an error it returns NULL.
 dim= the dimension of the decision variable in first stage*/
vector meanProblem(oneProblem *orig, stocType *stoc, int dim) {
	vector	xk;
	double	obj = 0.0;
	int 	n, status;

	/* setup problem in the solver */
	orig->lp = setupProblem(orig->name, orig->type, orig->mac, orig->mar, orig->objsen, orig->objx, orig->rhsx, orig->senx, orig->matbeg, orig->matcnt,
			orig->matind, orig->matval, orig->bdl, orig->bdu, NULL, orig->cname, orig->rname, orig->ctype);
	if ( orig->lp == NULL ) {
		errMsg("setup", "meanProblem", "failed to setup the mean problem", 0);
		return NULL;
	}

	/* change the coefficients and right-hand side to mean values */
    /*Jiajun: check stoc->mean*/
	for (n = 0; n < stoc->numOmega; n++ ) {
		status = changeCoef(orig->lp, stoc->row[n], stoc->col[n], stoc->mean[n]);
		if ( status ) {
			errMsg("setup", "meanProblem", "failed to change the coefficients with mean values", 0);
			return NULL;
		}
	}

	/* change the problem type to solve the relaxed mean value problem */
//    if ( orig->type == PROB_MILP || orig->type == PROB_MIQP ) {
//        status = changeProbType(orig->lp, PROB_LP);
//        if ( status ) {
//            errMsg("solver", "meanProblem", "failed to relax the mixed-integer program", 0);
//            return NULL;
//        }
//    }

	/* write the mean value problem */
	status = writeProblem(orig->lp, "original.lp");
	if ( status ) {
		errMsg("solver", "meanProblem", "failed to write the problem", 0);
		return NULL;
	}

	/* solve the mean value problem */
	changeLPSolverType(ALG_NET);
	status = solveProblem(orig->lp, orig->name, orig->type, &status);
	if ( status ) {
		errMsg("setup", "meanProblem", "failed to solve mean value problem", 0);
		return NULL;
	}

	/* obtain solution information and print */
	if ( !(xk = (vector) arr_alloc(dim+1, double)) )
		errMsg("allocation", "meanProblem", "sol", 0);

	/* print results */
	obj = getObjective(orig->lp, orig->type);
	printf("Optimal objective function value for mean value problem (Problem Type %d) = %lf\n", orig->type, obj);

	/* obtain the primal solution */
	getPrimal(orig->lp,	xk, dim);
    
	return xk;
}//END meanProblem()

vector calcLowerBound(oneProblem *orig, timeType *tim, stocType *stoc) {
	sparseVector	*bBar;
	sparseMatrix	*Cbar;
	vector		duals, vals, beta, lb;
	intvec		indices;
	double		alpha;
	int 		stat1, t, col, row, m, n;
	LPptr		lpClone;
	BOOL		zeroLB;

	if ( !(lb = (vector) arr_alloc(tim->numStages, double)) )
		errMsg("allocation", "getLowerBound", "lb", 0);
	if ( !(duals = (vector) arr_alloc(orig->mar+1, double)) )
		errMsg("allocation", "getLowerBound", "duals", 0);
	if (!(indices = (intvec) arr_alloc(orig->mac, int)))
		errMsg("allocation", "getLowerBound", "indices", 0);
	if (!(vals = (vector) arr_alloc(orig->mac, double)))
		errMsg("allocation", "getLowerBound", "vals", 0);
	if (!(bBar = (sparseVector *) mem_malloc(sizeof(sparseVector))) )
		errMsg("allocation", "getLowerBound", "bBar", 0);
	if (!(bBar->col = (intvec) arr_alloc(orig->mar+1, int)) )
		errMsg("allocation", "getLowerBound", "bBar->col", 0);
	if (!(bBar->val = (vector) arr_alloc(orig->mar+1, double)) )
		errMsg("allocation", "getLowerBound", "bBar->val", 0);
	if (!(Cbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
		errMsg("allocation", "getLowerBound", "Cbar", 0);
	if (!(Cbar->col = (intvec) arr_alloc(orig->matsz+1, int)) )
		errMsg("allocation", "getLowerBound", "Cbar->col", 0);
	if (!(Cbar->row = (intvec) arr_alloc(orig->matsz+1, int)) )
		errMsg("allocation", "getLowerBound", "Cbar->row", 0);
	if (!(Cbar->val = (vector) arr_alloc(orig->matsz+1, double)) )
		errMsg("allocation", "getLowerBound", "Cbar->val", 0);
	if ( !(beta = (vector) arr_alloc(orig->mac+3, double)) )
		errMsg("allocation", "getLowerBound", "beta", 0);

	/* obtain dual solutions from the mean value solve */
	if ( getDual(orig->lp, duals, orig->mar) ) {
		errMsg("setup", "getLowerBound", "failed to obtain dual for the mean value problem", 0);
		return NULL;
	}

	printf("Lower bounds computed = ");

	for ( t = 1; t < tim->numStages; t++ ) {
		zeroLB = TRUE;
		col = tim->col[t]; row = tim->row[t];
		bBar->cnt = 0; Cbar->cnt = 0;

		/* check to see if zero is a valid lower bound. If so, then we will use it as a lower bound */
		n = col;
		while ( n < orig->mac ) {
			if ( orig->objx[n]*orig->bdl[n] < 0 || orig->objx[n]*orig->bdu[n] < 0) {
				zeroLB = FALSE;
				break;
			}
			n++;
		}
		if ( !(zeroLB) ) {
			/* clone to problem to be used for computing the lower bound */
			lpClone = cloneProblem(orig->lp);

			/* extract bBar */
			m = 0;
			for (n = row; n < orig->mar; n++) {
				/* if the element has randomness in right-hand side, then make sure it is accounted */
				if ( n == stoc->row[m] && m < stoc->numOmega)
					bBar->val[bBar->cnt + 1] = stoc->mean[m++];
				else
					bBar->val[bBar->cnt + 1] = orig->rhsx[n];
				bBar->col[bBar->cnt + 1] = n - row + 1;
				++bBar->cnt;
			}

			/* extract Cbar */
			for (n = 0; n < col; n++) {
				for (m = orig->matbeg[n]; m < orig->matbeg[n] + orig->matcnt[n]; m++) {
					if (orig->matind[m] >= row) {
						/* The coefficient is part of the subproblem's T matrix */
						Cbar->val[Cbar->cnt + 1] = orig->matval[m];
						Cbar->row[Cbar->cnt + 1] = orig->matind[m] - row + 1;
						Cbar->col[Cbar->cnt + 1] = n + 1;
						++Cbar->cnt;
					}
				}
			}

			/* compute alpha and beta */
			alpha = 0.0;
			for ( n = 1; n <= bBar->cnt; n++ )
				alpha += duals[bBar->col[n]+row] * bBar->val[n];
			for (n = 0; n <= col + 2; n++)
				beta[n] = 0.0;
			for (n = 1; n <= Cbar->cnt; n++)
				beta[Cbar->col[n]] = beta[Cbar->col[n]] + duals[row + Cbar->row[n]] * Cbar->val[n];
			beta[col + 2] = -1;

			/* use mean-cut coefficients to create the lower-bounding problem */
			for (n = 0; n < col; n++) {
				indices[n] = n;
				vals[n] = -beta[n + 1];
			}
			for (n = col; n < orig->mac; n++) {
				indices[n] = n;
				vals[n] = 0.0;
			}

			/* Change the objective in solver to prepare for lower bound calculation */
			if ( changeObjx(lpClone, orig->mac, indices, vals) ) {
				errMsg("setup", "calcLowerBound", "failed to change objective coefficients in solver while computing lower bound", 0);
				return NULL;
			}

#ifdef SETUP_CHECK
			writeProblem(lpClone, "lowerBoundCalc.lp");
#endif

			/* solve the problem */
			changeLPSolverType(ALG_AUTOMATIC);
			if (solveProblem(lpClone, "lowerBoundCalc", PROB_LP, &stat1) ) {
				errMsg("setup", "calcLowerBound", "failed to solve problem computing lower bound", 0);
				return NULL;
			}

			/* get lower bound */
			lb[t-1] = min(0, getObjective(lpClone, PROB_LP) + alpha);

			/* release the problem */
			if ( freeProblem(lpClone) ) {
				errMsg("setup", "calcLowerBound", "failed to free problem", 0);
				return NULL;
			}
		}
		else
			lb[t-1] = 0.0;

		printf("%0.3lf\t", lb[t-1]);
	}
	lb[t-1] = 0.0;
	printf("\n");

	mem_free(indices);
	mem_free(vals);
	mem_free(beta);
	mem_free(duals);
	freeSparseVector(bBar); freeSparseMatrix(Cbar);

	/* change the problem type back to its original form */
	if ( orig->type == PROB_MILP ) {
		if ( changeProbType(orig->lp, PROB_MILP) ) {
			errMsg("solver", "calcLowerBound", "failed to relax the mean value problem", 0);
			return NULL;
		}
	}

	return lb;
}//END calcLowerBound()

/* free up the probType. Needs number of stages as input */
void freeProbType(probType **prob, int T) {
	int t;

	if ( prob ) {
		for ( t = 0; t < T; t++ ) {
			if (prob[t]) {
				if (prob[t]->Abar) freeSparseMatrix(prob[t]->Abar);
				if (prob[t]->Bbar) freeSparseMatrix(prob[t]->Bbar);
				if (prob[t]->Cbar) freeSparseMatrix(prob[t]->Cbar);
				if (prob[t]->Dbar) freeSparseMatrix(prob[t]->Dbar);
				if (prob[t]->aBar) freeSparseVector(prob[t]->aBar);
				if (prob[t]->bBar) freeSparseVector(prob[t]->bBar);
				if (prob[t]->cBar) freeSparseVector(prob[t]->cBar);
				if (prob[t]->dBar) freeSparseVector(prob[t]->dBar);
				if (prob[t]->sp) freeOneProblem(prob[t]->sp);
				if (prob[t]->name) mem_free(prob[t]->name);
				if (prob[t]->num) mem_free(prob[t]->num);
				if (prob[t]->coord) freeCoordType(prob[t]->coord);
				if (prob[t]->mean) mem_free(prob[t]->mean);
				mem_free(prob[t]);
			}
		}
		mem_free(prob);
	}

}//END freeProb()

/* free up the coordType */
void freeCoordType (coordType *coord) {

	if (coord->allRVCols) mem_free(coord->allRVCols);
	if (coord->allRVRows) mem_free(coord->allRVRows);
	if (coord->CCols) mem_free(coord->CCols);
	if (coord->CRows) mem_free(coord->CRows);
	if (coord->rvCols) mem_free(coord->rvCols);
	if (coord->rvRows) mem_free(coord->rvRows);
	if (coord->rvbOmRows) mem_free(coord->rvbOmRows);
	if (coord->rvdOmCols) mem_free(coord->rvdOmCols);
	if (coord->rvCOmCols) mem_free(coord->rvCOmCols);
	if (coord->rvCOmRows) mem_free(coord->rvCOmRows);
	if (coord->rvOffset) mem_free(coord->rvOffset);
	mem_free(coord);

}//END freeCoordType()

void printDecomposeSummary(FILE *fptr, string probName, timeType *tim, probType **prob) {
	int t;

	fprintf(fptr, "====================================================================================================================================\n");
	fprintf(fptr, "Stage optimization problem for given input s_t = (x_t, omega_t) is in the following form : \n\n");
	fprintf(fptr, "\t\t\t\t h_t(s_t) = c_t*x_t + min d_t*u_t\n");
	fprintf(fptr, "\t\t\t\t                      s.t. D_t u_t = b_t - C_t x_t,\n\n");
	fprintf(fptr, "with the following linear dynamics: x_{t+} = a_{t+} + A_{t+}x_t + B_{t+}u_t.\n\n");
	fprintf(fptr, "------------------------------------------------------------------------------------------------------------------------------------\n");
	fprintf(fptr, "\n====================================================================================================================================\n");
	fprintf(fptr, "-------------------------------------------------------- Problem Information -------------------------------------------------------\n");
	fprintf(fptr, "====================================================================================================================================\n");
	fprintf(fptr, "Problem                            : %s\n", probName);
	fprintf(fptr, "Number of stages                   : %d\n", tim->numStages);
	for ( t = 0; t < tim->numStages; t++ ) {
		fprintf(fptr,  "------------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(fptr,  "Stage %d\n", t+1);
		fprintf(fptr,  "Number of decision variables (u_t) = %d\t\t", prob[t]->sp->mac);
		fprintf(fptr,  "(Continuous = %d\tInteger = %d\tBinary = %d)\n", prob[t]->sp->mac - prob[t]->sp->numInt - prob[t]->sp->numBin, prob[t]->sp->numInt, prob[t]->sp->numBin);
		fprintf(fptr,  "Number of constraints              = %d\n", prob[t]->sp->mar);
		if ( prob[t]->num->numRV != 0 ) {
			fprintf(fptr,  "Number of random variables (omega) = %d\t\t", prob[t]->num->numRV);
			fprintf(fptr,  "(a_t = %d; b_t = %d; c_t = %d; d_t = %d; A_t = %d; B_t = %d; C_t = %d; D_t = %d)\n", prob[t]->num->rvaOmCnt, prob[t]->num->rvbOmCnt, prob[t]->num->rvcOmCnt, prob[t]->num->rvdOmCnt,
					prob[t]->num->rvAOmCnt, prob[t]->num->rvBOmCnt, prob[t]->num->rvCOmCnt, prob[t]->num->rvDOmCnt);
		}
		else
			fprintf(fptr,  "Number of random variables (omega) = 0\n");
	}
	fprintf(fptr, "====================================================================================================================================\n");

}//printDecomposeSummary()
