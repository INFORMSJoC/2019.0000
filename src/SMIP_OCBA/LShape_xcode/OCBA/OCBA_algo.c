/*
 * OCAB_algo.c
 *
 *  Created on: March 16, 2020
 *      Author: Jiajun Xu
 * Institution: University of Southern California
 *
 * Please send your comments or bug report to jiajunx (at) usc (dot) edu
 *
 */
#include <math.h>
#include <stdio.h>
#include "OCBA_algo.h"
#define DEBUG
extern string outputDir;
extern configType config;

#if defined(SAVE_DUALS)
dualsType *duals = NULL;
#endif

int algo (oneProblem *orig, timeType *tim, stocType *stoc, string probName) {
    probType **prob = NULL;
    cellType *cell = NULL;
    vector      meanSol;
    double   lb_mean=0, lb_std=0, ocba_time=0, naive_time=0, stdev=0, pr, temp;
    batchSummary *batch = NULL;
    int      rep, m, n, out_idx=0, idx;
    FILE     *sFile, *iFile = NULL, *fFile, *ocbaFile;
    char results_name[BLOCKSIZE];
    char incumb_name[BLOCKSIZE];
    clock_t    tic, tic_time;
    
    BOOL distinct;
    double   inverse_appearance[config.NUM_REPS];
    int Delta = 500;
    int     i, early_stop=0, bat=0;
    int ocba_samples=0, naive_samples=0;
    int candidates_per_batch=2;
    ocbaSummary *ocba = NULL;
    
    
    /* complete necessary initialization for the algorithm */
    if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol) )
        goto TERMINATE;

    if ( config.NUM_REPS > 1 )
        ocba  = newOcbaSummary(prob[0]->num->cols, config.NUM_REPS*candidates_per_batch);
        
    printf("Starting Benders decomposition.\n");
    sprintf(results_name, "rep_results%d.txt", out_idx);
    sFile = openFile(outputDir, results_name, "w");
    sprintf(results_name, "final_results%d.txt", out_idx);
    fFile = openFile(outputDir, results_name, "w");
    sprintf(incumb_name, "incumb%d.txt", out_idx);
    iFile = openFile(outputDir, incumb_name, "w");
    ocbaFile = openFile(outputDir, "ocba_results.txt", "w");
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
        
        tic = clock();
        /* Use two-stage algorithm to solve the problem */
        if ( solveBendersCell(stoc, prob, cell, ocba, rep, candidates_per_batch) ) {
            errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
            goto TERMINATE;
        }
        cell->time->repTime = ((double) (clock() - tic))/CLOCKS_PER_SEC;

        /* Write solution statistics for optimization process */
        writeStatistic(sFile, iFile, prob, cell);
        fclose(sFile);
        writeStatistic(stdout, NULL, prob, cell);

        /* check is the new solution is distict or not. */
        if ( config.MULTIPLE_REP )
            buildCompromise(prob[0], cell, batch);
    }//end all replciation run
    
    if ( config.MULTIPLE_REP ) {
        print_ocba_results(ocbaFile, ocba, prob, candidates_per_batch);
            
        for (idx = 0; idx < candidates_per_batch*config.NUM_REPS; idx++){
            if(ocba->objLB_list[idx] == -INF)
                continue;
            
            distinct = TRUE;
            for(i = 0; i < ocba->cnt; i++ ){
                if (equalVector(ocba->incumbX[i], ocba->incumbX_list[idx], prob[0]->num->cols, config.TOLERANCE)){
                    distinct = FALSE;
                    break;
                }
            }
            if(distinct == TRUE){
                ocba->objLB[ocba->cnt] = ocba->objLB_list[idx];
                ocba->incumbX[ocba->cnt] = duplicVector(ocba->incumbX_list[idx], prob[0]->num->cols);
                ocba->appearance[ocba->cnt] = 1;
                ocba->cnt++;
            }
            else{
                ocba->appearance[i]++;
                temp = ocba->objLB[i];
                ocba->objLB[i] = (ocba->appearance[i] - 1.0) / ocba->appearance[i] * temp + ocba->objLB_list[idx] / ocba->appearance[i];
            }
        }
        
        
        for (i = 0; i < ocba->cnt; i++){
            inverse_appearance[i] = 1.0 / ocba->appearance[i];
            printf("inverse_appearance[%d]=%f\n", i, inverse_appearance[i]);
        }
        
        
        tic = clock();
        
        //if only one solution
        if (ocba->cnt == 1){
            naive_samples += evaluate_samples(fFile, stoc, prob, cell, ocba->incumbX[0]);
            ocba_samples = naive_samples;
            naive_time = ((double) (clock() - tic)) / CLOCKS_PER_SEC;
            ocba_time = naive_time;
            ocba->idx = 0;
            
        }
        else{
            /*1. Naive SAA evaluate all solutions */
            for (i = 0; i < ocba->cnt; i++){
                naive_samples += evaluate_samples(fFile, stoc, prob, cell, ocba->incumbX[i]);
                printVectorInSparse(ocba->incumbX[i], prob[0]->num->cols, fFile);
            }
            naive_time = ((double) (clock() - tic)) / CLOCKS_PER_SEC;
            
            /* 2. Solve the ocba problem. */
            printf("\n ocba->cnt=%d", ocba->cnt);
            ocba->idx = solveOCBA(ocba->objLB, inverse_appearance, ocba->cnt, ocba->n, Delta, ocba->an);
            for (i = 0; i < ocba->cnt; i++)
                ocba->an[i] += 10;
            
            tic = clock();
            early_stop = eval_all(fFile, stoc, prob, cell, ocba, &pr);
            tic_time = clock() - tic;
            
            if (early_stop){
                ocba_samples += early_stop;
            }
            else{
                ocba_samples += Delta;
            }
            
            while(early_stop == 0 ){
               
                ocba->idx = solveOCBA(ocba->mean, ocba->var, ocba->cnt, ocba->n, Delta, ocba->an);
                
                
                tic = clock();
                early_stop = eval_all(fFile, stoc, prob, cell, ocba, &pr);

                tic_time += clock() - tic;
                
                if (early_stop){
                    ocba_samples += early_stop;
                    break;
                }
                else{
                    ocba_samples += Delta;
                }
                
            }
            
            ocba_time = ((double) tic_time) / CLOCKS_PER_SEC;
            
        }
        batch->time->repTime += ocba_time;

        fprintf(fFile, "\n----------------------------------------- Final solution --------------------------------------\n\n");
        /* Evaluate the compromise solution */
        fprintf(fFile, "OCBA solution, non-zero position: ");
        printVectorInSparse(ocba->incumbX[ocba->idx], prob[0]->num->cols, fFile);

        
        
        fprintf(fFile, "Time to solve ocba problem   : %f\n", ocba_time);
        fprintf(fFile, "Time for saive SAA evaluation   : %f\n", naive_time);
        fprintf(fFile, "Total time                         : %f\n", batch->time->repTime);
        fprintf(fFile, "Total time to solve master         : %f\n", batch->time->masterAccumTime);
        fprintf(fFile, "Total time to solve subproblems    : %f\n", batch->time->subprobAccumTime);
        
        
        LowerBoundVariance(batch, &lb_mean, &lb_std);
        fprintf(fFile, "Lower bound estimate               : %f\n", lb_mean);
        fprintf(fFile, "Lower bound estimation std               : %f\n", lb_std);
        
        
        fprintf(iFile, "\n----------------------------------------- Final solution by OCBA method(1-indexed) --------------------------------------\n\n");
        printVectorInSparse(ocba->incumbX[ocba->idx], prob[0]->num->cols, iFile);
        fclose(iFile);
        
        
        stdev = sqrt(ocba->var[ocba->idx] / ocba->n[ocba->idx]);
        
        fprintf(fFile, "\n ocba_samples=&%d, naive_samples=&%d", ocba_samples, naive_samples);
        fprintf(fFile, "\n print for latex\n");
        fprintf(fFile, "%d      %.2f%%  [%.3f, %.3f]    [%.3f, %.3f]     %.4f", ocba_samples, 100.0*(naive_samples - ocba_samples)/naive_samples, lb_mean-1.96 * lb_std,lb_mean+1.96 * lb_std, ocba->mean[ocba->idx]-1.96 * stdev, ocba->mean[ocba->idx]+1.96 * stdev, pr);
        printf("\n%d      %.2f%%  [%.3f, %.3f]    [%.3f, %.3f]     %.4f\n", ocba_samples, 100.0*(naive_samples - ocba_samples)/naive_samples, lb_mean-1.96 * lb_std,lb_mean+1.96 * lb_std, ocba->mean[ocba->idx]-1.96 * stdev, ocba->mean[ocba->idx]+1.96 * stdev, pr);
        
        fprintf(ocbaFile, "\n\n solution    sample_size mean    stdev\n");
        for (i = 0; i < ocba->cnt; i++){
            printVectorInSparse(ocba->incumbX[i], prob[0]->num->cols, ocbaFile);
            fprintf(ocbaFile, " %d  %f  %f\n\n", ocba->n[i], ocba->mean[i], sqrt(ocba->var[i] / ocba->n[i]));
        }
    }

    fclose(fFile);fclose(ocbaFile);
    


    /* free up memory before leaving */
    freeOCBA(ocba, config.NUM_REPS*candidates_per_batch);
    freeCellType(cell);
    freeProbType(prob, 2);
//    mem_free(meanSol);
    return 0;

    TERMINATE:
    if(cell) freeCellType(cell);
    if(prob) freeProbType(prob, 2);
    mem_free(meanSol);
    freeOCBA(ocba, config.NUM_REPS*candidates_per_batch);
    return 1;
}//END algo()

int solveBendersCell(stocType *stoc, probType **prob, cellType *cell, ocbaSummary *ocba, int rep, int candidates_per_batch) {
    int     candidCut;
    clock_t    tic, mainTic;
    int idx=0;
//    int lowest_obj_index;
//    double lowest_obj_val;
//    BOOL added;

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
        if (optimal(cell))
            break;
        
        ocba->objLB_list[rep*candidates_per_batch+idx] = cell->incumbEst;
        copyVector(cell->candidX, ocba->incumbX_list[rep*candidates_per_batch+idx], prob[0]->num->cols, TRUE);
        idx++;
        idx = idx % candidates_per_batch;
        

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

void print_ocba_results(FILE *ocbaFile, ocbaSummary *ocba, probType **prob, int candidates_per_batch){
    int rep, idx;
    fprintf(ocbaFile, "rep  sol1    obj1    sol2    obj2    sol3    obj3    \n");
    for ( rep = 0; rep < config.NUM_REPS; rep++){
        fprintf(ocbaFile, "%d   ", rep);
        for( idx = 0; idx < candidates_per_batch; idx++){
            if(ocba->objLB_list[rep*candidates_per_batch + idx] == -INF)
                fprintf(ocbaFile, "NaN  NaN ");
            else{
                printVectorInSparse(ocba->incumbX_list[rep*candidates_per_batch + idx], prob[0]->num->cols, ocbaFile);
                fprintf(ocbaFile, " %f   ", ocba->objLB_list[rep*candidates_per_batch + idx]);
            }
        }
        fprintf(ocbaFile, "\n");
    }
}

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

}//END WriteStat

ocbaSummary *newOcbaSummary(int first_stage_cols, int numCand) {
    ocbaSummary *ocba;
    int i;
    
    ocba = (ocbaSummary *) mem_malloc(sizeof(ocbaSummary));
    
    ocba->cnt = 0;
    ocba->idx = 0;
    
    ocba->ck = (intvec) arr_alloc(numCand, int);
    ocba->objLB = (vector) arr_alloc(numCand, double);
    ocba->objLB_list = (vector) arr_alloc(numCand, double);

    ocba->incumbX = (vector *) arr_alloc(numCand, vector);
    ocba->incumbX_list = (vector *) arr_alloc(numCand, vector);
    ocba->appearance = (intvec) arr_alloc(numCand, int);
    
    ocba->mean = (vector) arr_alloc(numCand, double);
    ocba->var = (vector) arr_alloc(numCand, double);
    
    for (i = 0; i < numCand; i++){
        ocba->incumbX[i] = (vector) arr_alloc(first_stage_cols + 1, double);
        ocba->incumbX_list[i] = (vector) arr_alloc(first_stage_cols + 1, double);
        ocba->objLB_list[i] = -INF;
    }
    
    ocba->n = (intvec) arr_alloc(numCand, int);
    ocba->an = (intvec) arr_alloc(numCand, int);
    return ocba;
}

void freeOCBA(ocbaSummary *ocba, int numCand){
    int i;
    
    if ( ocba ) {
        if (ocba->ck) mem_free(ocba->ck);
        if (ocba->objLB) mem_free(ocba->objLB);
        if (ocba->objLB_list) mem_free(ocba->objLB_list);
        if (ocba->mean) mem_free(ocba->mean);
        if (ocba->var) mem_free(ocba->var);
        if (ocba->appearance) mem_free(ocba->appearance);
        for (i = 0; i < numCand; i++){
            if (ocba->incumbX[i]) mem_free(ocba->incumbX[i]);
            if (ocba->incumbX_list[i]) mem_free(ocba->incumbX_list[i]);
        }
        if (ocba->incumbX) mem_free(ocba->incumbX);
        if (ocba->incumbX_list) mem_free(ocba->incumbX_list);
        if(ocba->n) mem_free(ocba->n);
        if(ocba->an) mem_free(ocba->an);
        mem_free(ocba);
    }
}


int solveOCBA(vector s_mean, vector s_var, int nd, intvec n, int add_budget, intvec an){
    /* This subroutine determines how many additional runs each design will should have for next iteration of simulation.
    s_mean[i]: sample mean of design i, i=0,1,..,nd-1
    s_var[i]: sample variance of design i, i=0,1,..,nd-1
    nd: the number of designs
    n[i]: number of simulation replications of design i, i=0,1,..,nd-1
    add_budget: the additional simulation budget
    an[i]: additional number of simulation replications assigned to design i, i=0,1,..,nd-1 */
    
    int i;
    int b, s;
    int t_budget, t1_budget;
    int morerun[nd], more_alloc; /* 1:Yes; 0:No */
    double t_s_mean[nd];
    double ratio[nd]; /* Ni/Ns */
    double ratio_s, temp;
    
    for(i = 0; i < nd; i++)
        t_s_mean[i] = s_mean[i];
    
    /*t_budget record total budget*/
    t_budget = add_budget;
    for(i = 0; i < nd; i++)
        t_budget += n[i];
    b = best(t_s_mean, nd);
    s = second_best(t_s_mean, nd, b);

    /* calculate ratio of Ni/Ns*/
    ratio[s] = 1.0;
    for(i = 0; i < nd; i++)
        if(i != s && i != b){
            temp = (t_s_mean[b] - t_s_mean[s]) / (t_s_mean[b] - t_s_mean[i]);
            ratio[i] = temp * temp * s_var[i] / s_var[s];
        }
    
    /* calculate Nb */
    temp = 0;
    for(i = 0; i < nd; i++)
        if(i != b)
            temp += (ratio[i] * ratio[i] / s_var[i]);
    ratio[b] = sqrt(s_var[b] * temp);
    
    
    for(i = 0; i < nd; i++)
        morerun[i] = 1;
    t1_budget = t_budget;
    
    do{
        more_alloc = 0;
        ratio_s = 0.0;
        for(i = 0; i < nd; i++)
            if(morerun[i])
                ratio_s += ratio[i];
        for(i = 0; i < nd; i++)
            if(morerun[i]) {
                an[i] = (int)(t1_budget / ratio_s * ratio[i]);
                
                /* disable those design which have been run too much */
                if(an[i] <= n[i]){
                    an[i] = n[i];
                    morerun[i] = 0;
                    more_alloc = 1;
                }
            }
        if(more_alloc) {
            t1_budget = t_budget;
            for(i = 0; i < nd; i++)
                if(!morerun[i])
                    t1_budget -= an[i];
        }
    } while(more_alloc); /* end of WHILE */
    
    /* calculate the difference */
    t1_budget = 0;
    for(i = 0; i < nd; i++)
        t1_budget += an[i];
    an[b] += (t_budget - t1_budget); /* give the difference to design b */
    for(i = 0; i < nd; i++)
        an[i] -= n[i];
    
#ifdef DEBUG
    for(i=0; i < nd; i++){
        printf("\n design=%d, s_mean=%.2f, s_var=%.2f, n=%d", i, s_mean[i], s_var[i], n[i]);
    }
#endif
    return b;
}



int eval_all(FILE *soln, stocType *stoc, probType **prob, cellType *cell, ocbaSummary *ocba, double *pr) {
    vector     observ, rhs, costTemp, cost;
    intvec    objxIdx;
    double     obj, mean, variance, stdev, temp;
    int        old_cnt, cnt, m, status, i, idx;

    if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
        errMsg("allocation", "evaluateOpt", "observ", 0);
    
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
    
    //evaluate cretical solution
    printf("\n Evaluate OCBA solution %d with obs %d.\n", ocba->idx, ocba->an[ocba->idx]);

    /* initialize parameters used for evaluations */
    cnt = ocba->n[ocba->idx]; mean = ocba->mean[ocba->idx] - vXvSparse(ocba->incumbX[ocba->idx], prob[0]->dBar); variance = ocba->var[ocba->idx]; stdev = INFBOUND; old_cnt = ocba->n[ocba->idx];

    /* change the right hand side with the solution */
    chgRHSwSoln(prob[1]->bBar, prob[1]->Cbar, rhs, ocba->incumbX[ocba->idx]);

    while (cnt < ocba->n[ocba->idx] + ocba->an[ocba->idx] ) {
        /* use the stoc file to generate observations */
        generateOmega(stoc, observ, config.TOLERANCE, &config.EVAL_SEED[0]);

        for ( m = 0; m < stoc->numOmega; m++ )
            observ[m] -= stoc->mean[m];

        /* Change right-hand side with random observation */
        if ( chgRHSwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, observ-1, rhs, ocba->incumbX[ocba->idx]) ) {
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
        
        //early stop
        if (cnt >= config.EVAL_MIN_ITER &&3.29 * stdev < config.EVAL_ERROR * DBL_ABS(mean)){
            *pr = 1.0;
            for (idx = 0; idx < ocba->cnt; idx++){
                if (idx != ocba->idx)
                    *pr -= 0.5 * erfc((ocba->mean[idx] - ocba->mean[ocba->idx]) / sqrt(2 * (ocba->var[idx]/ocba->n[idx] + ocba->var[ocba->idx]/ocba->n[ocba->idx])));
            }
            if (*pr > 0.7){
                mean += vXvSparse(ocba->incumbX[ocba->idx], prob[0]->dBar);
                ocba->n[ocba->idx] = cnt;
                ocba->mean[ocba->idx] = mean;
                ocba->var[ocba->idx] = variance;
                ocba->an[ocba->idx] = 0;
                
                stdev = sqrt(ocba->var[ocba->idx] / ocba->n[ocba->idx]);
                printf("\n mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", ocba->mean[ocba->idx], 3.92 * stdev / ocba->mean[ocba->idx],  ocba->mean[ocba->idx] - 1.96 * stdev, ocba->mean[ocba->idx] + 1.96 * stdev);
                mem_free(observ); mem_free(rhs); mem_free(objxIdx); mem_free(cost);
                return cnt - old_cnt;
            }
            else{
                config.EVAL_MIN_ITER += 2000;
            }
            
        }

        /* Print the results every once in a while for long runs */
        if (!(cnt % 100)) {
            printf(".");
            fflush(stdout);
        }
        if (!(cnt % 10000))
            printf("\nObs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", cnt, mean, 3.92 * stdev / mean,  mean - 1.96 * stdev, mean + 1.96 * stdev);
    }//END while loop
    mean += vXvSparse(ocba->incumbX[ocba->idx], prob[0]->dBar);
    ocba->n[ocba->idx] = cnt;
    ocba->mean[ocba->idx] = mean;
    ocba->var[ocba->idx] = variance;
    
    stdev = sqrt(ocba->var[ocba->idx] / ocba->n[ocba->idx]);
    printf("\n mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", ocba->mean[ocba->idx], 3.92 * stdev / ocba->mean[ocba->idx],  ocba->mean[ocba->idx] - 1.96 * stdev, ocba->mean[ocba->idx] + 1.96 * stdev);
    ocba->an[ocba->idx] = 0;
       
  //evaluate non-cretical solution
    for (i = 0; i < ocba->cnt; i++){
        
        if (i == ocba->idx)
            continue;
            
        printf("\nStarting evaluate solution %d with obs %d.\n", i, ocba->an[i]);

        /* initialize parameters used for evaluations */
        cnt = ocba->n[i]; mean = ocba->mean[i] - vXvSparse(ocba->incumbX[i], prob[0]->dBar); variance = ocba->var[i]; stdev = INFBOUND; old_cnt = ocba->n[i];


        /* change the right hand side with the solution */
        chgRHSwSoln(prob[1]->bBar, prob[1]->Cbar, rhs, ocba->incumbX[i]);

        while (cnt < ocba->n[i] + ocba->an[i] ) {
            /* use the stoc file to generate observations */
            generateOmega(stoc, observ, config.TOLERANCE, &config.EVAL_SEED[0]);

            for ( m = 0; m < stoc->numOmega; m++ )
                observ[m] -= stoc->mean[m];

            /* Change right-hand side with random observation */
            if ( chgRHSwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, observ-1, rhs, ocba->incumbX[i]) ) {
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
                printf("\nObs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", cnt, mean, 3.92 * stdev / mean,  mean - 1.96 * stdev, mean + 1.96 * stdev);
            
    
        }//END while loop
        mean += vXvSparse(ocba->incumbX[i], prob[0]->dBar);
        
        ocba->n[i] = cnt;
        ocba->mean[i] = mean;
        ocba->var[i] = variance;
        
        stdev = sqrt(ocba->var[i] / ocba->n[i]);
        printf("\n mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", ocba->mean[i], 3.92 * stdev / ocba->mean[i],  ocba->mean[i] - 1.96 * stdev, ocba->mean[i] + 1.96 * stdev);

    }
    mem_free(observ); mem_free(rhs); mem_free(objxIdx); mem_free(cost);
    return 0;

}//END de_all()

int best(vector t_s_mean, int nd){
    /*This function determines the best design based on current simulation results
     t_s_mean[i]: temporary array for sample mean of design i, i=0,1,..,ND-1
     nd: the number of designs */
    int i, min_index=0;
    for (i = 0; i < nd; i++){
        if(t_s_mean[i] < t_s_mean[min_index])
            min_index = i;
    }
    return min_index;
}

int second_best(vector t_s_mean, int nd, int b){
    int i, second_index;
    if (b==0)
        second_index = 1;
    else
        second_index = 0;
    for(i = 0;i < nd; i++){
        if(t_s_mean[i] < t_s_mean[second_index] && i != b)
        {
            second_index = i;
        }
    }
    return second_index;
}

