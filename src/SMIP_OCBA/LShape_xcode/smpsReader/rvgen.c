/*
 * rvgen.c
 *
 *  Created on: Feb 17, 2019
 *      Author: Jiajun Xu
 * Institution: University of Southern California
 *
 * Please send your comments or bug report to jiajunx (at) usc (dot) edu
 */

#include "utils.h"
#include "smps.h"

/* Currently works only for BLOCKS-DISCRETE. This subroutine generates the index of block according to the distribution provided in the stoch file. */
int generateOmegaIdx(stocType *stoc, long long *seed) {
	double 	val, cumm;
	int		m, groupID = 0;

	if ( strstr(stoc->type, "BLOCKS") != NULL ) {
		if ( strstr(stoc->type, "DISCRETE") != NULL  ) {
			val = scalit(0,1, seed);
			cumm = 0;
			for ( m = 0; val > cumm; m++ )
				cumm += stoc->probs[groupID][m];
			return (m-1);
		}
		else
			printf("Not yet\n");
	}
	else
		printf("Not yet\n");

	return 0;
}//END generateOmegaIdx()

void generateOmega(stocType *stoc, vector observ, double minVal, long long *seed) {
	int n, offset = 0;

	for ( n = 0; n < stoc->numGroups; n++ ) {
		if ( strstr(stoc->type, "INDEP") != NULL ) {
			if ( strstr(stoc->type, "DISCRETE") != NULL )
				generateIndep(stoc, observ+offset, n, seed);
			else if ( strstr(stoc->type, "NORMAL") != NULL )
				nnormal(stoc->mean+offset, stoc->vals[0]+offset, stoc->numPerGroup[n], observ+offset, seed);
			else
				errMsg("rvGeneration", "generateOmega", "random number generation for input type is missing",0);
			offset += stoc->numPerGroup[n];
		}
		else if ( strstr(stoc->type, "BLOCKS") != NULL ) {
			if ( strstr(stoc->type, "DISCRETE") != NULL  )
				generateBlocks(stoc, observ+offset, n, seed);
			else
				printf("Not yet\n");
			offset += stoc->numPerGroup[n];
		}
		else if ( strstr(stoc->type, "DISTRIB") != NULL ) {
			//	generateDistrib(observ+offset, n, seed);
			offset += stoc->numPerGroup[n];
		}
		else if ( strstr(stoc->type, "LINTRAN") != NULL) {
			generateLinTran(stoc, observ+offset, n, minVal, seed);
			n++; 			/* Linear transformation is associated with a block of residual random variables. */
		}
		else
			errMsg("rvgen", "generateOmega", "unknown section type in omegastuff", 0);
	}

}//END generateOmega()

void generateBlocks(stocType *stoc, vector observ, int groupID, long long *seed) {
	int 	n, m;
	double 	val, cumm;

	/* select which block */
	val = scalit(0,1, seed);
	cumm = 0;
	for ( m = 0; val > cumm; m++ )
		cumm += stoc->probs[groupID][m];

	/* read block realizations */
	for (n = 0; n < stoc->numPerGroup[groupID]; n++ )
		observ[n] = stoc->vals[stoc->groupBeg[groupID]+n][m-1];

}//END generateBlocks()

void generateIndep(stocType *stoc, vector observ, int groupID, long long *seed) {
	double 	val, cumm;
	int		n, m;

	for ( n = 0; n < stoc->numPerGroup[groupID]; n++) {
		val = scalit(0, 1, seed);
		cumm = 0;
		for ( m = 0; val > cumm && m < stoc->numVals[stoc->groupBeg[groupID]+n]; ++m )
			cumm += stoc->probs[stoc->groupBeg[groupID] + n][m];
		observ[stoc->groupBeg[groupID]+n] = stoc->vals[stoc->groupBeg[groupID]+n][m-1];
	}

}//END generateIndep()

/* Supporting only Wiener process */
void generateLinTran(stocType *stoc, vector observ, int groupID, double minVal, long long *seed) {
	vector eps;
	int n, offset, t;
	int numPeriods, epsLen, omegaLen;

	numPeriods = (int) stoc->numPerGroup[groupID]/stoc->mod->N;	/* Number of periods being simulated */
	epsLen     = stoc->mod->M;									/* Number of residual random variables */
	omegaLen   = stoc->mod->N;									/* Number of random variables per period */

	/* residual vector */
	eps = (vector) arr_alloc(numPeriods*epsLen, double);

	/* constant terms */
	offset = stoc->groupBeg[groupID];
	for (n = 0; n < stoc->numPerGroup[groupID]; n++ )
		observ[n] = stoc->vals[0][offset+n];

	for ( t = 0; t < numPeriods; t++ ) {
		nnormal(stoc->mod->muEps, stoc->mod->cvEps->val, epsLen, eps+epsLen*t, seed);
		if ( t == 0 ) {
			for (n = 0; n < stoc->mod->MA[0]->cnt; n++)
				observ[stoc->mod->MA[0]->row[n]] += stoc->mod->MA[0]->val[n] * eps[stoc->mod->MA[0]->col[n]];
		}
		else {
			for (n = 0; n < stoc->mod->MA[0]->cnt; n++)
				observ[omegaLen*t+stoc->mod->MA[0]->row[n]] += stoc->mod->MA[0]->val[n] * eps[epsLen*t+stoc->mod->MA[0]->col[n]];
			for (n = 0; n < stoc->mod->AR[0]->cnt; n++)
				observ[omegaLen*t+stoc->mod->AR[0]->row[n]] += stoc->mod->AR[0]->val[n] * observ[omegaLen*(t-1)+stoc->mod->AR[0]->col[n]];
		}
	}

	for (n = 0; n < stoc->numPerGroup[groupID]; n++ )
		observ[n] = max(observ[n], 0.01);

	mem_free(eps);
	return;
}//END generateLinTran()

/* The following inverse normal variate generator was published by Micheal J. Wichura, University of Chicago in Applied Statistics, as Algorithm AS 241.  The C function normal() was converted from the
 * Fortran function PPND7 and produces normal random variates for the lower tail of a normal distribution accurate to approx. 7 significant figures. */
int nnormal(vector mu, vector stdev, int numOmega, vector observ, long long *seed) {
	int i;
	float zero, one, half, split1, split2, const1, const2, a0, a1, a2, a3, b1;
	float b2, b3, c0, c1, c2, c3, d1, d2, e0, e1, e2, e3, f1, f2, p, q, r;
	float endval;

	for (i = 0; i < numOmega; i++) {
		p = scalit(0, 1, seed);

		zero = 0.0;
		one = 1.0;
		half = one / 2.0;
		split1 = 0.425;
		split2 = 5.0;
		const1 = 0.180625;
		const2 = 1.6;

		/* coefficients for p close to 1/2 */
		a0 = 3.3871327179;
		a1 = 50.434271938;
		a2 = 159.29113202;
		a3 = 59.109374720;
		b1 = 17.895169469;
		b2 = 78.775757664;
		b3 = 67.18756360;

		/* coefficients for p neither close to 1/2 nor 0 or 1 */
		c0 = 1.4234372777;
		c1 = 2.7568153900;
		c2 = 1.3067284816;
		c3 = .17023821103;
		d1 = .73700164250;
		d2 = .12021132975;

		/* coefficients for p near 0 or 1 */
		e0 = 6.6579051150;
		e1 = 3.0812263860;
		e2 = .42868294337;
		e3 = .017337203997;
		f1 = .24197894225;
		f2 = .012258202635;

		q = p - half;

		if (fabs(q) <= split1) {
			r = const1 - q * q;
			endval = q * (((a3 * r + a2) * r + a1) * r + a0)
									/ (((b3 * r + b2) * r + b1) * r + one);
			observ[i] = mu[i] + stdev[i] * endval;
			continue;
		}

		if (q < 0.0)
			r = p;
		else
			r = one - p;

		if (r <= zero)
			return 0;

		r = sqrt(-log(r));

		if (r <= split2) {
			r = r - const2;
			endval = (((c3 * r + c2) * r + c1) * r + c0)
									/ ((d2 * r + d1) * r + one);
			observ[i] = endval;
		}
		else {
			r = r - split2;
			endval = (((e3 * r + e2) * r + e1) * r + e0)
									/ ((f2 * r + f1) * r + one);
			observ[i] = endval;
		}
		if (q < 0)
			observ[i] = -1 * observ[i];
		observ[i] = mu[i] + stdev[i] * observ[i];
	}

	return (1);
}//nnormal()

int weibull(double scaleParam, double shapeParam, int numOmega, vector observ, long long *seed) {
	int 	n;
	double 	u;

	for (n = 0; n < numOmega; n++) {
		u = randUniform(seed);
		observ[n] = scaleParam*pow(log(u), 1/shapeParam);
	}

	return 0;
}//END weibull()

float scalit(float lower, float upper, long long *seed) {
	float val, wide;

	wide = upper - lower;
	val = randUniform(seed);

	return ((wide * val) + lower);
}//END scalit()

float randUniform(long long *SEED) {
	static int lo_bits, hi_bits;

	lo_bits = ((*SEED) & 0xFFFFL) * 16807;
	hi_bits = (int) (((*SEED) >> 16) * 16807) + (lo_bits >> 16);
	*SEED = ((lo_bits & 0xFFFFL) - 0x7FFFFFFFL) + ((hi_bits & 0x7FFFL) << 16) + (hi_bits >> 15);

	return ((*SEED) < 0 ? ((*SEED) += 0x7FFFFFFFL) : (*SEED)) * 4.656612875E-10;
}//END randUniform()

float randUniform_new(long long *SEED) {
	static int lo_bits, hi_bits;

	lo_bits = ((*SEED) & 0xFFFFL) * 16807;
	hi_bits = (int) (((*SEED) >> 16) * 16807) + (lo_bits >> 16);
	*SEED = ((lo_bits & 0xFFFFL) - 0x7FFFFFFFL) + ((hi_bits & 0x7FFFL) << 16) + (hi_bits >> 15);

	srand((unsigned) (*SEED));
	return (rand()/RAND_MAX+1);
}

int randInteger(long long *SEED, int iMax) {
	static int lo_bits, hi_bits;
	int val;

	lo_bits = ((*SEED) & 0xFFFFL) * 16807;
	hi_bits = (int) (((*SEED) >> 16) * 16807) + (lo_bits >> 16);
	*SEED = ((lo_bits & 0xFFFFL) - 0x7FFFFFFFL) + ((hi_bits & 0x7FFFL) << 16) + (hi_bits >> 15);

	srand((unsigned) (*SEED));
	val = rand() % iMax;

	return val;
}//END randInteger()

/* This function uses a sampling technique to set up a sample average approximation problem. The sampling procedure is conducted according to the continuous distribution and parameters provided in
 * stocType. The function takes number of samples as an input from the user. The function outputs the simulated observations as a matrix with each row corresponding to a random variable, and column corresponds to
 * a simulated observation. */
int setupSAA(stocType *stoc, long long *seed, vector **simObservVals, vector *probs, int *numSamples, double TOLERANCE) {
	int 	obs, n, offset, k;

	if ( (*numSamples) == 0 ) {
		/* number of samples in SAA */
		printf("Enter the number of samples used for setting up the SAA : ");
		scanf("%d", numSamples);
	}
	printf("Generating SAA with %d samples.\n", (*numSamples));

	(*simObservVals) = (vector *) mem_realloc((*simObservVals), (*numSamples)*sizeof(vector));
	(*probs) = (vector) mem_realloc((*probs), (*numSamples)*sizeof(double));

	if ( strstr(stoc->type, "BLOCKS_DISCRETE") != NULL ) {
		for (obs = 0; obs < (*numSamples); obs++ ) {
			(*simObservVals)[obs] = (vector) arr_alloc(stoc->numOmega+1, double);
           
//#ifdef STOCH_CHECK
//            for ( k = 0; k < stoc->numOmega+1; k++ )
//                printf("%4.6lf, ", (*simObservVals)[obs][k]);
//            printf("\n");
//#endif
            offset = 1;
            for ( n = 0; n < stoc->numGroups; n++ ) {
                generateBlocks(stoc, (*simObservVals)[obs]+offset, n, seed);
                offset += stoc->numPerGroup[n];
            }
#ifdef STOCH_CHECK
            printf("%d-th sample\n", obs);
            for ( k = 1; k <= stoc->numOmega; k++ )
                printf("%4.6lf, ", (*simObservVals)[obs][k]);
            printf("\n");
#endif
            (*probs)[obs] = 1.0/(double) (*numSamples);
		}
	}
	else if ( !strcmp(stoc->type, "INDEP_DISCRETE")) {
		for (obs = 0; obs < (*numSamples); obs++ ) {
			(*simObservVals)[obs] = (vector) arr_alloc(stoc->numOmega+1, double);
			generateIndep(stoc, (*simObservVals)[obs]+1, 0, seed);
			(*probs)[obs] = 1.0/(double) (*numSamples);
		}
	}
	else if ( !strcmp(stoc->type, "INDEP_NORMAL") ) {
		for (obs = 0; obs < (*numSamples); obs++ ) {
			(*simObservVals)[obs] = (vector) arr_alloc(stoc->numOmega+1, double);
			nnormal(stoc->mean, stoc->vals[0], stoc->numOmega, (*simObservVals)[obs]+1, seed);
			(*probs)[obs] = 1.0/(double) (*numSamples);
		}
	}
	else if ( !strcmp(stoc->type, "LINTRAN") ) {
		for (obs = 0; obs < (*numSamples); obs++ ) {
			(*simObservVals)[obs] = (vector) arr_alloc(stoc->numOmega+1, double);
			generateLinTran(stoc, (*simObservVals)[obs]+1, 0, TOLERANCE, seed);
			(*probs)[obs] = 1.0/(double) (*numSamples);
		}
	}
	else {
		errMsg("sampling", "setupSAA", "no procedure for simulating distribution type", 0);
		return 1;
	}

	return 0;
}//END setupSAA
