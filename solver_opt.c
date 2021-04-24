/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"

#define SAFE_ASSERT(condition, message) while(condition) {	\
	perror(message);		\
	goto failure;			\
}

double* my_solver(int N, double *A, double* B) {
	/* Base matrices aren't stored in registers since they are rarely referenced */
	double *C = NULL, *At = NULL, *BBt = NULL, *ABBt = NULL;

	/* Declaring indices as registers decreased the execution time */
	register int i, j, k;

	/* Partial sum */
	register double sum;

	/*
	 * all matrix operations below follow the next conventions:
	 *     O(utput) = L(eft) x R(ight)
	 * or: O(utput) = L'(eft) - for transposing
	 *
	 * register variables pO, pL and pR point the next element that will
	 * be used for the operation
	 *
	 * register variables pBaseX are helper base pointers, useful for
	 * pointer arithmetic
	 */
	register double *pO, *pL, *pR;
	register double *pBaseO, *pBaseL, *pBaseR;


	/* Matrices that need to be zeroed are alloc'd with calloc */
	C = calloc(sizeof(double), N * N);
	SAFE_ASSERT(C == NULL, "Failed calloc: C");

	At = calloc(sizeof(double), N * N);
	SAFE_ASSERT(At == NULL, "Failed calloc: At");

	BBt = malloc(sizeof(double) * N * N);
	SAFE_ASSERT(BBt == NULL, "Failed malloc: BBt");

	ABBt = calloc(sizeof(double), N * N);
	SAFE_ASSERT(ABBt == NULL, "Failed calloc: ABBt");


	/*
	 * Calculation: At = A'
	 * Why: later on I have to compute C = A' x A.
	 *      if done naively, this compuation roughly translates to:
	 *        C[i][j] = A[k][i] * A[k][j] (as seen in the unoptimized solution)
	 *      obviously this is bad because both accesses to A are non-sequential
	 * Solution: precompute A' and make better use of cache and the form of A'
	 *
	 * Optimizations for At = A':
	 *   - skip lower part of A since A is upper-triangular
	 *   - pointer arithmetic
	 *   - semi memory access optimization:
	 *       - access to A(L) is sequential;
	 *       - access to At(O) is non-sequential (N-by-N step)
	 *           (that's how the basic transpose algorithm works)
	 *   - inequality (!=) instead of less (<) comparisons
	 */

	pBaseO = At;
	pBaseL = A;
	for (i = 0; i != N; ++i) {
		pO = pBaseO;
		pL = pBaseL;
		for (j = i; j != N; ++j) {
			*pO = *pL;
			pL++;
			pO += N;
		}
		pBaseO += N + 1;
		pBaseL += N + 1;
	}

	/*
	 * Calculation: C = A' x A
	 * Naive: C[i][j] = At[i][k] * A[k][j]
	 *
	 * Two main mutually-exclusive optimizations can be applied here:
	 * 1. it is known that At = lower-triangular and A = upper-triangular
	 *    so we could skip a lot of scalar multiplications by limiting k range:
	 *      k = [0 .. min(i, j)]
	 *    however this means that k depends on both i and j, so swapping for's
	 *    order for better cache optimizations is not possible anymore
	 * 2. optimize memory accesses, however there will be some redundant products
	 *     (after some tests, I sticked with this approach)
	 *
	 * Optimizations:
	 *   - skip upper part of At since At is lower-triangular
	 *   - pointer arithmetic
	 *   - full memory access optimization:
	 *       - access to At(L) is constant;
	 *       - access to A(R) and C(O) is sequential
	 *   - inequality (!=) instead of less (<) comparisons
	 */

	pBaseO = C;
	pBaseL = At;
	pBaseR = A;
	for (i = 0; i != N; ++i) {
		pL = pBaseL;
		pR = pBaseR;
		for (k = 0; k <= i; ++k) {
			pO = pBaseO;
			for (j = 0; j != N; ++j) {
				*pO += *pL * *pR;
				pO++;
				pR++;
			}
			pL++;
		}
		pBaseO += N;
		pBaseL += N;
	}

	/*
	 * Calculation: BBt = B x B'
	 * Naive: BBt[i][j] = B[i][k] * B[j][k] (naive approach is pretty good already :) )
	 *
	 * Optimizations:
	 *   - pointer arithmetic
	 *   - full memory access optimization:
	 *       - access to BBt(O) is constant;
	 *       - access to B(L) and B(R) is sequential
	 *   - inequality (!=) instead of less (<) comparisons
	 */

	pO = BBt;
	pBaseL = B;
	pBaseR = B;
	for (i = 0; i != N; ++i) {
		pR = pBaseR;
		for (j = 0; j != N; ++j) {
			pL = pBaseL;
			sum = 0.0;
			for (k = 0; k != N; ++k) {
				sum += *pL * *pR;
				pL++;
				pR++;
			}
			*pO = sum;
			pO++;
		}
		pBaseL += N;
	}

	/*
	 * Calculation: ABBt = A x (B x B')
	 * Naive: ABBt[i][j] = A[i][k] * BBt[k][j]
	 *
	 * Optimizations:
	 *   - skip lower part of A since A is upper-triangular
	 *   - pointer arithmetic
	 *   - full memory access optimization:
	 *       - access to A(L) is constant;
	 *       - access to BBt(R) and ABBt(O) is sequential
	 *   - inequality (!=) instead of less (<) comparisons
	 */

	pBaseL = A;
	pBaseR = BBt;
	pBaseO = ABBt;
	for (i = 0; i != N; ++i) {
		pL = pBaseL;
		pR = pBaseR;
		for (k = i; k != N; ++k) {
			pO = pBaseO;
			for (j = 0; j != N; ++j) {
				*pO += *pL * *pR;
				pR++;
				pO++;
			}
			pL++;
		}
		pBaseL += N + 1;
		pBaseR += N;
		pBaseO += N;
	}

	/*
	 * Calculation: C = [A x (B x B')] + (A' x A)
	 * Code: C[i][j] += ABBt[i][j]   (C = A' x A already)
	 *
	 * Optimizations:
	 *   - skip lower part of A since A is upper-triangular
	 *   - pointer arithmetic
	 *   - full memory access optimization:
	 *       - access to C(O) and ABBt(L) is sequential
	 *   - inequality (!=) instead of less (<) comparisons
	 */

	pO = C;
	pL = ABBt;
	j = N * N;
	for (i = 0; i != j; ++i) {
		*pO++ += *pL++;
	}

	goto cleanup;

failure:
	free(C);
	C = NULL;

cleanup:
	free(At);
	free(BBt);
	free(ABBt);
	
	return C;
}
