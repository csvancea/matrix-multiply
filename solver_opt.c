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
	double *C = NULL, *At = NULL, *BBt = NULL, *ABBt = NULL;
	register int i, j, k;

	C = malloc(sizeof(double) * N * N);
	SAFE_ASSERT(C == NULL, "Failed calloc: C");

	At = calloc(sizeof(double), N * N);
	SAFE_ASSERT(At == NULL, "Failed calloc: At");

	BBt = malloc(sizeof(double) * N * N);
	SAFE_ASSERT(BBt == NULL, "Failed calloc: BBt");

	ABBt = calloc(sizeof(double), N * N);
	SAFE_ASSERT(ABBt == NULL, "Failed calloc: ABBt");

	double *orig_pa;
	double *orig_pb;
	double *orig_pc;
	register double *pa;
	register double *pb;
	register double *pc;

	orig_pa = A;
	orig_pc = At;
	for (i = 0; i < N; ++i) {
		pa = orig_pa;
		pc = orig_pc;
		for (j = i; j < N; ++j) {
			*pc = *pa;
			pa++;
			pc += N;
		}
		orig_pa += N + 1;
		orig_pc += N + 1;
	}

	/* C = A' x A */
	/* sau nu mai tin cont de L-U si modific ordinea for-urilor */

	pc = C;
	orig_pb = At;
	for (i = 0; i < N; ++i) {
		orig_pa = A;
		for (j = 0; j < N; ++j) {
			register double sum = 0.0;
			register int limit;
			if (i < j)
				limit = i;
			else
				limit = j;

			pb = orig_pb;
			pa = orig_pa;

			for (k = 0; k <= limit; ++k) {
				sum += *pb * *pa;
				pb++;
				pa += N;
			}
			*pc = sum;
			pc++;
			orig_pa++;
		}
		orig_pb += N;
	}	

	/* BBt = B x B' -- ok */
	double *pbbt = BBt;
	orig_pa = B;
	for (i = 0; i < N; ++i) {
		pb = &B[0];
		for (j = 0; j < N; ++j) {
			pa = orig_pa;
			register double sum = 0.0;
			for (k = 0; k < N; ++k) {
				sum += *pa * *pb;
				pa++;
				pb++;
			}
			*pbbt = sum;
			pbbt++;
		}
		orig_pa += N;
	}

	/* ABBt = A x (B x B') */
	orig_pa = A;
	orig_pb = BBt;
	orig_pc = ABBt;
	for (i = 0; i < N; ++i) {
		pa = orig_pa;
		pb = orig_pb;
		for (k = i; k < N; ++k) {
			pc = orig_pc;
			for (j = 0; j < N; ++j) {
				*pc += *pa * *pb;
				pb++;
				pc++;
			}
			pa++;
		}
		orig_pa += N + 1;
		orig_pb += N;
		orig_pc += N;
	}

	

	/* C = [A x (B x B')] + (A' x A) */
	register int limit = N*N;
	pc = C;
	pa = ABBt;
	for (i = 0; i < limit; ++i) {
		*pc++ += *pa++;
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
