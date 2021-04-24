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
	register double *C = NULL, *At = NULL, *BBt = NULL, *ABBt = NULL;
	register int i, j, k;

	C = calloc(sizeof(double), N * N);
	SAFE_ASSERT(C == NULL, "Failed calloc: C");

	At = calloc(sizeof(double), N * N);
	SAFE_ASSERT(At == NULL, "Failed calloc: At");

	BBt = malloc(sizeof(double) * N * N);
	SAFE_ASSERT(BBt == NULL, "Failed calloc: BBt");

	ABBt = calloc(sizeof(double), N * N);
	SAFE_ASSERT(ABBt == NULL, "Failed calloc: ABBt");

	register double *orig_pa;
	register double *orig_pb;
	register double *orig_pc;
	register double *pa;
	register double *pb;
	register double *pc;

	orig_pa = A;
	orig_pc = At;
	for (i = 0; i != N; ++i) {
		pa = orig_pa;
		pc = orig_pc;
		for (j = i; j != N; ++j) {
			*pc = *pa;
			pa++;
			pc += N;
		}
		orig_pa += N + 1;
		orig_pc += N + 1;
	}

	/* C = A' x A */
	orig_pc = C;
	orig_pa = At;
	for (i = 0; i != N; ++i) {
		pa = orig_pa;
		pb = A;
		for (k = 0; k <= i; ++k) {
			pc = orig_pc;
			for (j = 0; j != N; ++j) {
				*pc += *pa * *pb;
				pc++;
				pb++;
			}
			pa++;
		}
		orig_pc += N;
		orig_pa += N;
	}


	/* BBt = B x B' -- ok */
	pc = BBt;
	orig_pa = B;
	for (i = 0; i != N; ++i) {
		pb = &B[0];
		for (j = 0; j != N; ++j) {
			pa = orig_pa;
			register double sum = 0.0;
			for (k = 0; k != N; ++k) {
				sum += *pa * *pb;
				pa++;
				pb++;
			}
			*pc = sum;
			pc++;
		}
		orig_pa += N;
	}

	/* ABBt = A x (B x B') */
	orig_pa = A;
	orig_pb = BBt;
	orig_pc = ABBt;
	for (i = 0; i != N; ++i) {
		pa = orig_pa;
		pb = orig_pb;
		for (k = i; k != N; ++k) {
			pc = orig_pc;
			for (j = 0; j != N; ++j) {
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
	for (i = 0; i != limit; ++i) {
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
