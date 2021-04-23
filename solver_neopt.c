/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"

#define SAFE_ASSERT(condition, message) while(condition) {	\
	perror(message);		\
	goto failure;			\
}

static inline int min(int a, int b) {
	return a < b ? a : b;
}

double* my_solver(int N, double *A, double* B) {
	double *C = NULL, *AtA = NULL, *BBt = NULL, *ABBt = NULL;
	int i, j, k;

	C = calloc(sizeof(double), N * N);
	SAFE_ASSERT(C == NULL, "Failed calloc: C");

	AtA = calloc(sizeof(double), N * N);
	SAFE_ASSERT(AtA == NULL, "Failed calloc: AtA");

	BBt = calloc(sizeof(double), N * N);
	SAFE_ASSERT(BBt == NULL, "Failed calloc: BBt");

	ABBt = calloc(sizeof(double), N * N);
	SAFE_ASSERT(ABBt == NULL, "Failed calloc: ABBt");

	/* AtA = A' x A */
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			for (k = 0; k <= min(i, j); ++k) {
				AtA[i * N + j] += A[k * N + i] * A[k * N + j];
			}
		}
	}

	/* BBt = B x B' */
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			for (k = 0; k < N; ++k) {
				BBt[i * N + j] += B[i * N + k] * B[j * N + k];
			}
		}
	}

	/* ABBt = A x (B x B') */
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			for (k = i; k < N; ++k) {
				ABBt[i * N + j] += A[i * N + k] * BBt[k * N + j];
			}
		}
	}

	/* C = [A x (B x B')] + (A' x A) */
	for (i = 0; i < N * N; ++i) {
		C[i] = ABBt[i] + AtA[i];
	}

	goto cleanup;

failure:
	free(C);
	C = NULL;

cleanup:
	free(AtA);
	free(BBt);
	free(ABBt);
	
	return C;
}
