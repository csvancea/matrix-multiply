/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"

#include <string.h>
#include <cblas.h>

#define SAFE_ASSERT(condition, message) while(condition) {	\
	perror(message);		\
	goto failure;			\
}

double* my_solver(int N, double *A, double *B) {
	double *C = NULL, *AB = NULL;

	C = malloc(sizeof(double) * N * N);
	SAFE_ASSERT(C == NULL, "Failed malloc: C");

	AB = malloc(sizeof(double) * N * N);
	SAFE_ASSERT(AB == NULL, "Failed malloc: AB");

	/* C = A' x A */
	memcpy(C, A, sizeof(double) * N * N);
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
		N, N,
		1.0, A, N,
		C, N
	);

	/* AB = A x B */
	memcpy(AB, B, sizeof(double) * N * N);
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
		N, N,
		1.0, A, N,
		AB, N
	);

	/* C = [(A x B) x B'] + (A' x A) */
	/* C =    AB    x B'  +     C    */
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
		N, N, N,
		1.0, AB, N,
		B, N,
		1.0, C, N
	);

	goto cleanup;

failure:
	free(C);
	C = NULL;

cleanup:
	free(AB);
	
	return C;
}
