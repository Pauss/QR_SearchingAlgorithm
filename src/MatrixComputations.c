/*=========================================*/
/*
 * MatrixComputations.c
 *
 *  Created on: Nov 3, 2018
 *      Author: Paus
 */
/*=========================================*/
/*include*/
#include <stdio.h>
#include <math.h>
#include "MatrixComputations.h"
/*=========================================*/
/*private data*/
/*=========================================*/
/*global data*/
/*=========================================*/
/*private functions*/
/*=========================================*/

/*=========================================*/
/*Description*/
/*=========================================*/
void convert_to_gsl(T_FILE_DIM* file_dim, gsl_matrix *M) {

	for (uint8 i = 0; i < file_dim->lines; i++) {
		for (uint8 j = 0; j < file_dim->columns; j++) {
			gsl_matrix_set(M, i, j, file_dim->matrix[i][j]);
		}

	}

}

/*=========================================*/
/*Description*/
/*=========================================*/
void print_matrix(gsl_matrix *M) {

	for (uint8 i = 0; i < M->size1; i++) {
		for (uint8 j = 0; j < M->size2; j++) {
			printf("m(%d,%d) = %lf\t", i, j, gsl_matrix_get(M, i, j));
		}
		printf("\n");
	}
	printf("\n");
}

/*=========================================*/
/*This function compute inverse of input matrix*/
/*=========================================*/
boolean compute_matrix_inverse(gsl_matrix * M) {

	double ratio;
	uint8 i, j, k;
	uint8 sqaured_size, multiplied_size;
	gsl_matrix *M_local;
	gsl_matrix_view view;

	/*double matrix[M->size1*2][M->size2*2];*/

	/*In order to compute inverse, matrix should be n*n sized*/
	if (M->size1 != M->size2) {
		return FALSE;
	}

	/* In order to compute determinant value we can use QR decomposition
	 * A = Q * R
	 * determinant(A) = determinant(Q) * determinant(R)
	 * Q is an orthogonal matrix so determinant(Q) is whether 1 or -1
	 * determinant(A) = product(Rii)*/
	/*TODO*/

	/* Check this algorithm works also for [M] sized 15*15*/
	/*TODO*/

	/* set sqaured_size to be the size of input matrix*/
	sqaured_size = M->size1;
	/* set multiplied_size to be the size of input matrix multiplied to number 2 */
	multiplied_size = M->size1 * SCALE_2;

	/*set for the local matrix 2*size of input matrix*/
	M_local = gsl_matrix_alloc(multiplied_size, multiplied_size);

	/* local matrix used for computations*/
	M_local = gsl_matrix_alloc(multiplied_size, multiplied_size);

	for (i = 0; i < sqaured_size; i++) {

		for (j = 0; j < sqaured_size; j++) {

			/* copy input matrix data into local matrix*/
			gsl_matrix_set(M_local, i, j, gsl_matrix_get(M, i, j));

		}

		/* Attach to the local matrix identity matrix
		 * [local_matrix][0] changes to [local_matrix][I]*/
		gsl_matrix_set(M_local, i, (i + sqaured_size), IDENTITY_VALUE);

	}

	for (i = 0; i < sqaured_size; i++) {

		for (j = 0; j < sqaured_size; j++) {

			/*local_matrix[I] changes to [I][M inverse]*/
			if (i != j) {

				/*ratio is scaling factor for getting 0 for [I] matrix*/
				ratio = gsl_matrix_get(M_local, j, i)
						/ gsl_matrix_get(M_local, i, i);

				for (k = 0; k < multiplied_size; k++) {

					gsl_matrix_set(M_local, j, k,
							(gsl_matrix_get(M_local, j, k)
									- (ratio * gsl_matrix_get(M_local, i, k))));

				}

			}

		}

	}

	for (i = 0; i < sqaured_size; i++) {

		/*make a copy of element (i,i)*/
		double temp_i_i = gsl_matrix_get(M_local, i, i);

		for (j = 0; j < multiplied_size; j++) {

			/*getting value 1 for [I] matrix*/
			gsl_matrix_set(M_local, i, j,
					(gsl_matrix_get(M_local, i, j) / temp_i_i));

		}

	}


	gsl_matrix_set_zero(M);

	/*upper-left elements[0,squared_size], number of lines and number of columns of the view*/
	view = gsl_matrix_submatrix(M_local, INIT, sqaured_size, sqaured_size,
			sqaured_size);

	gsl_matrix_memcpy(M, &view.matrix);

	return TRUE;

}

