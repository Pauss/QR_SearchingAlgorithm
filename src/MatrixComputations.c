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
/*Function to convert double matrix to gsl_matrix type*/
/*=========================================*/
void convert_to_gsl(T_FILE_DIM* file_dim, gsl_matrix *M) {

	for (uint8 i = 0; i < file_dim->lines; i++) {
		for (uint8 j = 0; j < file_dim->columns; j++) {
			gsl_matrix_set(M, i, j, file_dim->matrix[i][j]);
		}

	}

}

/*=========================================*/
/*Function to print a gsl_matrix*/
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
/*Function to compute product of a gsl_matrix with a gsl_vector*/
/*=========================================*/
boolean product_matrix_vector(gsl_matrix* M, gsl_vector* v, gsl_vector* result) {

	/*if matrix's number of columns is not equal with vector size then product can't be done*/
	if (M->size2 != v->size) {
		return FALSE;

	} else {
		for (uint8 i = 0; i < M->size1; i++) {
			double S = 0;
			for (uint8 j = 0; j < M->size1; j++) {
				S += gsl_matrix_get(M, i, j) * gsl_vector_get(v, j);
			}
			gsl_vector_set(result, i, S);
		}
	}

	return TRUE;

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

		double temp_i_i = gsl_matrix_get(M_local, i, i);

		for (j = 0; j < multiplied_size; j++) {

			/*getting value 1 for [I] matrix*/
			gsl_matrix_set(M_local, i, j,
					(gsl_matrix_get(M_local, i, j) / temp_i_i));

		}

	}


	gsl_matrix_set_zero(M);

	/*extract [Inverse] matrix from local matrix*/
	view = gsl_matrix_submatrix(M_local, INIT, sqaured_size, sqaured_size,
			sqaured_size);

	/*set to the input matrix inverse matrix computed*/
	gsl_matrix_memcpy(M, &view.matrix);

	return TRUE;

}

/*=========================================*/
/*This function compute Euclidean norm*/
/*=========================================*/
double euclidean_norm(gsl_vector* V, uint8 i) {

	double sum = INIT;

	for (; i < V->size; i++) {
		sum += pow((gsl_vector_get(V, i)), 2);

	}

	return sqrt(sum);

}


