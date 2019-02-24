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

	printf("\n");
	for (uint8 i = 0; i < M->size1; i++) {
		for (uint8 j = 0; j < M->size2; j++) {
			printf("m(%d,%d) = %lf\t", i, j, gsl_matrix_get(M, i, j));
		}
		printf("\n");
	}
	printf("\n");
}

/*=========================================*/
/*Function to print a gsl_vector*/
/*=========================================*/
void print_vector (gsl_vector *V) {

	printf("\n");
	for (uint8 i = 0; i < V->size; i++) {

			printf("m(%d) = %lf\t", i, gsl_vector_get(V, i));
	}
	printf("\n");
}

/*=========================================*/
/*Function to compute product of a gsl_matrix with a gsl_vector*/
/*=========================================*/
boolean product_matrix_vector(gsl_matrix* M, gsl_vector* v, gsl_vector* result) {

	//if matrix's number of columns is not equal with vector size then product can't be done
	if ((M->size2 != v->size) || (M->size1 != result->size)) {
		return FALSE;

	} else {
		for (uint8 i = 0; i < M->size1; i++) {
			double S = 0;
			for (uint8 j = 0; j < M->size2; j++) {

				S += (double) gsl_matrix_get(M, i, j) * gsl_vector_get(v, j);
			}

			gsl_vector_set(result, i, S);

		}
	}

	return TRUE;

}

/*=========================================*/
/*Function to compute product of a gsl_matrix with a gsl_matrix*/
/*=========================================*/
boolean product_matrix(gsl_matrix* M, gsl_matrix* M2, gsl_matrix* result) {

	gsl_vector* product = gsl_vector_alloc(M->size1);
	gsl_vector* column = gsl_vector_alloc(M->size1);

	//if matrix's number of columns is not equal with vector size then product can't be done
	if (((M->size2 != M2->size1) || (M->size1 != result->size1))
			|| (M2->size2 != result->size2)) {
		return FALSE;

	} else {
		for (uint8 i = 0; i < M2->size2; i++) {
			gsl_matrix_get_col(column, M2, i);

			product_matrix_vector(M, column, product);

			gsl_matrix_set_col(result, i, product);
		}
	}

	gsl_vector_free(product);
	gsl_vector_free(column);

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

	gsl_matrix_free(M_local);

	return TRUE;

}

/*=========================================*/
/*This function compute Euclidean norm*/
/*=========================================*/
double euclidean_norm(gsl_vector* V) {

	double sum = INIT;

	for (uint8 i = INIT; i < V->size; i++) {
		sum += pow((gsl_vector_get(V, i)), SCALE_2);

	}

	return sqrt(sum);
	//return (sum);

}

/*=========================================*/
/*This function add a sub-matrix to a matrix*/
/*=========================================*/
void add_submatrix(gsl_matrix* R, gsl_matrix* sub_matrix, uint8 index1, uint8 index2)
{
	double l_element;

	for (uint8 i = 0; i < sub_matrix->size1; i++) {
		for (uint8 j = 0; j < sub_matrix->size2; j++) {

			l_element = gsl_matrix_get(sub_matrix, i, j);

			gsl_matrix_set(R, i + index1, j + index2, l_element);
		}
	}

}

/*=========================================*/
/*This function delete a column from a given matrix*/
/*=========================================*/
void delete_column(gsl_matrix* R, uint8 col) {

	gsl_matrix* l_matrix = gsl_matrix_alloc(R->size1, R->size2 - 1);
	gsl_vector* l_vector = gsl_vector_alloc(R->size1);

	if (col == R->size2) {

		/*only change the size of R*/
		R->size2--;

	} else {
		if (col > 0)

		{

			for (uint8 i = 0; i < col; i++) {
				gsl_matrix_get_col(l_vector, R, i);
				gsl_matrix_set_col(l_matrix, i, l_vector);
			}

		}

		for (uint8 i = col + 1; i < R->size2; i++) {
			gsl_matrix_get_col(l_vector, R, i);
			gsl_matrix_set_col(l_matrix, i - 1, l_vector);
		}

		R->size2--;
		gsl_matrix_memcpy(R, l_matrix);
	}

	gsl_matrix_free(l_matrix);
	gsl_vector_free(l_vector);

}
