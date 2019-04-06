/*=========================================*/
/*
 * MatrixComputations.c
 *
 *  Created on: Nov 3, 2018
 *      Author: Paus
 */
/*=========================================*/
/*include*/
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
	for (uint16 i = 0; i < M->size1; i++) {
		for (uint16 j = 0; j < M->size2; j++) {
			printf("m(%d,%d) = %f\t", i, j, gsl_matrix_get(M, i, j));
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
	for (uint16 i = 0; i < V->size; i++) {

			printf("v(%d) = %f ", i, gsl_vector_get(V, i));
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
		for (uint16 i = 0; i < M->size1; i++) {
			double S = 0;
			for (uint16 j = 0; j < M->size2; j++) {

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
		for (uint16 i = 0; i < M2->size2; i++) {
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
/*This function compute Euclidean norm*/
/*=========================================*/
double euclidean_norm(gsl_vector* V) {

	double sum = INIT;

	for (uint16 i = INIT; i < V->size; i++) {
		sum += pow((gsl_vector_get(V, i)), SCALE_2);

	}

	//return sqrt(sum);
	return (sum);

}

/*=========================================*/
/*This function add a sub-matrix to a matrix*/
/*=========================================*/
void add_submatrix(gsl_matrix* R, gsl_matrix* sub_matrix, uint8 index1, uint8 index2)
{
	double l_element;

	for (uint16 i = 0; i < sub_matrix->size1; i++) {
		for (uint16 j = 0; j < sub_matrix->size2; j++) {

			l_element = gsl_matrix_get(sub_matrix, i, j);

			gsl_matrix_set(R, i + index1, j + index2, l_element);
		}
	}

}

/*=========================================*/
/*This function delete Intercept column from the model after it adds Intercept to all columns*/
/*=========================================*/
gsl_matrix* add_intercept(gsl_matrix* R)
{

	gsl_vector* intercept_column = gsl_vector_alloc(R->size1);

    gsl_vector* temp_column = gsl_vector_alloc(R->size1);

	gsl_matrix* temp_matix = gsl_matrix_alloc(R->size1, R->size2+1);

	gsl_vector_set_all(intercept_column, INTERCEPT_VALUE);

	gsl_matrix_set_col(temp_matix, 0, intercept_column);

	for( uint16 i = 1; i < temp_matix->size2; i++){

		gsl_matrix_get_col(temp_column, R, i-1);

		gsl_matrix_set_col(temp_matix, i, temp_column );
	}

	return temp_matix;

}

/*=========================================*/
/*This function delete a column from a given matrix*/
/*=========================================*/
void delete_column(gsl_matrix* R, uint16 col) {

	gsl_matrix* l_matrix = gsl_matrix_alloc(R->size1, R->size2 - 1);
	gsl_vector* l_vector = gsl_vector_alloc(R->size1);

	if (col == R->size2) {

		/*only change the size of R*/
		R->size2--;

	} else {
		if (col > 0)

		{

			for (uint16 i = 0; i < col; i++) {
				gsl_matrix_get_col(l_vector, R, i);
				gsl_matrix_set_col(l_matrix, i, l_vector);
			}

		}

		for (uint16 i = col + 1; i < R->size2; i++) {
			gsl_matrix_get_col(l_vector, R, i);
			gsl_matrix_set_col(l_matrix, i - 1, l_vector);
		}

		R->size2--;
		gsl_matrix_memcpy(R, l_matrix);
	}

	gsl_matrix_free(l_matrix);
	gsl_vector_free(l_vector);

}
