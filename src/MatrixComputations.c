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
