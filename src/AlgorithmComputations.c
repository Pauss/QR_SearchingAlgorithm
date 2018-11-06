/*
 * AlgorithmComputations.c
 *
 *  Created on: Nov 4, 2018
 *      Author: Paus
 */

/*=========================================*/
/*include*/
#include "AlgorithmComputations.h"
/*=========================================*/
/*private data*/
static gsl_matrix* A;
static gsl_vector* y;
/*=========================================*/
/*global data*/
/*=========================================*/
/*private functions*/
/*=========================================*/

/*=========================================*/
/*Get 'A' in order to solve y=A*x */
/*=========================================*/
gsl_matrix* get_A_matrix()
{
	return A;
}

/*=========================================*/
/*Get 'y' vector in order to solve y=A*x */
/*=========================================*/
gsl_vector* get_y_vector()
{
	return y;
}

/*=========================================*/
/*Set 'A' matrix from input data in order to solve y=A*x */
/*=========================================*/
void set_A_matrix(T_FILE_DIM* file_dim)
{

	gsl_matrix* matrix_file = gsl_matrix_alloc(file_dim->lines, file_dim->columns);
	convert_to_gsl(file_dim, matrix_file);

	/*extract [A] matrix from input*/
	gsl_matrix_view A_view = gsl_matrix_submatrix(matrix_file, INIT, 1, file_dim->lines,
			file_dim->columns-1);

	/*set A matrix*/
	A = gsl_matrix_alloc(file_dim->lines, file_dim->columns -1);
	gsl_matrix_memcpy(A, &A_view.matrix);
}

/*=========================================*/
/*Set 'y' vector from input data in order to solve y=A*x */
/*=========================================*/
void set_y_vector(T_FILE_DIM* file_dim)
{
	gsl_matrix* matrix_file = gsl_matrix_alloc(file_dim->lines, file_dim->columns);
	convert_to_gsl(file_dim, matrix_file);

	/*extract [y] vector from input*/
	gsl_vector_view y_view = gsl_matrix_column(matrix_file, INIT);

	/*set A matrix*/
	y = gsl_vector_alloc(file_dim->lines);
	gsl_vector_memcpy(y, &y_view.vector);
}

/*=========================================*/
/*	/* 1.QR decomposition
	   2.unpack Q and R
	   3.Q transpose * y => y conjugate
	   4.compute RSS */
/*=========================================*/
/*TODO, change name or separate them, Q,R, to can use them in other functions*/
void QR_decomposition(void)
{

	/*QR used to decompose matrix A such that A = Q*R*/
	gsl_matrix* QR = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix_memcpy(QR, A);

	/*RSS = residual sum of squares*/
	double RSS = 0;

	gsl_vector* tau = gsl_vector_alloc(MIN_VALUE(A->size1, A->size2));
	gsl_matrix *Q = gsl_matrix_alloc(A->size1, A->size1);
	gsl_matrix *R = gsl_matrix_alloc(A->size1, A->size2);

	/*vector of approximated solution*/
	gsl_vector * x = gsl_vector_alloc(A->size2);
	/*vector of residual values*/
	gsl_vector * E = gsl_vector_alloc(A->size1);

	gsl_linalg_QR_decomp(QR, tau);

	/*This function unpacks the encoded QR decomposition (QR,tau) into the matrices Q and R, where Q is M-by-M and R is M-by-N. */
	gsl_linalg_QR_unpack(QR, tau, Q, R);

	/* This function finds the least squares solution to the over-determined system A * x = y,
	 * where the matrix A has more rows than columns.
	 * The least squares solution minimizes the Euclidean norm of the residual, ||Ax - b||.*/
	gsl_linalg_QR_lssolve(QR, tau, y, x, E);

	RSS = euclidean_norm(x);

	print_matrix(A);
	print_vector(y);

	print_matrix(QR);
	print_vector(tau);

	print_matrix(Q);
	print_matrix(R);

	print_vector(x);
	printf("RSS: %lf\n", RSS);

}



/*=========================================*/
/*Description*/
/*=========================================*/

/*=========================================*/
/*Description*/
/*=========================================*/

/*=========================================*/
/*Description*/
/*=========================================*/

/*=========================================*/
/*Description*/
/*=========================================*/
