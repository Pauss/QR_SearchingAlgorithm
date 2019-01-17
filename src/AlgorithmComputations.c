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
gsl_matrix* get_A_matrix() {
	return A;
}

/*=========================================*/
/*Get 'y' vector in order to solve y=A*x */
/*=========================================*/
gsl_vector* get_y_vector() {
	return y;
}

/*=========================================*/
/*Set 'A' matrix from input data in order to solve y=A*x */
/*=========================================*/
void set_A_matrix(T_FILE_DIM* file_dim) {

	gsl_matrix* matrix_file = gsl_matrix_alloc(file_dim->lines,
			file_dim->columns);
	convert_to_gsl(file_dim, matrix_file);

	/*extract [A] matrix from input*/
	gsl_matrix_view A_view = gsl_matrix_submatrix(matrix_file, INIT, 1,
			file_dim->lines, file_dim->columns - 1);

	/*set A matrix*/
	A = gsl_matrix_alloc(file_dim->lines, file_dim->columns - 1);
	gsl_matrix_memcpy(A, &A_view.matrix);
}

/*=========================================*/
/*Set 'y' vector from input data in order to solve y=A*x */
/*=========================================*/
void set_y_vector(T_FILE_DIM* file_dim) {
	gsl_matrix* matrix_file = gsl_matrix_alloc(file_dim->lines,
			file_dim->columns);
	convert_to_gsl(file_dim, matrix_file);

	/*extract [y] vector from input*/
	gsl_vector_view y_view = gsl_matrix_column(matrix_file, INIT);

	/*set y vector*/
	y = gsl_vector_alloc(file_dim->lines);
	gsl_vector_memcpy(y, &y_view.vector);
}

/*=========================================*/
/* 1.QR decomposition
 * 2.unpack Q and R
 * 3.calculates approximate solution, residual and RSS*/
/*=========================================*/
void QR_decomposition(gsl_matrix* matrix_input,
		Model_QR_components* matrix_components) {

	matrix_components->model = gsl_matrix_alloc(matrix_input->size1,
			matrix_input->size2);
	gsl_matrix_memcpy(matrix_components->model, matrix_input);

	/*QR used to decompose matrix A such that A = Q*R*/
	matrix_components->QR = gsl_matrix_alloc(matrix_input->size1,
			matrix_input->size2);
	gsl_matrix_memcpy(matrix_components->QR, matrix_input);

	matrix_components->tau = gsl_vector_alloc(
			MIN_VALUE(matrix_components->QR->size1,
					matrix_components->QR->size2));
	matrix_components->Q = gsl_matrix_alloc(matrix_components->QR->size1,
			matrix_components->QR->size1);
	matrix_components->R = gsl_matrix_alloc(matrix_components->QR->size1,
			matrix_components->QR->size2);
	/*vector of approximated solution*/
	matrix_components->solution = gsl_vector_alloc(
			matrix_components->QR->size2);
	/*vector of residual values*/
	matrix_components->residual = gsl_vector_alloc(
			matrix_components->QR->size1);

	gsl_linalg_QR_decomp(matrix_components->QR, matrix_components->tau);

	/*This function unpacks the encoded QR decomposition (QR,tau) into the matrices Q and R, where Q is M-by-M and R is M-by-N.*/
	gsl_linalg_QR_unpack(matrix_components->QR, matrix_components->tau,
			matrix_components->Q, matrix_components->R);

	/*This function calculates approximate solution of equation and residual */
	gsl_linalg_QR_lssolve(matrix_components->QR, matrix_components->tau, y,
			matrix_components->solution, matrix_components->residual);

	matrix_components->RSS = euclidean_norm(matrix_components->residual);

}

/*=========================================*/
/* Function to compute compute RSS
 * shorter version, without interest of other components
 * used for compute RSS of all sub-models*/
/*=========================================*/
double RSS_compute(gsl_matrix* QR) {

	/*vector of approximated solution*/
	gsl_vector * x = gsl_vector_alloc(QR->size2);
	/*vector of residual values*/
	gsl_vector * E = gsl_vector_alloc(QR->size1);

	gsl_vector* tau = gsl_vector_alloc(MIN_VALUE(QR->size1, QR->size2));
	double RSS = INIT;

	gsl_linalg_QR_decomp(QR, tau);

	/* This function finds the least squares solution to the over-determined system A * x = y,
	 * where the matrix A has more rows than columns.
	 * The least squares solution minimizes the Euclidean norm of the residual, ||Ax - b||.*/
	gsl_linalg_QR_lssolve(QR, tau, y, x, E);

	/*RSS = euclidean_norm(x);*/

	RSS = euclidean_norm(E);

	return (RSS);
}

/*=========================================*/
/*This function compute all possible subsets of ModelMatrix*/
/*=========================================*/
void compute_transitions_QR(void) {
	/*ToDo*/
	/*Error handling*/

	remove("Output_Steps");

	/*computing all possible combinations of columns*/
	for (uint8 i = 1; i <= A->size2; i++) {

		gsl_combination * columns_transitions = gsl_combination_alloc(A->size2, i);

		/*set combinations with values from 0 to n-1*/
		gsl_combination_init_first(columns_transitions);

		/*do combination while reaching all possible subsets*/
		do {

			/*get matrix model*/
			gsl_matrix* my_model = sub_model_matrix(columns_transitions);
			/*print_steps(RSS_compute(my_model), columns_transitions);*/
			printf("\n%lf\n", RSS_compute(my_model));
			gsl_combination_fprintf(stdout, columns_transitions,"%u ");


		} while (GSL_SUCCESS == gsl_combination_next(columns_transitions));
	}
}

/*=========================================*/
/*This function form a sub-model matrix*/
/*=========================================*/
gsl_matrix* sub_model_matrix(gsl_combination* matrix_combination) {

	/*make a copy of full model matrix*/
	gsl_matrix* matrix_transitions = gsl_matrix_alloc(A->size1,
			matrix_combination->k);
	gsl_matrix_set_zero(matrix_transitions);

	for (uint8 i = 0; i < matrix_transitions->size2; i++) {

		/* take each column index from full model matrix that is present in combination
		 * and add it to sub-model*/
		gsl_combination_get(matrix_combination, i);
		gsl_vector * v = gsl_vector_alloc(A->size1);
		gsl_matrix_get_col(v, A, gsl_combination_get(matrix_combination, i));
		gsl_matrix_set_col(matrix_transitions, i, v);

	}

	return matrix_transitions;

}

/*=========================================*/
/*Description*/
/*TODO*/
/*add a function to choose for new different scores (RSS)*/
/*=========================================*/

/*=========================================*/
/*This function calculates R* by decomposing (A|Y)*/
/*=========================================*/
void get_base_R(gsl_matrix* R, gsl_vector* x, gsl_matrix* R_result)
{

	gsl_vector* C;
	C = gsl_vector_alloc(R->size1);

    product_matrix_vector(R, x, C);

    gsl_vector_view column;

    for(uint8 i=0; i< (R->size2); i++)
    {
    	column = gsl_matrix_column(R, i);
    	gsl_matrix_set_col(R_result, i, &column.vector);
    }

    gsl_matrix_set_col(R_result, (R_result->size2)-1, C);

    gsl_vector_free(C);

}

/*=========================================*/
/*This function compute RSS using elements of C*/
/*=========================================*/
void efficient_RSS(gsl_vector* C, gsl_vector* RSS_models, const double RSS)
{

	uint8 n = RSS_models->size;
	double result = RSS;

	for(uint8 i=0; i < n; i++)
	{
		result += pow(gsl_vector_get(C,(n-i)),2);
		gsl_vector_set(RSS_models,i, result);
	}

}

/*=========================================*/
/*This function add a sub-matrix to a matrix*/
/*=========================================*/
void add_submatrix(gsl_matrix* R, gsl_matrix* sub_matrix, uint8 index1, uint8 index2)
{

	for(uint8 i = 0; i< sub_matrix->size1; i++)
	{
		for(uint8 j = 0; j<sub_matrix->size2; j++ )
		{
			gsl_matrix_set(R, i+index1, j+index2, gsl_matrix_get(sub_matrix, i, j));
		}
	}

}

/*=========================================*/
/*This function performs re - triangularization after switching 2 columns in R*/
/*=========================================*/

void retriangularization_R(gsl_matrix* R, uint8 column1, uint8 column2)
{
	double a, b, s, c, result;
	gsl_matrix* rotation = gsl_matrix_alloc(NR_ELEMENTS, NR_ELEMENTS);

	gsl_matrix_swap_columns(R, column1, column2);

	a = gsl_matrix_get(R, column1, column1);
	b = gsl_matrix_get(R, column2, column1);

	gsl_linalg_givens(a, b, &c, &s);

	c = ABS_VALUE(c);
	s = ABS_VALUE(s);

	/*todo*/
	/*SIGN of R and case when r = 0*/

	gsl_matrix_set_all(rotation, c);
	gsl_matrix_set(rotation, 0, 1, s);
	gsl_matrix_set(rotation, 1, 0, -s);

	result = sqrt(pow(a,2) + pow(b,2));

	gsl_matrix_set(R, column1, column1, result);
	gsl_matrix_set(R, column2, column1, 0);

	gsl_matrix_view matrix_substract;

	matrix_substract = gsl_matrix_submatrix(R, column1, column2, 2, R->size2-2);

	gsl_matrix* M1 = gsl_matrix_alloc(2, R->size2-2);

	product_matrix(rotation, &matrix_substract.matrix, M1);

	add_submatrix(R, M1, column1, column2);

}


/*=========================================*/
/*This function makes all computations in order to perform efficient search for all sub-models*/
/*=========================================*/
void efficient_alg()
{
	Model_QR_components matrix_components;

	gsl_matrix* base_R;
	gsl_vector* RSS_models;
	gsl_vector_view view_C;

	QR_decomposition(A, &matrix_components);

	base_R = gsl_matrix_alloc(matrix_components.model->size1, matrix_components.model->size1);


	get_base_R(matrix_components.R, matrix_components.solution, base_R);

	view_C = gsl_matrix_column(base_R, (base_R->size2)-1);

	RSS_models = gsl_vector_alloc(matrix_components.R->size2 - 1);

	efficient_RSS(&view_C.vector,RSS_models, matrix_components.RSS);

	print_vector(RSS_models);

	/*Applying Columns mutations*/

	gsl_matrix* swap_R = gsl_matrix_alloc(base_R->size1, base_R->size2);

	gsl_matrix_memcpy(swap_R, base_R);

	retriangularization_R(swap_R, 1, 2 );

	view_C = gsl_matrix_column(swap_R, (swap_R->size2)-1);

	efficient_RSS(&view_C.vector,RSS_models, matrix_components.RSS);

	printf("\nAfter Switching column 2 with 1\n");
	print_vector(RSS_models);

}





