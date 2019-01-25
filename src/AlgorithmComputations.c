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
static gsl_matrix* main_model_A;
static gsl_vector* solution_y;
/*=========================================*/
/*global data*/
/*=========================================*/
/*private functions*/
/*=========================================*/
static void column_removal_retriangularization_R(gsl_matrix* R, uint8 column1);
static void columns_transition_retriangularization_R(gsl_matrix* R, uint8 column1);

/*=========================================*/
/*Get 'A' in order to solve y=A*x */
/*=========================================*/
gsl_matrix* get_A_matrix() {
	return main_model_A;
}

/*=========================================*/
/*Get 'y' vector in order to solve y=A*x */
/*=========================================*/
gsl_vector* get_y_vector() {
	return solution_y;
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
	main_model_A = gsl_matrix_alloc(file_dim->lines, file_dim->columns - 1);
	gsl_matrix_memcpy(main_model_A, &A_view.matrix);
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
	solution_y = gsl_vector_alloc(file_dim->lines);
	gsl_vector_memcpy(solution_y, &y_view.vector);
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
	gsl_linalg_QR_lssolve(matrix_components->QR, matrix_components->tau, solution_y,
			matrix_components->solution, matrix_components->residual);

	matrix_components->RSS = euclidean_norm(matrix_components->residual);

}

/*=========================================*/
/* Function to compute compute RSS
 * shorter version, without interest of other components
 * used for compute RSS of all sub-models*/
/*=========================================*/
double RSS_compute(gsl_matrix* QR) {

	gsl_vector * x = gsl_vector_alloc(QR->size2);
	gsl_vector * E = gsl_vector_alloc(QR->size1);
	gsl_vector* tau = gsl_vector_alloc(MIN_VALUE(QR->size1, QR->size2));
	double RSS = INIT;

	gsl_linalg_QR_decomp(QR, tau);

	gsl_linalg_QR_lssolve(QR, tau, solution_y, x, E);

	/* The least squares solution minimizes the Euclidean norm of the residual, ||Ax - b||.*/
	RSS = euclidean_norm(E);

	return (RSS);
}

/*=========================================*/
/*This function compute all possible subsets of ModelMatrix*/
/*=========================================*/
void compute_transitions_QR(void) {
	/*ToDo*/
	/*Error handling*/

	//remove("Output_Steps");

	/*computing all possible combinations of columns*/
	for (uint8 i = 1; i <= main_model_A->size2; i++) {

		gsl_combination * columns_transitions = gsl_combination_alloc(
				main_model_A->size2, i);

		/*set combinations with values from 0 to n-1*/
		gsl_combination_init_first(columns_transitions);

		/*do combination while reaching all possible subsets*/
		do {

			/*get matrix model*/
			gsl_matrix* my_model = sub_model_matrix(columns_transitions);
			/*print_steps(RSS_compute(my_model), columns_transitions);*/

			gsl_combination_fprintf(stdout, columns_transitions, "%u ");
			printf("\nRSS: %lf\n", RSS_compute(my_model));

		} while (GSL_SUCCESS == gsl_combination_next(columns_transitions));
	}
}

/*=========================================*/
/*This function form a sub-model matrix*/
/*=========================================*/
gsl_matrix* sub_model_matrix(gsl_combination* matrix_combination) {

	/*make a copy of full model matrix*/
	gsl_matrix* matrix_transitions = gsl_matrix_alloc(main_model_A->size1,
			matrix_combination->k);

	gsl_matrix_set_zero(matrix_transitions);

	for (uint8 i = 0; i < matrix_transitions->size2; i++) {

		/* take each column index from full model matrix that is present in combination
		 * and add it to sub-model*/
		gsl_combination_get(matrix_combination, i);
		gsl_vector * v = gsl_vector_alloc(main_model_A->size1);
		gsl_matrix_get_col(v, main_model_A,
				gsl_combination_get(matrix_combination, i));
		gsl_matrix_set_col(matrix_transitions, i, v);

	}

	return matrix_transitions;

}

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
	double result = pow(RSS,2);
	double temp_result = 0;

	for(uint8 i=0; i < n; i++)
	{
		result += pow(gsl_vector_get(C,(n-i)),2);

		temp_result = sqrt(result);

		gsl_vector_set(RSS_models,i, temp_result);
	}

}

/*=========================================*/
/*This function calculates Sign of new element when Re-triangularization is applied*/
/*=========================================*/
int8 sign_r_element(double a, double b)
{
	int8 sign = 1;

	if (abs(a) > abs(b)) {
		if (a < 0) {
			sign = -1;
		}
	} else {
		if (b < 0) {
			sign = -1;
		}
	}

	return (sign);
}

/*=========================================*/
/*This function calculates the two component cos and sign of rotation*/
/*=========================================*/
void compute_rotation_angle(gsl_matrix* R, gsl_matrix* rotation, uint8 column)
{
	double l_a, l_b, l_sin, l_cos, l_result;
	uint8 row = column + 1;

	l_a = gsl_matrix_get(R, column, column);
	l_b = gsl_matrix_get(R, row, column);

	l_result = (double) sqrt(pow(l_a, 2) + pow(l_b, 2));
	l_result = l_result * sign_r_element(l_a, l_b);

	if (l_result == 0) {
		l_cos = 1;
		l_sin = 0;
	} else {

		gsl_linalg_givens(l_a, l_b, &l_cos, &l_sin);

		l_cos = ABS_VALUE(l_cos);
		l_sin = ABS_VALUE(l_sin);
	}

	gsl_matrix_set_all(rotation, l_cos);
	gsl_matrix_set(rotation, 0, 1, l_sin);
	gsl_matrix_set(rotation, 1, 0, -l_sin);

	//update matrix R with the elements of vector [r,0]
	gsl_matrix_set(R, column, column, l_result);
	gsl_matrix_set(R, row, column, 0);
}


/*=========================================*/
/*This function performs re - triangularization after switching 2 columns in R, columns number difference must be 1*/
/*=========================================*/
static void columns_transition_retriangularization_R(gsl_matrix* R, uint8 column1)
{
	uint8 column2 = column1 + 1;
	gsl_matrix* rotation = gsl_matrix_alloc(NR_ELEMENTS, NR_ELEMENTS);;
	gsl_matrix* M1;
	gsl_matrix_view matrix_substract;

	gsl_matrix_swap_columns(R, column1, column2);

	//computes new r element of R, where rotation is applied in order to get 0 under diagonal
	compute_rotation_angle(R, rotation, column1);

	/*update the rest of 2 rows of R that were affected by re-triangularization*/
	matrix_substract = gsl_matrix_submatrix(R, column1, column2, 2, R->size2-column2);

	M1 = gsl_matrix_alloc(2, matrix_substract.matrix.size2);

	product_matrix(rotation, &matrix_substract.matrix, M1);

	add_submatrix(R, M1, column1, column2);

}

/*=========================================*/
/*This function performs re - triangularization after deleting a column from R*/
/*=========================================*/
static void column_removal_retriangularization_R(gsl_matrix* R, uint8 column1)
{

	uint8 column2 = column1 + 1;
	uint8 max_number_retriangularization = R->size2 - column1;

	gsl_matrix* rotation = gsl_matrix_alloc(NR_ELEMENTS, NR_ELEMENTS);
	gsl_matrix* M1;
	gsl_matrix_view matrix_substract;

	delete_column(R, column1);

	for (uint8 i = 0; i < max_number_retriangularization - 2; i++)

	{
		//computes new r element of R, where rotation is applied in order to get 0 under diagonal
		compute_rotation_angle(R, rotation, column1);

		//update the rest of 2 rows of R that were affected by re-triangularization
		matrix_substract = gsl_matrix_submatrix(R, column1, column2, 2,
				R->size2 - column2);

		M1 = gsl_matrix_alloc(2, matrix_substract.matrix.size2);

		product_matrix(rotation, &matrix_substract.matrix, M1);

		add_submatrix(R, M1, column1, column2);

		column1++;
		column2++;
	}

}

/*=========================================*/
/*This function makes all computations in order to perform efficient search for all sub-models*/
/*=========================================*/
void efficient_alg(void)
{
	Model_QR_components matrix_components;

	gsl_matrix* base_R;
	gsl_matrix* swap_R;
	gsl_vector* RSS_models;
	gsl_vector_view view_C;

	QR_decomposition(main_model_A, &matrix_components);

	base_R = gsl_matrix_alloc(matrix_components.model->size1,
			matrix_components.model->size1);
	RSS_models = gsl_vector_alloc(matrix_components.R->size2 - 1);

	get_base_R(matrix_components.R, matrix_components.solution, base_R);

	swap_R = gsl_matrix_alloc(base_R->size1, base_R->size2);
	gsl_matrix_memcpy(swap_R, base_R);

	/*Applying re-triangularization after switching 2 columns*/
	for (uint8 i = 0; i < 1/*swap_R->size2 - 2*/; i++) {

		columns_transition_retriangularization_R(swap_R, i);

		view_C = gsl_matrix_column(swap_R, (swap_R->size2) - 1);

		efficient_RSS(&view_C.vector, RSS_models, matrix_components.RSS);

		print_matrix(swap_R);

		printf("\nAfter Switching column %d with %d\n", i, i + 1);

		print_vector(RSS_models);
	}

	swap_R = gsl_matrix_alloc(base_R->size1, base_R->size2);
	gsl_matrix_memcpy(swap_R, base_R);

	//Applying re-triangularization after deleting a column
	/*	for (uint8 i = 0; i < base_R->size1 - 2; i++) {

	 column_removal_retriangularization_R(swap_R, 1);

	 view_C = gsl_matrix_column(swap_R, (swap_R->size2) - 1);

	 efficient_RSS(&view_C.vector, RSS_models, matrix_components.RSS);

	 printf("\nAfter Deleting column %d\n", i);

	 print_matrix(swap_R);

	 print_vector(RSS_models);

	 }*/

}

