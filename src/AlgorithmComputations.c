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
/*=========================================*/
/*global data*/
gsl_vector* solution_y;
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

	if(intercept)
	{
		add_intercept(main_model_A);
	}
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

	//QR used to decompose matrix A such that A = Q*R
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

	//vector of approximated solution
	matrix_components->solution = gsl_vector_alloc(
			matrix_components->QR->size2);

	//vector of residual values
	matrix_components->residual = gsl_vector_alloc(
			matrix_components->QR->size1);

	//vector resulting after product:  Q-transpose*Y
	matrix_components->QtransposeY = gsl_vector_alloc(matrix_input->size1);

	gsl_linalg_QR_decomp(matrix_components->QR, matrix_components->tau);

	/*This function unpacks the encoded QR decomposition (QR,tau) into the matrices Q and R, where Q is M-by-M and R is M-by-N.*/
	gsl_linalg_QR_unpack(matrix_components->QR, matrix_components->tau,
			matrix_components->Q, matrix_components->R);

	/*This function calculates approximate solution of equation and residual */
	gsl_linalg_QR_lssolve(matrix_components->QR, matrix_components->tau, solution_y,
			matrix_components->solution, matrix_components->residual);

	matrix_components->RSS = euclidean_norm(matrix_components->residual);

	gsl_vector_memcpy(matrix_components->QtransposeY, solution_y);

	gsl_linalg_QR_QTvec(matrix_components->QR, matrix_components->tau, matrix_components->QtransposeY);

	matrix_components->QtransposeY->size = matrix_input->size2;

}

/*=========================================*/
/* Function to compute compute RSS
 * shorter version, without interest of other components
 * used for compute RSS of all sub-models*/
/*=========================================*/
double RSS_compute(gsl_matrix* QR) {

	double RSS = INIT;

	gsl_vector * x = gsl_vector_alloc(QR->size2);
	gsl_vector * E = gsl_vector_alloc(QR->size1);
	gsl_vector* tau = gsl_vector_alloc(MIN_VALUE(QR->size1, QR->size2));
	gsl_matrix* temp_QR = gsl_matrix_alloc(QR->size1, QR->size2);

	if(QR->size2 > 0)
	{

	gsl_matrix_memcpy(temp_QR, QR);


	gsl_linalg_QR_decomp(temp_QR, tau);

	gsl_linalg_QR_lssolve(temp_QR, tau, solution_y, x, E);

	/* The least squares solution minimizes the Euclidean norm of the residual, ||Ax - b||.*/
	RSS = euclidean_norm(E);

	}
	gsl_vector_free(x);
	gsl_vector_free(E);
	gsl_vector_free(tau);
	gsl_matrix_free(temp_QR);

	return (RSS);
}


/*=========================================*/
/*This function compute all possible subsets of ModelMatrix*/
/*=========================================*/
void naive_alg(void) {

	//remove("Output_Steps");

	gsl_vector* columns_selected;
	T_INDIVIDUAL solution;
	solution.columns = gsl_vector_alloc(main_model_A->size2);

	double result  = 1;
	double new_RSS;


	/*computing all possible combinations of columns*/
	for (uint8 i = 1; i <= main_model_A->size2; i++) {

		gsl_combination * columns_transitions = gsl_combination_alloc(
				main_model_A->size2, i);

		/*set combinations with values from 0 to n-1*/
		gsl_combination_init_first(columns_transitions);

		/*do combination while reaching all possible subsets*/
		do {

			columns_selected = gsl_vector_alloc(columns_transitions->k);
			for(uint8 i = 0; i < columns_selected->size; i++)
			{
				columns_selected->data[i] = columns_transitions->data[i];
			}

			/*get matrix model*/
			gsl_matrix* my_model = sub_model_matrix(columns_selected);

			gsl_combination_fprintf(stdout, columns_transitions, "%u ");

			new_RSS = RSS_compute(my_model);
			printf("\nRSS: %lf\n", new_RSS);

			//if best solution is find then update model
			if (FALSE != criterion(new_RSS, main_model_A->size2, columns_selected->size, &result)) {
				solution.columns->size = columns_selected->size;
				solution.RSS = new_RSS;
				solution.fitness_value = result;
				gsl_vector_memcpy(solution.columns, columns_selected);
			}

		} while (GSL_SUCCESS == gsl_combination_next(columns_transitions));


	}

	printf("\nValue of best solution based on criterion: %lf\nRSS: %lf",
			solution.fitness_value, solution.RSS);
	print_vector(solution.columns);

}

/*=========================================*/
/*This function form a sub-model matrix*/
/*=========================================*/
gsl_matrix* sub_model_matrix(gsl_vector* matrix_combination) {

	/*make a copy of full model matrix*/
	gsl_matrix* matrix_transitions = gsl_matrix_alloc(main_model_A->size1,
			matrix_combination->size);
	gsl_vector * v;

	gsl_matrix_set_zero(matrix_transitions);

	for (uint8 i = 0; i < matrix_transitions->size2; i++) {

		/* take each column index from full model matrix that is present in combination
		 * and add it to sub-model*/

		v = gsl_vector_alloc(main_model_A->size1);

		gsl_matrix_get_col(v, main_model_A,
				gsl_vector_get(matrix_combination, i));

		gsl_matrix_set_col(matrix_transitions, i, v);

	}

	//gsl_vector_free(v);

	return matrix_transitions;

}

/*=========================================*/
/*This function calculates R* using R*x(approximated solution) */
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
/*This function calculates R* using Q-transpose*Y */
/*=========================================*/
void get_base_R_byQ( gsl_matrix* R, gsl_vector* QtY, gsl_matrix* R_result)
{


    gsl_vector_view column;

    for(uint8 i=0; i< (R->size2); i++)
    {
    	column = gsl_matrix_column(R, i);
    	gsl_matrix_set_col(R_result, i, &column.vector);
    }

    QtY->size = R_result->size1;

    gsl_matrix_set_col(R_result, (R_result->size2)-1, QtY);

}

/*=========================================*/
/*This function compute RSS using elements of C*/
/*=========================================*/
void efficient_RSS(gsl_vector* C, gsl_vector* RSS_models, const double RSS)
{

	uint8 n = RSS_models->size;
	double result = RSS;
	double temp_result = 0;

	printf("\nRSS: %lf\n", result);

	for (uint8 i = 0; i < n; i++) {
		double el = gsl_vector_get(C, (n - i));

		printf("\nElement: %lf\n", el);

		result += (double) pow(el, 2);

		temp_result = result;

		gsl_vector_set(RSS_models, i, temp_result);
	}

}

/*=========================================*/
/*This function compute RSS using elements of C in case of columns are deleted*/
/*=========================================*/
void efficient_RSS_D(gsl_vector* C, gsl_vector* RSS_models, const double RSS, uint8 SubModel_size, uint8 Model_size)
{

	double result = RSS;
	double temp_result = 0;

	RSS_models->size = Model_size;

	for (uint8 i = 0; i < Model_size; i++) {
		double el = gsl_vector_get(C, (Model_size - i));

		printf("\nElement: %lf\n", el);

		result += (double) pow(el, 2);

		temp_result = result;

		gsl_vector_set(RSS_models, i, temp_result);

	}

	gsl_vector_reverse(RSS_models);
	RSS_models->size = SubModel_size;
	gsl_vector_reverse(RSS_models);

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
	double a, b, sin, cos, result;
	uint8 row = column + 1;

	a = gsl_matrix_get(R, column, column);
	b = gsl_matrix_get(R, row, column);

	result = (double) sqrt(pow(a, 2) + pow(b, 2));
	result = result * sign_r_element(a, b);

	printf("\na: %lf\n", a);
	printf("b: %lf\n", b);

	if (result == 0) {
		cos = 1;
		sin = 0;
	} else {

		cos = (double) a / result;
		sin = (double) b / result;

	}

	printf("\nsin: %lf\n", sin);
	printf("cos: %lf\n", cos);
	printf("\nres: %lf\n", result);

	gsl_matrix_set_all(rotation, cos);
	gsl_matrix_set(rotation, 0, 1, sin);
	gsl_matrix_set(rotation, 1, 0, (-1) * sin);

	print_matrix(rotation);

	//update matrix R with the elements of vector [r,0]
	gsl_matrix_set(R, column, column, result);
	gsl_matrix_set(R, row, column, 0);
}


/*=========================================*/
/*This function performs re - triangularization after switching 2 columns in R, columns number difference must be 1*/
/*=========================================*/
static void columns_transition_retriangularization_R(gsl_matrix* R, uint8 column1)
{

	uint8 column2 = column1 + 1;

	printf("\nSwap between columns %d and %d .\n", column1, column2);

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

	gsl_matrix_free(rotation);
	gsl_matrix_free(M1);

}

/*=========================================*/
/*This function performs re - triangularization after deleting a column from R*/
/*=========================================*/
static void column_removal_retriangularization_R(gsl_matrix* R, uint8 column1)
{

	printf("\nElimination of column %d.\n", column1);

	uint8 column2 = column1 + 1;

	uint8 max_number_retriangularization;

	gsl_matrix* rotation = gsl_matrix_alloc(NR_ELEMENTS, NR_ELEMENTS);
	gsl_matrix* M1;
	gsl_matrix_view matrix_substract;

	delete_column(R, column1);

	//number of total rotations: number of columns of R* minus 1(columns of C)- index of column removed
	max_number_retriangularization = R->size2 - column1 - 1;

	print_matrix(R);

	for (uint8 i = 0; i < max_number_retriangularization; i++)

	{

		printf("\n################################################");
		printf("\n#################Iteration%d#####################\n",i);
		printf("################################################\n");

		//computes new r element of R*, where rotation is applied in order to get 0 under diagonal
		compute_rotation_angle(R, rotation, column1);

		print_matrix(R);

		//update the rest of 2 rows of R that were affected by re-triangularization
		matrix_substract = gsl_matrix_submatrix(R, column1, column2, 2,
				R->size2 - column2);

		print_matrix(&matrix_substract.matrix);

		M1 = gsl_matrix_alloc(2, matrix_substract.matrix.size2);


		product_matrix(rotation, &matrix_substract.matrix, M1);

		print_matrix(M1);

		add_submatrix(R, M1, column1, column2);

		column1++;
		column2++;
		gsl_matrix_free(M1);

		print_matrix(R);
	}

	gsl_matrix_free(rotation);

}

void get_submodels_Rss(gsl_matrix* Model, void (*f_method)(gsl_matrix* M, uint8 index) , double RSS, gsl_vector* RSS_models)
{
	uint8 Model_columns = Model->size2;

	gsl_vector_view view_C;

	for (uint8 i = 0; i < 2/*Model->size2 - 2*/; i++) {

		(*f_method)(Model, i);

		view_C = gsl_matrix_column(Model, (Model->size2) - 1);

		if (Model_columns != Model->size2) {

			efficient_RSS_D(&view_C.vector, RSS_models, RSS, Model->size2 - 1,
					Model_columns - 1);

		} else {

			efficient_RSS(&view_C.vector, RSS_models, RSS);
		}

		printf("\nRSS of new model.\n");

		print_vector(RSS_models);

	}
}

/*=========================================*/
/*This function makes all computations in order to perform efficient search for all sub-models*/
/*=========================================*/
void efficient_alg(T_EFFICIENT_METHOD method)
{
	Model_QR_components matrix_components;

	gsl_matrix* base_R;
	gsl_matrix* base_R_byQ;
	gsl_vector* RSS_models;

	print_matrix(main_model_A);

	QR_decomposition(main_model_A, &matrix_components);

	const uint8 nRows = matrix_components.R->size1;
	const uint8 nColumns = matrix_components.R->size2;

	base_R = gsl_matrix_alloc(nRows, nColumns+1);
	base_R_byQ = gsl_matrix_alloc(nRows, nColumns+1);

	get_base_R(matrix_components.R, matrix_components.solution, base_R);
	get_base_R_byQ(matrix_components.R, matrix_components.QtransposeY, base_R_byQ);

	printf("\nColumn of data for RSS computation, using C = R*x(approximated solution)\n");
	print_matrix(base_R);

	gsl_matrix* elimination_R = gsl_matrix_alloc(base_R->size1, base_R->size2);
	gsl_matrix_memcpy(elimination_R, base_R);

	gsl_vector_view view_C = gsl_matrix_column(elimination_R, (elimination_R->size2) - 1);

	uint8 nSubmodels = elimination_R->size2 - 2;
			RSS_models = gsl_vector_alloc(nSubmodels);

	efficient_RSS(&view_C.vector, RSS_models, matrix_components.RSS);

	printf("\nRSS of Main model.\n");

	print_vector(RSS_models);

	//
	method = 1;

	switch (method) {
	case columns_transitions: {
		gsl_matrix* swap_R = gsl_matrix_alloc(base_R_byQ->size1, base_R_byQ->size2);
		gsl_matrix_memcpy(swap_R, base_R_byQ);

		//print_matrix(swap_R);

		uint8 nSubmodels = swap_R->size2 - 2;
		RSS_models = gsl_vector_alloc(nSubmodels);

		printf("\n################################################");
		printf("\n##########Applying Columns Transitions##########\n");
		printf("################################################\n");

		//Applying re-triangularization after switching 2 columns
		get_submodels_Rss(swap_R,columns_transition_retriangularization_R, matrix_components.RSS, RSS_models);

		break;
	}
	case columns_removal: {
		gsl_matrix* elimination_R = gsl_matrix_alloc(base_R->size1, base_R->size2);
		gsl_matrix_memcpy(elimination_R, base_R);

		//print_matrix(elimination_R);

		uint8 nSubmodels = elimination_R->size2 - 1;
		RSS_models = gsl_vector_alloc(nSubmodels);

		printf("\n################################################");
		printf("\n##########Applying Columns Removal##############\n");
		printf("################################################\n");

		//Applying re-triangularization after deleting a column
		get_submodels_Rss(elimination_R,column_removal_retriangularization_R, matrix_components.RSS, RSS_models);

		break;
	}
	default: {
		//Nothing To Do
		break;
	}
	}

	gsl_matrix_free(base_R);
}

