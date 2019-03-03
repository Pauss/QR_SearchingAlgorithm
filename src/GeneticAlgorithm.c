/*
 * GeneticAlgorithm.c
 *
 *  Created on: Mar 3, 2019
 *      Author: Paus
 */
/*=========================================*/
/*include*/
#include "GeneticAlgorithm.h"
/*=========================================*/
/*private data*/
/*=========================================*/
/*global data*/
/*=========================================*/
/*private functions*/
/*=========================================*/


/*=========================================*/
/*This function compute a possible subset of ModelMatrix*/
/*=========================================*/
boolean get_random_model(gsl_vector* my_random_model, uint8 n, uint8 size_A)
{
	uint8 i = n;
	uint8 found = 0;

	gsl_vector_set_zero(my_random_model);

	if (n >= my_random_model->size) {
		return FALSE;
	} else {
		while (i) {

			uint8 r = rand() % (size_A);

			//for getting n random columns
			while (gsl_vector_get(my_random_model, r) == 1) {
				r = rand() % (size_A);
			}

			gsl_vector_set(my_random_model, r, 1);
			i--;
		}
	}

	for (uint8 j = 0; j < my_random_model->size; j++) {
		//if column was selected for sub-model set it with column number
		if (gsl_vector_get(my_random_model, j) == 1) {

			gsl_vector_set(my_random_model, found, j);
			found++;
		}
	}

	my_random_model->size = found;

	return TRUE;

}

/*=========================================*/
/*This function applies an criterion over sub-models*/
/*=========================================*/
boolean criterion(double RSS, uint8 k, double* best_solution) {
	double result;

	result = (double) ((-2) * log(RSS)) + 2 * k;

	if (result < *best_solution) {
		*best_solution = result;

		return TRUE;

	}

	return FALSE;

}


/*=========================================*/
/*This function applies a GA for getting best solution*/
/*=========================================*/
void GA_alg(void)
{
	gsl_matrix * main_model = get_A_matrix();
	gsl_vector* random_model = gsl_vector_alloc(main_model->size2);
	gsl_matrix* sub_model;
	uint8 model_size_k = NUMBER_OF_GENES;
	uint8 iterations = NUMBER_OF_CHROMOSOMES;
	double RSS_model;

	Model_QR_components matrix_components;
	QR_decomposition(main_model, &matrix_components);

	T_BEST_SUBMODEL best_solution;
	best_solution.columns = gsl_vector_alloc(main_model->size2);
	double new_value  = 100;

	do {

		if (FALSE
				!= get_random_model(random_model, model_size_k,
						main_model->size2)) {
			sub_model = sub_model_matrix(random_model);

			print_vector(random_model);

			RSS_model = RSS_compute(sub_model);
			printf("RSS: %lf\n", RSS_model);

			//if best solution is find then update model
			if (criterion(RSS_model, model_size_k, &new_value)) {

				best_solution.value = new_value;
				best_solution.columns->size = random_model->size;
				best_solution.RSS = RSS_model;
				gsl_vector_memcpy(best_solution.columns, random_model);
			}

			//if size changes
			random_model->size = main_model->size2;
		}

		iterations--;

	} while (iterations > 0);

	printf("\nValue of best solution based on criterion: %lf\nRSS: %lf",
			best_solution.value, best_solution.RSS);
	print_vector(best_solution.columns);



}

