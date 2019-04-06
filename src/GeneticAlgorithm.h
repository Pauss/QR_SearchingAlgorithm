/*
 * GeneticAlgorithm.h
 *
 *  Created on: Mar 3, 2019
 *      Author: Paus
 */

#ifndef SRC_GENETICALGORITHM_H_
#define SRC_GENETICALGORITHM_H_

/*=========================================*/
/*include*/
#include "QR_SearchingAlgorithm.h"
#include "MatrixComputations.h"
#include "AlgorithmComputations.h"
#include "Types.h"
#include <gsl/gsl_combination.h>
#include <Math.h>
#include <stdlib.h>
#include <time.h>
/*=========================================*/
/*define*/
/*=========================================*/
/*enumerations*/

typedef enum
{
	roulette_wheel = 0,
	tournament
}T_SELECTION_METHOD;

typedef enum
{
	//mutations
	flip = 0u,
	interchanging = 1u,
	interchanging_abs = 2u,
	reversing = 3u,
	//crossovers
	_1point_simple = 4u,
	_1point = 5u,
	RRC = 6u,
	uniform = 7u,
	//NONE
	no_operator = 8u
}T_OPERATOR_METHOD;

typedef enum
{
	data_no_error = 0,
	data_index_out_of_range,
	data_nr_of_columns_0,
	data_nr_of_columns_high,
	data_nr_of_individuals_0,
	data_nr_of_individuals_high

}DATA_ERRORS;

/*=========================================*/
/*typedef*/

typedef struct{
	boolean selected;
	double fitness_value;
	double RSS;
	double selection_probability;
	gsl_vector* columns;
	gsl_vector* bit_columns;
	gsl_matrix* submodel;
}T_INDIVIDUAL;
/*=========================================*/
/*external functions*/
boolean get_random_model(gsl_vector* my_random_model, uint16 n, uint16 size_A);
boolean criterion (double RSS, uint16 n, uint16 k, double* result);
void fitness_func(T_INDIVIDUAL* individual, uint16 model_size_n, uint16 model_size_k, double* result);
void GA_naive_alg(T_SELECTION_METHOD method, T_OPERATOR_METHOD op1,  T_OPERATOR_METHOD op2);
void GA_simulated_annealing(T_OPERATOR_METHOD op1, T_OPERATOR_METHOD op2);
/*=========================================*/

#endif /* SRC_GENETICALGORITHM_H_ */
