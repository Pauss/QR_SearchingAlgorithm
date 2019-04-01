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
/*===========Genetic Algorithm=============*/
#define PERCENTAJE_OF_CHROMOSOMES 70u
#define PERCENTAJE_OF_GENES 55u//55u
#define CONVERGE 30u
#define MAX_FITNESS 100000
/*==========Tournament Selection===========*/
#define PERCENTAJE_OF_TOURNAMENT_K 30u
/*==========Simulated Annealing============*/
#define TEMP 100000 //10000 100000 100000
#define COOLING_RATE 0.86//0.003 0.98 0.330 0.729 0.85
#define NR_ITERATIONS 200u
/*================Crossover================*/
#define REFERENCE_PROBABILITY 0.5
#define NR_OF_ATTEMPTS 10u
/*macro definition function*/
#define PERCENTAJE(x, y) ( ((x) > 0 && (y) > 0) ? ((float)(x)/100) * (y) : 0)
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
boolean get_random_model(gsl_vector* my_random_model, uint8 n, uint8 size_A);
boolean criterion (double RSS, uint8 n, uint8 k, double* result);
void fitness_func(T_INDIVIDUAL* individual, uint8 model_size_n, uint8 model_size_k, double* result);
void GA_naive_alg(T_SELECTION_METHOD method, T_OPERATOR_METHOD op1,  T_OPERATOR_METHOD op2);
void GA_simulated_annealing(T_OPERATOR_METHOD op1, T_OPERATOR_METHOD op2);
/*=========================================*/

#endif /* SRC_GENETICALGORITHM_H_ */
