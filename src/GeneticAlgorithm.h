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
#define NUMBER_OF_CHROMOSOMES 10u
#define NUMBER_OF_GENES 7u
#define TEMP 30u
#define MIN_FITNESS 1000u
#define TOURNAMENT_K 4u
#define CONVERGE 200u
#define REFERENCE_PROBABILITY 0.5


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
	interchanghing = 1u,
	reversing = 2u,
	//crossovers
	_1point = 3u,
	RRC = 4u,
	uniform = 5u,
	//NONE
	no_operator = 6u
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
void 	GA_alg(void);
boolean criterion (double RSS, uint8 n, uint8 k, double* result);
void naive_GA_alg(T_SELECTION_METHOD method, T_OPERATOR_METHOD op1,  T_OPERATOR_METHOD op2);
void print_individual(T_INDIVIDUAL individual);
void copy_individual(T_INDIVIDUAL* dest, T_INDIVIDUAL* src);
void copy_individual_into_population(T_INDIVIDUAL* dest, T_INDIVIDUAL* src, uint8 index1,uint8 index2);
void copy_population(T_INDIVIDUAL* dest, T_INDIVIDUAL* src, uint8 new_size);
void new_population_computed(T_INDIVIDUAL* temp_population, uint8 new_size, uint8 n, T_OPERATOR_METHOD op);
void probability_selection(T_INDIVIDUAL* population, uint8 size);
void population_selected(T_INDIVIDUAL* population,
		T_INDIVIDUAL* temp_population, uint8 size, uint8* new_size);
void print_population(T_INDIVIDUAL* population, uint8 size);
void crossover_1point(gsl_vector* bit_individ_1, gsl_vector *bit_individ_2);
void crossover_uniform(gsl_vector* bit_individ_1, gsl_vector *bit_individ_2);
void crossover_RRC(gsl_vector* bit_individ_1, gsl_vector *bit_individ_2);
/*=========================================*/

#endif /* SRC_GENETICALGORITHM_H_ */
