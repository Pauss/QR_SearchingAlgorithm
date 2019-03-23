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
#define TEMP 10u
#define MIN_FITNESS 1000
#define TOURNAMENT_K 4
#define CONVERGE 20

/*=========================================*/
/*enumerations*/

typedef enum
{
	roulette_wheel = 0,
	tournament
}T_SELECTION_METHOD;

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
void naive_GA_alg(T_SELECTION_METHOD method);
void print_individual(T_INDIVIDUAL individual);
void copy_individual(T_INDIVIDUAL* dest, T_INDIVIDUAL* src);
void copy_individual_into_population(T_INDIVIDUAL* dest, T_INDIVIDUAL* src, uint8 index1,uint8 index2);
void copy_population(T_INDIVIDUAL* dest, T_INDIVIDUAL* src, uint8 new_size);
void new_population_computed(T_INDIVIDUAL* temp_population, uint8 new_size, uint8 n);
void probability_selection(T_INDIVIDUAL* population, uint8 size);
void population_selected(T_INDIVIDUAL* population,
		T_INDIVIDUAL* temp_population, uint8 size, uint8* new_size);
/*=========================================*/

#endif /* SRC_GENETICALGORITHM_H_ */
