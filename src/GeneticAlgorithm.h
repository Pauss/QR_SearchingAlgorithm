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
#define TOURNAMENT_K 6
#define CONVERGE 50

/*=========================================*/
/*enumerations*/
/*=========================================*/

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
void 	naive_GA_alg(void);
void print_individual(T_INDIVIDUAL individual);
/*=========================================*/

#endif /* SRC_GENETICALGORITHM_H_ */
