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
#include "MatrixComputations.h"
#include "AlgorithmComputations.h"
#include "Types.h"
#include <gsl/gsl_combination.h>
#include <Math.h>
/*=========================================*/
/*define*/
#define NUMBER_OF_CHROMOSOMES 20u
#define NUMBER_OF_GENES 5u

/*=========================================*/
/*enumerations*/
/*=========================================*/
/*typedef*/
typedef struct{
	double value;
	double RSS;
	gsl_vector* columns;
}T_BEST_SUBMODEL;
/*=========================================*/
/*external functions*/
boolean get_random_model(gsl_vector* my_random_model, uint8 n, uint8 size_A);
void 	GA_alg(void);
boolean criterion (double RSS, uint8 k, double* best_solution);
/*=========================================*/

#endif /* SRC_GENETICALGORITHM_H_ */
