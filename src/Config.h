/*
 * Config.h
 *
 *  Created on: Apr 6, 2019
 *      Author: Paus
 */

#ifndef SRC_CONFIG_H_
#define SRC_CONFIG_H_
/*=========================================*/
/*define*/
/*===========Genetic Algorithm=============*/
#define PERCENTAJE_OF_CHROMOSOMES 70u //70
#define PERCENTAJE_OF_GENES 55u//55u
#define CONVERGE 100u
#define MAX_FITNESS 100000
#define FIXED_NR_GENES 0u
/*=========================================*/
/*==========Tournament Selection===========*/
/*=========================================*/
#define PERCENTAJE_OF_TOURNAMENT_K 30u
/*=========================================*/
/*==========Simulated Annealing============*/
/*=========================================*/
#define TEMP 100000 //10000 100000 100000
#define COOLING_RATE 0.95//0.003 0.98 0.330 0.729 0.85/ 0.86
#define NR_ITERATIONS 500u
/*=========================================*/
/*================Crossover================*/
/*=========================================*/
#define REFERENCE_PROBABILITY 0.5
#define NR_OF_ATTEMPTS 10u
/*=========================================*/
/*macro definition function*/
/*=========================================*/
#define PERCENTAJE(x, y) ( ((x) > 0 && (y) > 0) ? ((float)(x)/100) * (y) : 0)
/*=========================================*/

#endif /* SRC_CONFIG_H_ */
