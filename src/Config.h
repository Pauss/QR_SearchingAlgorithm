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
#define PERCENTAJE_OF_CHROMOSOMES 100u //70
#define CONVERGE 50u
#define MAX_FITNESS 100000
#define FIXED_NR_GENES 0u
#define PERCENTAJE_OF_GENES 30u//55u
/*=========================================*/
/*===========BUilding BLocks===============*/
/*=========================================*/
#define PERCENTAJE_OF_SCHEMAS 70u //70
/*=========================================*/
/*==========Tournament Selection===========*/
/*=========================================*/
#define PERCENTAJE_OF_TOURNAMENT_K 40u
/*=========================================*/
/*==========Simulated Annealing============*/
/*=========================================*/
#define TEMP 1000000 //10000 100000 100000 // 10000000
#define COOLING_RATE 0.75//0.003 0.98 0.330 0.729 0.85/ 0.86 /0.95
/*This parameter is also used in Hill Climbing*/
#define NR_ITERATIONS 1000u
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
