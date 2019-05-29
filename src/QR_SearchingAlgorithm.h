#ifndef QR_SEARCHINGALGORITHM_H_
#define QR_SEARCHINGALGORITHM_H_

/*=========================================*/
/*include*/
#include <string.h>
#include <ctype.h>
#include "Types.h"
#include "Config.h"
#include "GeneticAlgorithm.h"
#include "MatrixComputations.h"
#include "AlgorithmComputations.h"
#include "FileOperations.h"
#include <time.h>
#include <sys/time.h>
/*=========================================*/
/*define*/
#define SEED_RAND 1000u
#define NR_PARAMETERS 6u
/*=========================================*/
/*enumerations*/
typedef enum
{
	not_installed = 0,
	installed
}T_Intercept;
/*=========================================*/
/*typedef*/
/*=========================================*/
/*external functions*/
extern T_Intercept intercept;
/*=========================================*/
#endif /* SRC_QR_SEARCHINGALGORITHM_H_ */
