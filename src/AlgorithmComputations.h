#ifndef SRC_ALGORITHMCOMPUTATIONS_H_
#define SRC_ALGORITHMCOMPUTATIONS_H_

/*=========================================*/
/*include*/
#include "FileOperations.h"
#include "MatrixComputations.h"
#include "Types.h"
#include <gsl/gsl_combination.h>
#include "GeneticAlgorithm.h"
/*=========================================*/
/*define*/
#define MIN_VALUE(X, Y)  ((X) < (Y) ? (X) : (Y))
#define MAX_VALUE(X, Y)  ((X) < (Y) ? (Y) : (X))
#define ABS_VALUE(X) ((X) < 0 ? -(X) : (X))
#define NR_ELEMENTS 2u

/*=========================================*/
/*enumerations*/

typedef enum
{
	columns_transitions = 0,
	columns_removal
}T_EFFICIENT_METHOD;
/*=========================================*/
/*typedef*/

typedef struct{
	gsl_matrix* model;
	gsl_matrix* QR;
	gsl_matrix* Q;
	gsl_matrix* R;
	gsl_vector* tau;
	gsl_vector* residual;
	gsl_vector* solution;
	gsl_vector* QtransposeY;
	double RSS;
}Model_QR_components;
/*=========================================*/
/*external functions*/
void 				 QR_decomposition(gsl_matrix* matrix_input, Model_QR_components* matrix_components);
void 				 naive_alg(void);
void 				 set_y_vector(T_FILE_DIM* file_dim);
void 			     set_A_matrix(T_FILE_DIM* file_dim);
void 			     efficient_alg(T_EFFICIENT_METHOD method);
gsl_matrix* 		 get_A_matrix();
gsl_matrix* 		 sub_model_matrix(gsl_vector* matrix_combination);
double 				 RSS_compute(gsl_matrix* QR);
extern gsl_vector* solution_y;
/*=========================================*/

#endif /* SRC_ALGORITHMCOMPUTATIONS_H_ */
