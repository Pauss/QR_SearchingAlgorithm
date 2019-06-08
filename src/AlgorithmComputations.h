#ifndef SRC_ALGORITHMCOMPUTATIONS_H_
#define SRC_ALGORITHMCOMPUTATIONS_H_

/*=========================================*/
/*include*/
#include "FileOperations.h"
#include "MatrixComputations.h"
#include "Types.h"
#include <gsl/gsl_combination.h>
#include "GeneticAlgorithm.h"
#include "QR_SearchingAlgorithm.h"
/*=========================================*/
/*define*/
#define MIN_VALUE(X, Y)  ((X) < (Y) ? (X) : (Y))
#define MAX_VALUE(X, Y)  ((X) < (Y) ? (Y) : (X))
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
/*external variables*/
extern gsl_vector* solution_y;
/*=========================================*/
/*external functions*/
void 				 QR_decomposition(gsl_matrix* matrix_input, Model_QR_components* matrix_components);
void 				 naive_alg(void);
void 				 set_y_vector(T_FILE_DIM* file_dim);
void 			     set_A_matrix(T_FILE_DIM* file_dim);
void 			     efficient_alg(T_EFFICIENT_METHOD method);
gsl_matrix* 		 get_A_matrix();
gsl_matrix* 		 submodel_matrix(uint16* matrix_combination, uint16 size);
void 				 submodel_matrix2(size_t* matrix_combination, gsl_matrix* matrix_transitions);
double 				 RSS_compute(gsl_matrix* QR);
Model_QR_components* get_model_elements(void);
void 				 set_model_elements(void);
void 				 column_removal_retriangularization_R(gsl_matrix* R, uint16 column1);
void				 get_submodels_Rss(gsl_matrix* Model, void (*f_method)(gsl_matrix* M, uint16 index) , double RSS, gsl_vector* RSS_models, uint16 column, uint16 Model_columns);
/*=========================================*/

#endif /* SRC_ALGORITHMCOMPUTATIONS_H_ */
