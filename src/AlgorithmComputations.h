#ifndef SRC_ALGORITHMCOMPUTATIONS_H_
#define SRC_ALGORITHMCOMPUTATIONS_H_

/*=========================================*/
/*include*/
#include "FileOperations.h"
#include "MatrixComputations.h"
#include "Types.h"
#include <gsl/gsl_combination.h>
/*=========================================*/
/*define*/
#define MIN_VALUE(X, Y)  ((X) < (Y) ? (X) : (Y))

/*=========================================*/
/*enumerations*/
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
	double RSS;
}Model_QR_components;

/*=========================================*/
/*external functions*/
void 				 QR_decomposition(gsl_matrix* matrix_input, Model_QR_components* matrix_components);
void 				 compute_transitions_QR(void);
void 				 set_y_vector(T_FILE_DIM* file_dim);
void 			     set_A_matrix(T_FILE_DIM* file_dim);
gsl_matrix* 		 get_A_matrix();
gsl_matrix* 		 sub_model_matrix(gsl_combination* matrix_combination);

/*=========================================*/

#endif /* SRC_ALGORITHMCOMPUTATIONS_H_ */
