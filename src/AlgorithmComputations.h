#ifndef SRC_ALGORITHMCOMPUTATIONS_H_
#define SRC_ALGORITHMCOMPUTATIONS_H_

/*=========================================*/
/*include*/
#include "FileOperations.h"
#include "MatrixComputations.h"
#include "Types.h"
/*=========================================*/
/*define*/
#define MIN_VALUE(X, Y)  ((X) < (Y) ? (X) : (Y))

/*=========================================*/
/*enumerations*/
/*=========================================*/
/*typedef*/

/*=========================================*/
/*external functions*/
void QR_decomposition(void);
void compute_transitions(void);
gsl_matrix* sub_model_matrix(gsl_combination* matrix_combination);
void set_y_vector(T_FILE_DIM* file_dim);
void set_A_matrix(T_FILE_DIM* file_dim);
/*=========================================*/

#endif /* SRC_ALGORITHMCOMPUTATIONS_H_ */
