#ifndef MATRIXCOMPUTATIONS_H
#define MATRIXCOMPUTATIONS_H
/*=========================================*/
/*include*/
#include<gsl\gsl_linalg.h>
#include<gsl\gsl_matrix.h>
#include<gsl\gsl_vector.h>
#include<gsl/gsl_combination.h>
#include "Types.h"
#include "FileOperations.h"
/*=========================================*/
/*define*/
#define IDENTITY_VALUE 1u
#define INVERSE_MULTIPLI(value) (1/(value))
#define SCALE_2 2u
/*=========================================*/
/*enumerations*/
/*=========================================*/
/*typedef*/
/*=========================================*/
/*external functions*/
void 	print_matrix(gsl_matrix *M);
void 	convert_to_gsl(T_FILE_DIM* file_dim, gsl_matrix *M);
void 	print_vector (gsl_vector *V);
boolean compute_matrix_inverse(gsl_matrix * M);
boolean product_matrix_vector(gsl_matrix* M, gsl_vector* v, gsl_vector* result);
boolean product_matrix(gsl_matrix* M, gsl_matrix* M2, gsl_matrix* result);
void    add_submatrix(gsl_matrix* R, gsl_matrix* sub_matrix, uint8 index1, uint8 index2);
void    delete_column(gsl_matrix* R, uint8 col);
double  euclidean_norm(gsl_vector* V);
gsl_matrix* 	add_intercept(gsl_matrix* R);

/*=========================================*/
#endif
