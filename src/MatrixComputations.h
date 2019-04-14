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
#include <stdio.h>
#include <math.h>
/*=========================================*/
/*define*/
#define SCALE_2 2u
#define INTERCEPT_VALUE 1u
/*=========================================*/
/*enumerations*/
/*=========================================*/
/*typedef*/
/*=========================================*/
/*external functions*/
void 	print_matrix(gsl_matrix *M);
void 	convert_to_gsl(T_FILE_DIM* file_dim, gsl_matrix *M);
void    copy_vector8(uint8 *dest, uint8* src, uint16 size);
void    copy_vector16(uint16 *dest, uint16* src, uint16 size);
void 	print_vector16 (uint16 *v, uint16 size);
void 	print_vector8 (uint8 *v, uint16 size);
void 	print_vector (gsl_vector *v);
void    add_submatrix(gsl_matrix* R, gsl_matrix* sub_matrix, uint8 index1, uint8 index2);
void    delete_column(gsl_matrix* R, uint16 col);
void 	vector_swap(uint8 *vector, uint8 *vector2, uint16 size);
void 	vector_swap_elements(uint8 *vector, uint16 index1, uint16 index2);
double  euclidean_norm(gsl_vector* V);
boolean compute_matrix_inverse(gsl_matrix * M);
boolean product_matrix_vector(gsl_matrix* M, gsl_vector* v, gsl_vector* result);
boolean product_matrix(gsl_matrix* M, gsl_matrix* M2, gsl_matrix* result);
gsl_matrix* add_intercept(gsl_matrix* R);

/*=========================================*/
#endif
