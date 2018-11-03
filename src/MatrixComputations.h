#ifndef MATRIXCOMPUTATIONS_H
#define MATRIXCOMPUTATIONS_H
/*=========================================*/
/*include*/
#include<gsl\gsl_linalg.h>
#include<gsl\gsl_matrix.h>
#include<gsl\gsl_vector.h>
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
void print_matrix(gsl_matrix *M);
void convert_to_gsl(T_FILE_DIM* file_dim, gsl_matrix *M);
boolean compute_matrix_inverse(gsl_matrix * M);
double euclidean_norm(gsl_vector* V, uint8 i);
/*=========================================*/
#endif
