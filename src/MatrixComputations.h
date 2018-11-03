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
/*=========================================*/
/*enumerations*/
/*=========================================*/
/*typedef*/
/*=========================================*/
/*external functions*/
void print_matrix(gsl_matrix *M);
void convert_to_gsl(T_FILE_DIM* file_dim, gsl_matrix *M);
/*=========================================*/
#endif
