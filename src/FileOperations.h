#pragma once
#ifndef FILEOPERATIONS_H
#define FILEOPERATIONS_H
/*=========================================*/
/*include*/
#include"Types.h"
#include <gsl/gsl_combination.h>
/*=========================================*/
/*define*/
#define MAX_LINE 100000000
#define MAX_STRING 256
#define MAX_NR_LENGTH 15
#define DEFAULT_FILE "house.txt"
/*=========================================*/
/*enumerations*/
typedef enum
{
	file_not_found = 0u,
	file_dim_invalid,
	file_data_invalid, /*not numbers*/
	file_no_error
}T_FILE_ERRORS;
/*=========================================*/
/*typedef*/
typedef struct
{
	uint16 lines;
	uint16 columns;
	double** matrix;

}T_FILE_DIM;
/*=========================================*/
/*external functions*/
boolean 	  open_file(int8 *name_file);
boolean 	  open_file_w(int8* in);
T_FILE_ERRORS fileIsValid(int8 *name_file);
T_FILE_DIM*   get_file_dimensions(void);
void 		  clean_file(void);
void 		  remove_file(const int8* file_name);
/*=========================================*/

#endif
