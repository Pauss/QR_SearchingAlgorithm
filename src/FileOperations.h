#pragma once
#ifndef FILEOPERATIONS_H
#define FILEOPERATIONS_H
/*=========================================*/
/*include*/
#include"Types.h"
/*=========================================*/
/*define*/
#define MAX_LINE 500
#define MAX_STRING 20
#define MAX_NR_LENGTH 15
#define DEFAULT_FILE "house.txt"
/*=========================================*/
/*enumerations*/
typedef enum
{
	FILE_NOT_FOUND = 0u,
	FILE_DIM_INVALID,
	FILE_DATA_INVALID, /*not numbers*/
	FILE_NO_ERROR
}FILE_ERRORS;
/*=========================================*/
/*typedef*/
typedef struct
{
	uint8 lines;
	uint8 columns;
	double** matrix;

}FILE_DIM;
/*=========================================*/
/*external functions*/
boolean open_file(int8 *name_file);
boolean var_f(int8* in);
FILE_ERRORS fileIsValid(int8 *name_file, FILE_DIM* file_dim);
void clean_file(void);
/*=========================================*/

#endif
