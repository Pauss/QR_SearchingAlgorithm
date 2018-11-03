/*=========================================*/
/*
 * MatrixComputations.c
 *
 *  Created on: Nov 2, 2018
 *      Author: Paus
 */
/*=========================================*/
/*include*/
#include <stdio.h>
#include <stdlib.h>
#include "FileOperations.h"
/*=========================================*/
/*private data*/
static FILE *file;
/*=========================================*/
/*global data*/
/*=========================================*/
/*private functions*/
static boolean get_dimension_file(T_FILE_DIM* file_dim);
static boolean get_matrix(T_FILE_DIM* file_dim);

/*=========================================*/
/*This function open a file and rise error if file can't be opened*/
/*=========================================*/
boolean open_file(int8 *name_file) {

	file = fopen(name_file, "r");

	if (file == NULL) {

		return FALSE;
	}

	return TRUE;

}

/*=========================================*/
/*This function is a wrapper for open_file
 *it's purpose is to set a default file if the input string is empty*/
/*=========================================*/
boolean open_file_w(int8* in) {

	/* wrapper to set default file to be opened*/
	int8* i_out = (int8*) malloc(sizeof(int8) * MAX_STRING);

	if (strlen(in) != INIT) {
		strcpy(i_out, in);
	} else {
		strcpy(i_out, DEFAULT_FILE);
	}

	return open_file(i_out);
}

/*=========================================*/
/*special define*/
/*=========================================*/
#define open_file(...) open_file_w((int8*){__VA_ARGS__});

/*=========================================*/
/*This function is getting the size of file: number of columns and lines
 *Also is checking if dimensions are valid*/
/*=========================================*/
static boolean get_dimension_file(T_FILE_DIM* file_dim) {

	int8 *l_line = (int8 *) malloc(sizeof(int8) * MAX_LINE);
	boolean l_inword = FALSE;
	uint8 curent_columns = INIT;

	file_dim->columns = INIT;

	/*set pointer to the start of file*/
	rewind(file);

	while (fgets(l_line, MAX_LINE, file) != NULL) {

		/*get number of lines*/
		file_dim->lines++;
		/*get number of columns of each line*/
		curent_columns = INIT;
		/*get number of columns, from current line*/
		do
			switch (*l_line) {
			case '\0': case ' ':case '\t':case '\n':case '\r':
				/*if delimiter found, number of columns is increased*/
				if (l_inword) {
					l_inword = FALSE;
					curent_columns++;
				}
				break;
			default:
				/*if any other character*/
				l_inword = TRUE;
			}
		/*check each character until reach NULL*/
		while (*l_line++);

		/*set number of columns to be checked*/
		if (INIT == file_dim->columns) {
			file_dim->columns = curent_columns;
		}

		/*if there are 2 consecutive lines with different number of columns than file is corrupt */
		if (file_dim->columns != curent_columns) {
			return FALSE;
		}

	}

	if ((file_dim->columns == INIT) || (file_dim->lines == INIT)) {
		/* one or both dimensions are equal with 0
		 * Any processing is irrelevant for this program purpose*/
		return FALSE;
	}

	return TRUE;
}

/*=========================================*/
/*This function extract data from file into a matrix*/
/*=========================================*/
static boolean get_matrix(T_FILE_DIM* matrix_dim) {

	/*line index*/
	uint8 i = INIT;
	/*column index*/
	uint8 j = INIT;
	char* line1 = (char*) malloc(sizeof(char) * MAX_LINE);

	/*set pointer to the start of file*/
	rewind(file);
	/*ignore first line with description of data*/
	fgets(line1, MAX_LINE, file);

	for (; i < matrix_dim->lines; i++) {
		for (j = INIT; j < matrix_dim->columns; j++) {

			if (!fscanf(file, "%lf", &matrix_dim->matrix[i][j])) {

				/*if data can't be read in double format than is invalid*/
				return FALSE;
			}

		}

	}

	free(line1);
	return TRUE;
}

/*=========================================*/
/*This function checks possible errors of an input file*/
/*=========================================*/
T_FILE_ERRORS fileIsValid(int8* name_file, T_FILE_DIM *file_dim) {

	boolean this_file_checks = FALSE;

	/*initialize only number of columns and lines*/
	file_dim->columns = INIT;
	file_dim->lines = INIT;

	/* check that file is existent on the disk*/
	this_file_checks = open_file(name_file)
	;

	if (FALSE == this_file_checks) {
		return file_not_found;
	}

	/* check that number of columns is equivalent for all lines
	 * get file dimensions: number of lines and number of columns*/
	if (FALSE == get_dimension_file(file_dim)) {
		return file_dim_invalid;
	}

	/* first line is documented with data representation
	 * not used in further computations*/
	file_dim->lines--;

	/* dynamic memory allocation for the matrix
	 * in case there is a big amount of data to compute*/
	file_dim->matrix = malloc((file_dim->lines) * sizeof(double*));
	for (uint8 i = 0; i < file_dim->lines; ++i)
		file_dim->matrix[i] = malloc((file_dim->columns) * sizeof(double));

	/*populate the matrix with data from file*/
	if (FALSE == get_matrix(file_dim)) {
		return file_data_invalid;
	}

	return file_no_error;

}

/*=========================================*/
/*Function to close file*/
/*=========================================*/
void clean_file(void)
{
	fclose(file);
}

